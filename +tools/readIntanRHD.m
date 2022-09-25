function [output, frequency_parameters, notes, spike_triggers, ... 
          varargout] = readIntanRHD(filename, varargin)
% Function to read intan .RHD files. Adapted from read_Intan_RHD2000_file.m
% Version 3.0 as can be downloaded from the manufacuter's website: 
%         http://www.intantech.com/downloads.html
% 
% This function requires user inputs, opens no dialogues and generates no
% command line output (except warnings). 
% 
% Original version (2013) by Laurens Witter
% Updated version (2022) by Calvin Eiber 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% filename : Filename (including path if necessary)
% 'range' : declares that the following input is a range of data (in sec)
%           to load.
% 'ampch' : Lists the amplifier channnels to load
% 'auxch' : Lists the auxilary channels to load
% 'volt'  : Declares whether the voltage channel should be loaded or no
%           ('yes' or 'no')
% 'ADCch' : Lists the onboard ADC channels to load
% 'DIGin' : Lists the onboard digital input channels to load
% 'DIGout' : Lists the onboard digital output channels to load
% 'tempch': Declares whether the temperature data should be loaded ('yes'
%           or 'no')
% 'notch' : Assign '50Hz' or '60Hz' to enable notch filtering. If left
%           undefined then no filter is applied. If the file was saved with
%           notch filtering enabled (v3 or higher), it is not recomputed.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS
% The first output is always the overview of all outputs.
% The second output is always a structure with frequency settings
% The third output is always the notes structure
% The fourth output is always the spike_triggers structure
% Then the function will provide the requested output as data (in
% timeseries format) and channel information pairs. This is always in a fixed
% order:
% 1. amplifier channels
% 2. auxiliary channels
% 3. ADC channels
% 4. Digital channels
% 5. volts
% 6. temperature
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLES:
% [output,freqs,notes,spike_triggers,ampl_data,ampl_ch,aux_data,aux_ch]...
%    = read_Intan(fn,'auxch',[1 2 3],'ampch',[1]);
% This loads the three auxiliary streams, and channel 1 from the general
% amplifier. No range is specified, so the whole stream is loaded.

% #TODO - CE - convert the code to use my named/ 

% First we parse the variable inputs to check whether the list makes sense:
p = inputParser;
p.addRequired('Filename',@(n) ischar(n));
p.addParameter('notch','notchfreq',@(nf) ...
                  strcmp(nf,'50Hz') || strcmp(nf,'60Hz'));
p.addParameter('range','roi',@(rng) length(rng)==2);
p.addParameter('ampch','ch',@(chlist) ~isempty(chlist));
p.addParameter('auxch','aux',@(aux) ~isempty(aux));
p.addParameter('volt','Vch', @(Vch) ...
                       strcmpi(Vch,'yes') || strcmpi(Vch,'no'));
p.addParameter('ADCch','ADC',@(ADC) ~isempty(ADC));
p.addParameter('DIGin','DI',@(DIG) ~isempty(DIG)>0);
p.addParameter('DIGout','DO',@(DIG) ~isempty(DIG)>0);
p.addParameter('tempch','temp',@(temp) ...
                        strcmpi(temp,'yes') || strcmpi(temp,'no'));

p.parse(filename,varargin{:,:})
    
% Now we assign the values defined in the inputs to the respective
% variables and set some defaults.
opts.nrstreamsrequested = 0;
opts.notchfreq = 0;
opts.channellist = 0;
opts.range = 0;
opts.do_aux = 0;
opts.do_supply_volts = 0;
opts.do_ADC_channels = 0;
opts.do_DI_channels = 0;
opts.do_DO_channels = 0;
opts.do_temperature = 0;
for arg=1:2:length(varargin)
    switch (varargin{arg})
        case 'notch'
            opts.notchfreq = varargin{arg+1};
        case 'ampch'
            opts.channellist = varargin{arg+1};
            opts.nrstreamsrequested=opts.nrstreamsrequested+1;
        case 'range'
            opts.range = varargin{arg+1};
        case 'auxch'
            opts.do_aux = varargin{arg+1};
            opts.nrstreamsrequested=opts.nrstreamsrequested+1;
        case 'volt'
            opts.do_supply_volts = 1;
            opts.nrstreamsrequested=opts.nrstreamsrequested+1;
        case 'ADCch'
            opts.do_ADC_channels = varargin{arg+1};
            opts.nrstreamsrequested=opts.nrstreamsrequested+1;
        case 'DIGin'
            opts.do_DI_channels = varargin{arg+1};
            opts.nrstreamsrequested=opts.nrstreamsrequested+1;
        case 'DIGout'
            opts.do_DO_channels = varargin{arg+1};
            opts.nrstreamsrequested=opts.nrstreamsrequested+1;            
        case 'tempch'
            opts.do_temperature = 1;
            opts.nrstreamsrequested=opts.nrstreamsrequested+1;
    end
end
    
% The actual start of the file opening. Here we first open the file to gain
% access to the parts in it.
fid = fopen(filename, 'r');
s = dir(filename);
filesize = s.bytes;
    
when_done_fclose = onCleanup(@() fclose(fid)); 

% Check 'magic number' at beginning of file to make sure this is an Intan
% Technologies RHD2000 data file.
magic_number = fread(fid, 1, 'uint32');
if magic_number ~= hex2dec('c6912702')
    error('Unrecognized file type.');
end

% Read version number.
data_file_main_version_number = fread(fid, 1, 'int16');
data_file_secondary_version_number = fread(fid, 1, 'int16');

fprintf('Reading Intan Technologies RHD2000 Data File, Version %d.%d\n', ...
         data_file_main_version_number, data_file_secondary_version_number);
printInfo(); % Added CE 12 Sep 2022


if (data_file_main_version_number == 1), num_samples_per_data_block = 60;
else                                     num_samples_per_data_block = 128;
end

% Read information of sampling rate and amplifier frequency settings.
sample_rate = fread(fid, 1, 'single');
dsp_enabled = fread(fid, 1, 'int16');
actual_dsp_cutoff_frequency = fread(fid, 1, 'single');
actual_lower_bandwidth = fread(fid, 1, 'single');
actual_upper_bandwidth = fread(fid, 1, 'single');

desired_dsp_cutoff_frequency = fread(fid, 1, 'single');
desired_lower_bandwidth = fread(fid, 1, 'single');
desired_upper_bandwidth = fread(fid, 1, 'single');

% This tells us if a software 50/60 Hz notch filter was enabled during
% the data acquisition.
notch_filter_mode = fread(fid, 1, 'int16');
notch_filter_frequency = 0;
if (notch_filter_mode == 1),     notch_filter_frequency = 50;
elseif (notch_filter_mode == 2), notch_filter_frequency = 60;
end

desired_impedance_test_frequency = fread(fid, 1, 'single');
actual_impedance_test_frequency = fread(fid, 1, 'single');
    
% Place notes in data strucure
notes = struct( 'note1', fread_QString(fid), ...
                'note2', fread_QString(fid), ...
                'note3', fread_QString(fid) );
    
% If data file is from GUI v1.1 or later, 
% see if temperature sensor data was saved.
num_temp_sensor_channels = 0;
if ((data_file_main_version_number == 1 && data_file_secondary_version_number >= 1) ...
    || (data_file_main_version_number > 1))
    num_temp_sensor_channels = fread(fid, 1, 'int16');
end

% If data file is from GUI v1.3 or later, load board mode.
board_mode = 0;
if ((data_file_main_version_number == 1 && data_file_secondary_version_number >= 3) ...
    || (data_file_main_version_number > 1))
    board_mode = fread(fid, 1, 'int16');
end

% If data file is from v2.0 or later (Intan Recording Controller),
% load name of digital reference channel.
if (data_file_main_version_number > 1)
    notes.reference_channel = fread_QString(fid);
end


% Place frequency-related information in data structure.
frequency_parameters = struct( ...
        'amplifier_sample_rate', sample_rate, ...
        'aux_input_sample_rate', sample_rate / 4, ...
        'supply_voltage_sample_rate', sample_rate / 60, ...
        'board_adc_sample_rate', sample_rate, ...
        'board_dig_in_sample_rate', sample_rate, ...
        'desired_dsp_cutoff_frequency', desired_dsp_cutoff_frequency, ...
        'actual_dsp_cutoff_frequency', actual_dsp_cutoff_frequency, ...
        'dsp_enabled', dsp_enabled, ...
        'desired_lower_bandwidth', desired_lower_bandwidth, ...
        'actual_lower_bandwidth', actual_lower_bandwidth, ...
        'desired_upper_bandwidth', desired_upper_bandwidth, ...
        'actual_upper_bandwidth', actual_upper_bandwidth, ...
        'notch_filter_frequency', notch_filter_frequency, ...
        'desired_impedance_test_frequency', desired_impedance_test_frequency, ...
        'actual_impedance_test_frequency', actual_impedance_test_frequency );
    
% Define data structure for spike trigger settings.
spike_trigger_struct = struct( ...
        'voltage_trigger_mode', {}, ...
        'voltage_threshold', {}, ...
        'digital_trigger_channel', {}, ...
        'digital_edge_polarity', {} );

new_trigger_channel = struct(spike_trigger_struct);
spike_triggers = struct(spike_trigger_struct);

% Define data structure for data channels.
channel_struct = struct( ...
        'native_channel_name', {}, ...
        'custom_channel_name', {}, ...
        'native_order', {}, ...
        'custom_order', {}, ...
        'board_stream', {}, ...
        'chip_channel', {}, ...
        'port_name', {}, ...
        'port_prefix', {}, ...
        'port_number', {}, ...
        'electrode_impedance_magnitude', {}, ...
        'electrode_impedance_phase', {} );

new_channel = struct(channel_struct);

% Create structure arrays for each type of data channel.
amplifier_channels = struct(channel_struct);
aux_input_channels = struct(channel_struct);
supply_voltage_channels = struct(channel_struct);
board_adc_channels = struct(channel_struct);
board_dig_in_channels = struct(channel_struct);
board_dig_out_channels = struct(channel_struct);

amplifier_index = 1;
aux_input_index = 1;
supply_voltage_index = 1;
board_adc_index = 1;
board_dig_in_index = 1;
board_dig_out_index = 1;
    
% Read signal summary from data file header.
number_of_signal_groups = fread(fid, 1, 'int16');

for signal_group = 1:number_of_signal_groups

    signal_group_name = fread_QString(fid);
    signal_group_prefix = fread_QString(fid);
    signal_group_enabled = fread(fid, 1, 'int16');
    signal_group_num_channels = fread(fid, 1, 'int16');
    signal_group_num_amp_channels = fread(fid, 1, 'int16');

    % printInfo added CE 12 Sep 2022
    printInfo('(%d/%d) Reading Signal Group %s [%d chanels] ', ...
               signal_group, number_of_signal_groups, ...
               signal_group_name, signal_group_num_amp_channels)

    if (signal_group_num_channels <= 0 || signal_group_enabled <= 0)
        continue
    end

    new_channel(1).port_name = signal_group_name;
    new_channel(1).port_prefix = signal_group_prefix;
    new_channel(1).port_number = signal_group;

    for signal_channel = 1:signal_group_num_channels

        new_channel(1).native_channel_name = fread_QString(fid);
        new_channel(1).custom_channel_name = fread_QString(fid);
        new_channel(1).native_order = fread(fid, 1, 'int16');
        new_channel(1).custom_order = fread(fid, 1, 'int16');
        signal_type = fread(fid, 1, 'int16');
        channel_enabled = fread(fid, 1, 'int16');
        new_channel(1).chip_channel = fread(fid, 1, 'int16');
        new_channel(1).board_stream = fread(fid, 1, 'int16');
        new_trigger_channel(1).voltage_trigger_mode = fread(fid, 1, 'int16');
        new_trigger_channel(1).voltage_threshold = fread(fid, 1, 'int16');
        new_trigger_channel(1).digital_trigger_channel = fread(fid, 1, 'int16');
        new_trigger_channel(1).digital_edge_polarity = fread(fid, 1, 'int16');
        new_channel(1).electrode_impedance_magnitude = fread(fid, 1, 'single');
        new_channel(1).electrode_impedance_phase = fread(fid, 1, 'single');

        if ~channel_enabled, continue, end

        switch signal_type
          case 0
            amplifier_channels(amplifier_index) = new_channel;
            spike_triggers(amplifier_index) = new_trigger_channel;
            amplifier_index = amplifier_index + 1;
          case 1
            aux_input_channels(aux_input_index) = new_channel;
            aux_input_index = aux_input_index + 1;
          case 2
            supply_voltage_channels(supply_voltage_index) = new_channel;
            supply_voltage_index = supply_voltage_index + 1;
          case 3
            board_adc_channels(board_adc_index) = new_channel;
            board_adc_index = board_adc_index + 1;
          case 4
            board_dig_in_channels(board_dig_in_index) = new_channel;
            board_dig_in_index = board_dig_in_index + 1;
          case 5
            board_dig_out_channels(board_dig_out_index) = new_channel;
            board_dig_out_index = board_dig_out_index + 1;
          otherwise
            error('Unknown channel type');
        end
    end
end

fprintf('\n'); % terminate printInfo ... 

% Summarize contents of data file.
num_amplifier_channels = amplifier_index - 1;
num_aux_input_channels = aux_input_index - 1;
num_supply_voltage_channels = supply_voltage_index - 1;
num_board_adc_channels = board_adc_index - 1;
num_board_dig_in_channels = board_dig_in_index - 1;
num_board_dig_out_channels = board_dig_out_index - 1;

% Determine how many samples the data file contains.

% Each data block contains num_samples_per_data_block amplifier samples.
bytes_per_block = num_samples_per_data_block * 4;  % timestamp data
bytes_per_block = bytes_per_block + num_samples_per_data_block * 2 * num_amplifier_channels;
% Auxiliary inputs are sampled 4x slower than amplifiers
bytes_per_block = bytes_per_block + (num_samples_per_data_block / 4) * 2 * num_aux_input_channels;
% Supply voltage is sampled once per data block
bytes_per_block = bytes_per_block + 1 * 2 * num_supply_voltage_channels;
% Board analog inputs are sampled at same rate as amplifiers
bytes_per_block = bytes_per_block + num_samples_per_data_block * 2 * num_board_adc_channels;
% Board digital inputs are sampled at same rate as amplifiers
if (num_board_dig_in_channels > 0)
    bytes_per_block = bytes_per_block + num_samples_per_data_block * 2;
end
% Board digital outputs are sampled at same rate as amplifiers
if (num_board_dig_out_channels > 0)
    bytes_per_block = bytes_per_block + num_samples_per_data_block * 2;
end
% Temp sensor is sampled once per data block
if (num_temp_sensor_channels > 0)
   bytes_per_block = bytes_per_block + 1 * 2 * num_temp_sensor_channels; 
end

% How many data blocks remain in this file?
bytes_remaining = filesize - ftell(fid);
data_present = (bytes_remaining > 0);
num_data_blocks = bytes_remaining / bytes_per_block;

num_amplifier_samples = num_samples_per_data_block * num_data_blocks;
num_aux_input_samples = (num_samples_per_data_block / 4) * num_data_blocks;
num_supply_voltage_samples = 1 * num_data_blocks;
num_board_adc_samples = num_samples_per_data_block * num_data_blocks;
num_board_dig_in_samples = num_samples_per_data_block * num_data_blocks;
num_board_dig_out_samples = num_samples_per_data_block * num_data_blocks;

% record_time = num_amplifier_samples / sample_rate;

% Pre-allocate memory for data.
t_amplifier = zeros(1, num_amplifier_samples);

assert(data_present,'this file does not contain data to load!')

amplifier_data = zeros(num_amplifier_channels, num_amplifier_samples);
aux_input_data = zeros(num_aux_input_channels, num_aux_input_samples);
supply_voltage_data = zeros(num_supply_voltage_channels, num_supply_voltage_samples);
temp_sensor_data = zeros(num_temp_sensor_channels, num_supply_voltage_samples);
board_adc_data = zeros(num_board_adc_channels, num_board_adc_samples);
board_dig_in_data = zeros(num_board_dig_in_channels, num_board_dig_in_samples);
board_dig_in_raw = zeros(1, num_board_dig_in_samples);
board_dig_out_data = zeros(num_board_dig_out_channels, num_board_dig_out_samples);
board_dig_out_raw = zeros(1, num_board_dig_out_samples);

% Read sampled data from file.

amplifier_index = 1;
aux_input_index = 1;
supply_voltage_index = 1;
board_adc_index = 1;
board_dig_in_index = 1;

for i=1:num_data_blocks
    % In version 1.2, we moved from saving timestamps as unsigned
    % integeters to signed integers to accomidate negative (adjusted)
    % timestamps for pretrigger data.
    if (data_file_main_version_number > 1) || ...
       (data_file_main_version_number == 1 && ...
        data_file_secondary_version_number >= 2)
         t_amplifier(amplifier_index:(amplifier_index + num_samples_per_data_block - 1)) = fread(fid, num_samples_per_data_block, 'int32');
    else t_amplifier(amplifier_index:(amplifier_index + num_samples_per_data_block - 1)) = fread(fid, num_samples_per_data_block, 'uint32');
    end
    if (num_amplifier_channels > 0)
        amplifier_data(:, amplifier_index:(amplifier_index + num_samples_per_data_block - 1)) = fread(fid, [num_samples_per_data_block, num_amplifier_channels], 'uint16')';
    end
    if (num_aux_input_channels > 0)
        aux_input_data(:, aux_input_index:(aux_input_index + (num_samples_per_data_block / 4) - 1)) = fread(fid, [(num_samples_per_data_block / 4), num_aux_input_channels], 'uint16')';
    end
    if (num_supply_voltage_channels > 0)
        supply_voltage_data(:, supply_voltage_index) = fread(fid, [1, num_supply_voltage_channels], 'uint16')';
    end
    if (num_temp_sensor_channels > 0)
        temp_sensor_data(:, supply_voltage_index) = fread(fid, [1, num_temp_sensor_channels], 'int16')';
    end
    if (num_board_adc_channels > 0)
        board_adc_data(:, board_adc_index:(board_adc_index + num_samples_per_data_block - 1)) = fread(fid, [num_samples_per_data_block, num_board_adc_channels], 'uint16')';
    end
    if (num_board_dig_in_channels > 0)
        board_dig_in_raw(board_dig_in_index:(board_dig_in_index + num_samples_per_data_block - 1)) = fread(fid, num_samples_per_data_block, 'uint16');
    end
    if (num_board_dig_out_channels > 0)
        board_dig_out_raw(board_dig_out_index:(board_dig_out_index + num_samples_per_data_block - 1)) = fread(fid, num_samples_per_data_block, 'uint16');
    end

    amplifier_index = amplifier_index + num_samples_per_data_block;
    aux_input_index = aux_input_index + (num_samples_per_data_block / 4);
    supply_voltage_index = supply_voltage_index + 1;
    board_adc_index = board_adc_index + num_samples_per_data_block;
    board_dig_in_index = board_dig_in_index + num_samples_per_data_block;
    board_dig_out_index = board_dig_out_index + num_samples_per_data_block;

end


% Make sure we have read exactly the right amount of data.
bytes_remaining = filesize - ftell(fid);
if (bytes_remaining ~= 0)
    warning('End of file not reached.');
end

% Close data file.
clear when_done_fclose % fclose(fid);


% Extract digital input channels to separate variables.
for i=1:num_board_dig_in_channels
    mask = 2^(board_dig_in_channels(i).native_order) * ones(size(board_dig_in_raw));
    board_dig_in_data(i, :) = (bitand(board_dig_in_raw, mask) > 0);
end
for i=1:num_board_dig_out_channels
    mask = 2^(board_dig_out_channels(i).native_order) * ones(size(board_dig_out_raw));
    board_dig_out_data(i, :) = (bitand(board_dig_out_raw, mask) > 0);
end

% Scale voltage levels appropriately.
amplifier_data = 0.195 * (amplifier_data - 32768); % units = microvolts
aux_input_data = 37.4e-6 * aux_input_data; % units = volts
supply_voltage_data = 74.8e-6 * supply_voltage_data; % units = volts
temp_sensor_data = temp_sensor_data / 100; % units = deg C

if (board_mode == 1)
     board_adc_data = 152.59e-6 * (board_adc_data - 32768); % units = volts
elseif (board_mode == 13) % Intan Recording Controller
     board_adc_data = 312.5e-6 * (board_adc_data - 32768); % units = volts    
else board_adc_data = 50.354e-6 * board_adc_data; % units = volts
end

% Check for gaps in timestamps.
num_gaps = sum(diff(t_amplifier) ~= 1);
if (num_gaps == 0)
else fprintf('Warning: %d gaps in timestamp data found. %s\n', ...
                       num_gaps, 'Time scale will not be uniform!');
end

% Scale time steps (units = seconds).
t_amplifier = t_amplifier / sample_rate;
t_aux_input = t_amplifier(1:4:end);
t_supply_voltage = t_amplifier(1:60:end);
t_board_adc = t_amplifier;
t_dig = t_amplifier;
t_temp_sensor = t_supply_voltage;

% If the software notch filter was selected during the recording, apply the
% same notch filter to amplifier data here.  But don't do this for v3.0+ 
% files (from Intan RHX software) because RHX saves notch-filtered data.
if (notch_filter_frequency > 0 && data_file_main_version_number < 3)
  for i=1:num_amplifier_channels
      amplifier_data(i,:) = notch_filter(amplifier_data(i,:), ...
                                         sample_rate, ...
                                         notch_filter_frequency, 10);
  end
end


% Now that we have the whole file we select the parts we need.
% First we select the datastreams requested, then select the data range
% requested. With every step we check whether we do have the data in
% the file actually. If something is requested, but not included in the
% file, a warning is displayed in the commandline.
    
if opts.nrstreamsrequested==0  % Nothing specified, all output requested
    % Output order:
    % amplifier channels => auxiliary channels => ADC channels => Digital
    % channels => volts => temperature
    output = {'ampch', 1:size(amplifier_data,1), ...
              'auxch', 1:size(aux_input_data,1), ...
              'ADCch', 1:size(board_adc_data,1), ...
              'DIGin', 1:size(board_dig_in_data,1), ...
              'DIGout', 1:size(board_dig_out_data,1), ...
              'volt', 'yes', ...
              'tempch', 'yes' };
else
    output = cell(0);
    arg_named = @(x) ismember(varargin(1:2:end),x);
    if any(arg_named('ampch'))
        output(1,end+(1:2)) = {'ampch', opts.channellist};
    end
    if any(arg_named('auxch'))
        output(1,end+(1:2)) = {'auxch', opts.do_aux};
    end
    if any(arg_named('ADCch'))
        output(1,end+(1:2)) = {'ADCch', opts.do_ADC_channels};
    end
    if any(arg_named('DIGin'))
        output(1,end+(1:2)) = {'DIGin', opts.do_DI_channels};
    end
    if any(arg_named('DIGout'))
        output(1,end+(1:2)) = {'DIGout', opts.do_DO_channels};
    end
    if any(arg_named('volt'))
        output(1,end+(1:2)) = {'volt', opts.do_supply_volts};
    end
    if any(arg_named('tempch'))
        output(1,end+(1:2)) = {'tempch', opts.do_supply_volts};
    end
end

for arg=1:2:length(output)
  switch (output{arg})
    case 'ampch'
      if size(amplifier_data,1) > 0
        amplifier_channels = amplifier_channels(1,output{arg+1});
        amplifier_data = amplifier_data(output{arg+1},:);
        ts_amp = timeseries;
        ts_amp.name = 'ampch';
        ts_amp.time = t_amplifier;
        ts_amp.data = amplifier_data';
        if opts.range ~=0
            ts_amp = select_time(ts_amp,opts.range);
        end
        ts_amp.DataInfo.Units = 'uV';
        ts_amp.TimeInfo.Units = 'seconds';
        varargout{arg} = ts_amp;
        varargout{arg+1} = amplifier_channels;
      elseif opts.nrstreamsrequested > 0
        disp('WARNING: No amplifier channels found in file!')
      end
    case 'auxch'
      if size(aux_input_data,1) > 0
        aux_input_channels = aux_input_channels(1,output{arg+1});
        aux_input_data = aux_input_data(output{arg+1},:);
        ts_aux = timeseries;
        ts_aux.name = 'auxch';
        ts_aux.time = t_aux_input;
        ts_aux.data = aux_input_data';
        if opts.range~=0
            ts_aux = select_time(ts_aux,opts.range);
        end
        ts_aux.DataInfo.Unit = 'uV';
        ts_aux.TimeInfo.Units = 'seconds';
        varargout{arg} = ts_aux;
        varargout{arg+1} = aux_input_channels;
      elseif opts.nrstreamsrequested > 0
        disp('WARNING: No auxiliary inputs found in file!')
      end
    case 'ADCch'
      if size(board_adc_data,1)>0
        board_adc_channels = board_adc_channels(1,output{arg+1});
        board_adc_data = board_adc_data(output{arg+1},:);
        ts_board_adc = timeseries;
        ts_board_adc.name = 'ADCch';
        ts_board_adc.time = t_board_adc;
        ts_board_adc.data = board_adc_data';
        if opts.range ~=0
            ts_board_adc = select_time(ts_board_adc,opts.range);
        end
        ts_board_adc.DataInfo.Unit = 'uV';
        ts_board_adc.TimeInfo.Units = 'seconds';
        varargout{arg} = ts_board_adc;
        varargout{arg+1} = board_adc_channels;
      elseif opts.nrstreamsrequested > 0
        disp('WARNING: No board ADC channels found in file!');
      end
    case 'DIGin'
      if size(board_dig_in_data,1)>0
        board_dig_in_channels = board_dig_in_channels(1,output{arg+1});
        board_dig_in_data = board_dig_in_data(output{arg+1},:);
        ts_dig_in = timeseries;
        ts_dig_in.name = 'DIGin';
        ts_dig_in.time = t_amplifier;
        ts_dig_in.data = board_dig_in_data';
        if opts.range ~=0
            ts_dig_in = select_time(ts_dig_in,opts.range);
        end
        ts_dig_in.TimeInfo.Units = 'seconds';
        varargout{arg} = ts_dig_in;
        varargout{arg+1} = board_dig_in_channels;
      elseif opts.nrstreamsrequested > 0
        disp('WARNING: No board digital in channels found in file!');
      end
    case 'DIGout'
      if size(board_dig_out_data,1)>0
        board_dig_out_channels = board_dig_out_channels(1,output{arg+1});
        board_dig_out_data = board_dig_out_data(output{arg+1},:);
        ts_dig_out = timeseries;
        ts_dig_out.name = 'DIGout';
        ts_dig_out.time = t_amplifier;
        ts_dig_out.data = board_dig_out_data';
        if opts.range ~=0
            ts_dig_out = select_time(ts_dig_out,opts.range);
        end
        ts_dig_out.TimeInfo.Units = 'seconds';
        varargout{arg} = ts_dig_out;
        varargout{arg+1} = board_dig_out_channels;
      elseif opts.nrstreamsrequested > 0
        disp('WARNING: No board digital out channels found in file!');
      end
    case 'volt'
      if size(supply_voltage_data,1)>0
        ts_volt = timeseries;
        ts_volt.name = 'supply_voltage';
        ts_volt.time = t_supply_voltage;
        ts_volt.data = supply_voltage_data';
        if opts.range ~=0
            ts_volt = select_time(ts_volt,opts.range);
        end
        ts_volt.TimeInfo.Units = 'seconds';
        varargout{arg} = ts_volt;
        varargout{arg+1} = supply_voltage_channels;
      elseif opts.nrstreamsrequested > 0
        disp('WARNING: No supply voltage data found in file!');
      end
    case 'tempch'
      if size(temp_sensor_data,1)>0
        ts_temp = timeseries;
        ts_temp.name = 'temperature';
        ts_temp.time = t_temp_sensor;
        ts_temp.data = temp_sensor_data';
        if opts.range ~= 0
            ts_temp = select_time(ts_temp,opts.range);
        end
        ts_temp.TimeInfo.Units = 'seconds';
        varargout{arg} = ts_temp;
        varargout{arg+1} = '1';
      else
        disp('WARNING: No temperature data found in file!');
      end
  end
end
 
return


function new = select_time(old,time_range)
% Enter a timeseries and a range.
% Returns a new timeseries according to range.
% Range should be a 1x2 matrix with [start end]
new = timeseries;
data = old.data(old.time>time_range(1,1),:);
time = old.time(old.time>time_range(1,1));
data = data(time<time_range(1,2),:);
time = time(time<time_range(1,2));
new.data = data;
new.time = time;

return

function a = fread_QString(fid)

% a = read_QString(fid)
%
% Read Qt style QString.  The first 32-bit unsigned number indicates
% the length of the string (in 16-bit unicode characters).  If this
% number equals 0xFFFFFFFF, the string is null.

a = '';
length = fread(fid, 1, 'uint32');
if length == hex2num('ffffffff')
    return;
end

% convert length from bytes to 16-bit unicode words
length = length / 2;
a = char(fread(fid, length, 'uint16'))'; %#ok<FREAD> 

% for i=1:length
%     a(i) = fread(fid, 1, 'uint16');
% end

return


function out = notch_filter(in, fSample, fNotch, Bandwidth)
% out = notch_filter(in, fSample, fNotch, Bandwidth)
%
% Implements a notch filter (e.g., for 50 or 60 Hz) on vector 'in'.
% fSample = sample rate of data (in Hz or Samples/sec)
% fNotch = filter notch frequency (in Hz)
% Bandwidth = notch 3-dB bandwidth (in Hz).  A bandwidth of 10 Hz is
%   recommended for 50 or 60 Hz notch filters; narrower bandwidths lead to
%   poor time-domain properties with an extended ringing response to
%   transient disturbances.
%
% Example:  If neural data was sampled at 30 kSamples/sec
% and you wish to implement a 60 Hz notch filter:
%
% out = notch_filter(in, 30000, 60, 10);

tstep = 1/fSample;
Fc = fNotch*tstep;

L = length(in);

% Calculate IIR filter parameters
d = exp(-2*pi*(Bandwidth/2)*tstep);
b = (1 + d*d)*cos(2*pi*Fc);

a0 = 1; a1 = -b; a2 = d*d;
a = (1 + d*d)/2;
b0 = 1; b1 = -2*cos(2*pi*Fc); b2 = 1;

out = zeros(size(in));
out(1) = in(1);  
out(2) = in(2);

% (If filtering a continuous data stream, change out(1) and out(2) to the
%  previous final two values of out.)

% Run filter
for i=3:L
    out(i) = (a*b2*in(i-2) + a*b1*in(i-1) + a*b0*in(i) - a2*out(i-2) - a1*out(i-1))/a0;
end
return





% Compressed text output function added CE 12 Sep 2022
function printInfo(varargin)

persistent s;
if nargin == 0, s = ''; return, end
fprintf('%s',char(8*ones(size(s)))); 
s = sprintf(varargin{:});
fprintf('%s',s)

return