function output = readIntanRHD(filename, varargin)
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
% 'roi'   : Return a specified range of data (in sec)
% 'amp'   : List the amplifier channnels to load
% 'aux'   : List the auxilary channels to load
% 'volt'  : Should the supply voltage channel should be loaded? (logical)
% 'ADC'   : List the on-board ADC channels to load
% 'DI'    : List the onboard digital input channels to load
% 'DO'    : List the onboard digital output channels to load
% 'temp'  : Should temperature data should be loaded? (logical)
% 'notch' : Assign '50Hz' or '60Hz' to enable notch filtering. If left
%           undefined then no filter is applied. If the file was saved with
%           notch filtering enabled (v3 or higher), it is not recomputed.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS
% 
% The output is a structure containing the following fields: 
% 
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

named = @(n) strncmpi(varargin,n,length(n));
get_ = @(v) varargin{find(named(v))+1};

if nargin < 1, filename = '?'; end
if ~exist(filename,'file')
    [fn,fp] = uigetfile('*.rhd',[],filename);
    if ~any(fn), error('Cancelled'), end
    filename = [fp fn];
end

% Now we assign the values defined in the inputs to the respective
% variables and set some defaults.
opts.n_streams_req = 0;
opts.notchfreq = 0;
opts.time_roi = 0;

% opts.do_AMP_channels = 0;
% opts.do_AUX_channels = 0;
% opts.do_VOLT_channels = 0;
% opts.do_ADC_channels = 0;
% opts.do_DI_channels = 0;
% opts.do_DO_channels = 0;
% opts.do_TEMP_channels = 0;

if any(named('notch')), opts.notchfreq = get_('notch'); end
if any(named('roi')),       opts.time_roi = get_('roi');
elseif any(named('range')), opts.time_roi = get_('range');
end

channel_types = {'amp','aux','volt','adc','DI','DO','temp'};

for chan_type = channel_types

  opt_id = ['do_' upper(chan_type{1}) '_channels'];
  opts.(opt_id) = false; 

  if any(named(chan_type{1})), opts.(opt_id) = get_('amp');  
          opts.n_streams_req = opts.n_streams_req+1;
  end
end
    
% The actual start of the file opening. Here we first open the file to gain
% access to the parts in it.
fid = fopen(filename, 'r');
s = dir(filename);
filesize = s.bytes;

output.filename = filename; 
output.notes = ''; 
output.config = struct;

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
output.notes = struct( 'note1', fread_QString(fid), ...
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
    output.notes.reference_channel = fread_QString(fid);
end

output.notes.version = [data_file_main_version_number ... 
                        data_file_secondary_version_number];

% Place frequency-related information in data structure.
output.config.general = struct( ...
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

channels = struct;
chan_ = @(ty) [upper(ty{1}) '_channels'];
sams_ = @(ty) [upper(ty{1}) '_samples'];

for ty = channel_types
    channels.(chan_(ty)) = struct(channel_struct);
end

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
          case 0, channels.AMP_channels(end+1) = new_channel;
            spike_triggers(end+1) = new_trigger_channel;
          case 1, channels.AUX_channels(end+1) = new_channel;
          case 2, channels.VOLT_channels(end+1) = new_channel;
          case 3, channels.ADC_channels(end+1) = new_channel;
          case 4, channels.DI_channels(end+1) = new_channel;
          case 5, channels.DO_channels(end+1) = new_channel;
          otherwise
            error('Unknown channel type');
        end
    end
end

fprintf('\n'); % terminate printInfo ... 

% Summarize contents of data file.
num = struct;
for ty = channel_types
    num.(chan_(ty)) = numel(channels.(chan_(ty)));
end

% Determine how many samples the data file contains.
nspd = num_samples_per_data_block;

% Each data block contains num_samples_per_data_block amplifier samples.
bytes_per_block = nspd * 4;  % timestamp data
bytes_per_block = bytes_per_block + nspd * 2 * num.AMP_channels;
% Auxiliary inputs are sampled 4x slower than amplifiers
bytes_per_block = bytes_per_block + (nspd / 4) * 2 * num.AUX_channels;
% Supply voltage is sampled once per data block
bytes_per_block = bytes_per_block + 1 * 2 * num.VOLT_channels;
% Board analog inputs are sampled at same rate as amplifiers
bytes_per_block = bytes_per_block + nspd * 2 * num.ADC_channels;
% Board digital inputs are sampled at same rate as amplifiers
if (num.DI_channels > 0)
    bytes_per_block = bytes_per_block + nspd * 2;
end
% Board digital outputs are sampled at same rate as amplifiers
if (num.DO_channels > 0)
    bytes_per_block = bytes_per_block + nspd * 2;
end
% Temp sensor is sampled once per data block
bytes_per_block = bytes_per_block + 1 * 2 * num.TEMP_channels; 

% How many data blocks remain in this file?
bytes_remaining = filesize - ftell(fid);
data_present = (bytes_remaining > 0);
num_data_blocks = bytes_remaining / bytes_per_block;

num.AMP_samples  = nspd * num_data_blocks;
num.AUX_samples  = (nspd / 4) * num_data_blocks;
num.VOLT_samples = 1 * num_data_blocks;
num.ADC_samples  = nspd * num_data_blocks;
num.DI_samples   = nspd * num_data_blocks;
num.DO_samples   = nspd * num_data_blocks;
num.TEMP_samples = num_temp_sensor_channels;

% record_time = num_amplifier_samples / sample_rate;

% Pre-allocate memory for data.
t_amplifier = zeros(1, num.AMP_samples);

assert(data_present,'this file does not contain data to load!')

data = struct; 
for ty = channel_types
    data.(chan_(ty)) = zeros(num.(chan_(ty)), num.(sams_(ty)));
    data.(sams_(ty)) = 1; 
end

data.DI_raw = zeros(1,num.DI_channels);
data.DO_raw = zeros(1,num.DO_channels);

% Read sampled data from file.

for ii = 1:num_data_blocks
    % In version 1.2, we moved from saving timestamps as unsigned
    % integeters to signed integers to accomidate negative (adjusted)
    % timestamps for pretrigger data.
    if (data_file_main_version_number > 1) || ...
       (data_file_main_version_number == 1 && ...
        data_file_secondary_version_number >= 2)
         t_amplifier(amplifier_index:(amplifier_index + nspd - 1)) = fread(fid, nspd, 'int32');
    else t_amplifier(amplifier_index:(amplifier_index + nspd - 1)) = fread(fid, nspd, 'uint32');
    end
    

    if (num.AMP_channels > 0)
      data.AMP_chanels(:, data.AMP_samples + (0 : nspd-1)) = ...
        fread(fid, [nspd, num.AMP_channels], 'uint16')';
    end
    if (num.AUX_channels > 0)
      data.AUX_channels(:, data.AUX_samples + (0:(nspd/4)-1)) = ... 
        fread(fid, [(nspd / 4), num.AUX_channels], 'uint16')';
    end
    if (num.VOLT_channels > 0)
      data.VOLT_channels(:, data.VOLT_samples) = ...
        fread(fid, [1, num.VOLT_channels], 'uint16')';
    end
    if (num.TEMP_channels > 0)
      data.TEMP_channels(:, data.TEMP_samples) = ...
        fread(fid, [1, num.TEMP_channels], 'int16')';
    end
    if (num.ADC_channels > 0)
      data.ADC_channels(:, data.ADC_samples + (0 : nspd-1)) = ...
        fread(fid, [nspd, num.ADC_channels], 'uint16')';
    end
    if (num.DI_channels > 0)
      data.DI_raw(data.DI_samples + (0 : nspd-1)) = ...
        fread(fid, nspd, 'uint16');
    end
    if (num.DO_channels > 0)
      data.DO_raw(data.DO_samples + (0: nspd-1)) = ...
        fread(fid, nspd, 'uint16');
    end
    
    % update indices 
    data.AMP_samples = data.AMP_samples + nspd;
    data.AUX_samples = data.AUX_samples + (nspd / 4);
    data.VOLT_samples = data.VOLT_samples + 1;
    data.TEMP_samples = data.TEMP_samples + 1;
    data.ADC_samples = data.ADC_samples + nspd;
    data.DI_samples = data.DI_samples + nspd;
    data.DO_samples = data.DO_samples + nspd;
end


% Make sure we have read exactly the right amount of data.
bytes_remaining = filesize - ftell(fid);
if (bytes_remaining ~= 0)
    warning('End of file not reached.');
end

% Close data file.
clear when_done_fclose % fclose(fid);

% Extract digital input channels to separate variables.
for ii = 1:num.DI_channels
    mask = 2^(channels.DI_channels(ii).native_order) * ones(size(data.DI_raw));
    data.DI_channels(ii, :) = (bitand(data.DI_raw, mask) > 0);
end
for ii=1:num.DO_channels
    mask = 2^(channels.DO_channels(ii).native_order) * ones(size(data.DO_raw));
    data.DO_channels(ii, :) = (bitand(data.DO_raw, mask) > 0);
end

% Scale voltage levels appropriately.
data.AMP_channels = 0.195 * (data.AMP_channels - 32768); % units = microvolts
data.AUX_channels = 37.4e-6 * data.AUX_channels; % units = volts
data.VOLT_channels = 74.8e-6 * data.VOLT_channels; % units = volts
data.TEMP_channels = data.TEMP_channels / 100; % units = deg C

if (board_mode == 1)
     data.ADC_channels = 152.59e-6 * (data.ADC_channels - 32768); % units = volts
elseif (board_mode == 13) % Intan Recording Controller
     data.ADC_channels = 312.5e-6 * (data.ADC_channels - 32768); % units = volts    
else data.ADC_channels = 50.354e-6 * data.ADC_channels; % units = volts
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
  for ii = 1:num.AMP_channels
      data.AMP_channels(ii,:) = notch_filter(data.AMP_channels(ii,:), ...
                                             sample_rate, ...
                                             notch_filter_frequency, 10);
  end
end

% Now that we have the whole file we select the parts we need.
% First we select the datastreams requested, then select the data range
% requested. With every step we check whether we do have the data in
% the file actually. If something is requested, but not included in the
% file, a warning is displayed in the commandline.
    
error I_think_redundant

if opts.n_streams_req == 0  % Nothing specified, all output requested
    % Output order:
    % amplifier channels => auxiliary channels => ADC channels => Digital
    % channels => volts => temperature
    selected_output = {'ampch', 1:size(amplifier_data,1), ...
                       'auxch', 1:size(aux_input_data,1), ...
                       'ADCch', 1:size(board_adc_data,1), ...
                       'DIGin', 1:size(board_dig_in_data,1), ...
                       'DIGout', 1:size(board_dig_out_data,1), ...
                       'volt', 'yes', ...
                       'tempch', 'yes' };
else
    selected_output = cell(0);
    arg_named = @(x) ismember(varargin(1:2:end),x);
    if any(arg_named('ampch'))
        selected_output(1,end+(1:2)) = {'ampch', opts.do_AMP_channels};
    end
    if any(arg_named('auxch'))
        selected_output(1,end+(1:2)) = {'auxch', opts.do_AUX_channels};
    end
    if any(arg_named('ADCch'))
        selected_output(1,end+(1:2)) = {'ADCch', opts.do_ADC_channels};
    end
    if any(arg_named('DIGin'))
        selected_output(1,end+(1:2)) = {'DIGin', opts.do_DI_channels};
    end
    if any(arg_named('DIGout'))
        selected_output(1,end+(1:2)) = {'DIGout', opts.do_DO_channels};
    end
    if any(arg_named('volt'))
        selected_output(1,end+(1:2)) = {'volt', opts.do_VOLT_channels};
    end
    if any(arg_named('tempch'))
        selected_output(1,end+(1:2)) = {'tempch', opts.do_VOLT_channels};
    end
end

for arg=1:2:length(selected_output)
  switch (selected_output{arg})
    case 'ampch'
      if size(amplifier_data,1) > 0
        amplifier_channels = amplifier_channels(1,selected_output{arg+1});
        amplifier_data = amplifier_data(selected_output{arg+1},:);
        ts_amp = timeseries;
        ts_amp.name = 'ampch';
        ts_amp.time = t_amplifier;
        ts_amp.data = amplifier_data';
        if opts.time_roi ~=0
            ts_amp = select_time(ts_amp,opts.time_roi);
        end
        ts_amp.DataInfo.Units = 'uV';
        ts_amp.TimeInfo.Units = 'seconds';
        varargout{arg} = ts_amp;
        varargout{arg+1} = amplifier_channels;
      elseif opts.n_streams_req > 0
        disp('WARNING: No amplifier channels found in file!')
      end
    case 'auxch'
      if size(aux_input_data,1) > 0
        aux_input_channels = aux_input_channels(1,selected_output{arg+1});
        aux_input_data = aux_input_data(selected_output{arg+1},:);
        ts_aux = timeseries;
        ts_aux.name = 'auxch';
        ts_aux.time = t_aux_input;
        ts_aux.data = aux_input_data';
        if opts.time_roi~=0
            ts_aux = select_time(ts_aux,opts.time_roi);
        end
        ts_aux.DataInfo.Unit = 'uV';
        ts_aux.TimeInfo.Units = 'seconds';
        varargout{arg} = ts_aux;
        varargout{arg+1} = aux_input_channels;
      elseif opts.n_streams_req > 0
        disp('WARNING: No auxiliary inputs found in file!')
      end
    case 'ADCch'
      if size(board_adc_data,1)>0
        board_adc_channels = board_adc_channels(1,selected_output{arg+1});
        board_adc_data = board_adc_data(selected_output{arg+1},:);
        ts_board_adc = timeseries;
        ts_board_adc.name = 'ADCch';
        ts_board_adc.time = t_board_adc;
        ts_board_adc.data = board_adc_data';
        if opts.time_roi ~=0
            ts_board_adc = select_time(ts_board_adc,opts.time_roi);
        end
        ts_board_adc.DataInfo.Unit = 'uV';
        ts_board_adc.TimeInfo.Units = 'seconds';
        varargout{arg} = ts_board_adc;
        varargout{arg+1} = board_adc_channels;
      elseif opts.n_streams_req > 0
        disp('WARNING: No board ADC channels found in file!');
      end
    case 'DIGin'
      if size(board_dig_in_data,1)>0
        board_dig_in_channels = board_dig_in_channels(1,selected_output{arg+1});
        board_dig_in_data = board_dig_in_data(selected_output{arg+1},:);
        ts_dig_in = timeseries;
        ts_dig_in.name = 'DIGin';
        ts_dig_in.time = t_amplifier;
        ts_dig_in.data = board_dig_in_data';
        if opts.time_roi ~=0
            ts_dig_in = select_time(ts_dig_in,opts.time_roi);
        end
        ts_dig_in.TimeInfo.Units = 'seconds';
        varargout{arg} = ts_dig_in;
        varargout{arg+1} = board_dig_in_channels;
      elseif opts.n_streams_req > 0
        disp('WARNING: No board digital in channels found in file!');
      end
    case 'DIGout'
      if size(board_dig_out_data,1)>0
        board_dig_out_channels = board_dig_out_channels(1,selected_output{arg+1});
        board_dig_out_data = board_dig_out_data(selected_output{arg+1},:);
        ts_dig_out = timeseries;
        ts_dig_out.name = 'DIGout';
        ts_dig_out.time = t_amplifier;
        ts_dig_out.data = board_dig_out_data';
        if opts.time_roi ~=0
            ts_dig_out = select_time(ts_dig_out,opts.time_roi);
        end
        ts_dig_out.TimeInfo.Units = 'seconds';
        varargout{arg} = ts_dig_out;
        varargout{arg+1} = board_dig_out_channels;
      elseif opts.n_streams_req > 0
        disp('WARNING: No board digital out channels found in file!');
      end
    case 'volt'
      if size(supply_voltage_data,1)>0
        ts_volt = timeseries;
        ts_volt.name = 'supply_voltage';
        ts_volt.time = t_supply_voltage;
        ts_volt.data = supply_voltage_data';
        if opts.time_roi ~=0
            ts_volt = select_time(ts_volt,opts.time_roi);
        end
        ts_volt.TimeInfo.Units = 'seconds';
        varargout{arg} = ts_volt;
        varargout{arg+1} = supply_voltage_channels;
      elseif opts.n_streams_req > 0
        disp('WARNING: No supply voltage data found in file!');
      end
    case 'tempch'
      if size(temp_sensor_data,1)>0
        ts_temp = timeseries;
        ts_temp.name = 'temperature';
        ts_temp.time = t_temp_sensor;
        ts_temp.data = temp_sensor_data';
        if opts.time_roi ~= 0
            ts_temp = select_time(ts_temp,opts.time_roi);
        end
        ts_temp.TimeInfo.Units = 'seconds';
        varargout{arg} = ts_temp;
        varargout{arg+1} = '1';
      elseif opts.n_streams_req > 0
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

% Run filter`
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