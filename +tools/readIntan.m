function output = readIntan(filename, varargin)
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
% 'roi' [t_start t_end] : Return a specified range of data (in sec)
% 'amp'  [list] : Select amplifier channnels to load
% 'aux'  [list] : Select auxilary channels to load
% 'volt' [true] : Load the supply voltage channel?
% 'ADC'  [list] : Select on-board ADC channels to load
% 'DI'   [list] : Select onboard digital input channels to load
% 'DO'   [list] : Select onboard digital output channels to load
% 'temp' [true] : Should temperature data should be loaded? (logical)
% 'notch' [50]: Assign '50Hz' or '60Hz' to enable notch filtering. If left
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
opts.time_roi = [];
opts.export_ts = any(named('-timeseries')) || any(named('-ts'));

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

do_ch_ = @(t) ['do_' upper(t{1}) '_channels'];

for chan_type = channel_types
  opts.(do_ch_(chan_type)) = false; 
  if any(named(chan_type{1}))
    opts.(do_ch_(chan_type)) = get_(chan_type{1});  
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
time_ = @(ty) [upper(ty{1}) '_time'];

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
assert(data_present,'this file does not contain data to load!')

data = struct; 
for ty = channel_types
    data.(chan_(ty)) = zeros(num.(chan_(ty)), num.(sams_(ty)));
    data.(sams_(ty)) = 1; 
end

data.AMP_time = zeros(1, num.AMP_samples);

data.DI_raw = zeros(1,num.DI_channels);
data.DO_raw = zeros(1,num.DO_channels);

% Read sampled data from file.

printInfo(); 

for ii = 1:num_data_blocks

    % Make console progress bar
    console_progress_bar(ii / num_data_blocks);
    
    % In version 1.2, we moved from saving timestamps as unsigned
    % integeters to signed integers to accomidate negative (adjusted)
    % timestamps for pretrigger data.
    if (data_file_main_version_number > 1) || ...
       (data_file_main_version_number == 1 && ...
        data_file_secondary_version_number >= 2)
         data.AMP_time( data.AMP_samples + (0 : nspd-1) ) = fread(fid, nspd, 'int32');
    else data.AMP_time( data.AMP_samples + (0 : nspd-1) ) = fread(fid, nspd, 'uint32');
    end

    if (num.AMP_channels > 0)
      data.AMP_channels(:, data.AMP_samples + (0 : nspd-1)) = ...
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
else disp('Done.')
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
num_gaps = sum(diff(data.AMP_time) ~= 1);
if (num_gaps == 0)
else fprintf('Warning: %d gaps in timestamp data found. %s\n', ...
                       num_gaps, 'Time scale will not be uniform!');
end

% Scale time steps (units = seconds).
data.AMP_time  = data.AMP_time / sample_rate;
data.AUX_time  = data.AMP_time(1:4:end);
data.VOLT_time = data.AMP_time(1:60:end);
data.ADC_time  = data.AMP_time;
data.DI_time   = data.AMP_time;
data.DO_time   = data.AMP_time;
data.TEMP_time = data.VOLT_time;

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

% as defined above:
% channel_types = {'amp','aux','volt','adc','DI','DO','temp'};
% opt_id = ['do_' upper(chan_type{1}) '_channels'];

for ty = channel_types

    sel = opts.(do_ch_(ty)); % selection
    n_c = numel(channels.(chan_(ty))); 

    % First, check to see if we skip this output type 
    if opts.n_streams_req > 0 && ~any(opts.(do_ch_(ty))), continue, end
    if opts.n_streams_req == 0, sel = 1:n_c;
      if ~any(sel), continue, end
    elseif n_c == 0
      sprintf('WARNING: No %s channels found in file!\n', ty{1})
      continue
    end
 
    % parse selection applied as necessary
    if numel(sel) == 1 && n_c > 1 && islogical(sel) && sel
        sel = 1:n_c;
    elseif isnumeric(sel), sel(sel>n_c | sel<1) = []; 
    elseif ischar(sel) 
        if strcmpi(sel,'off'), continue, end
        sel = 1:n_c;
    end

    output.config.(chan_(ty)) = channels.(chan_(ty))(sel);

    if opts.export_ts, this = timeseries;
    else this = struct; 
    end

    this.name = upper(ty{1});
    this.time = data.(time_(ty));
    this.data = data.(chan_(ty))(sel,:)';

    if ~isempty(opts.time_roi)
        this = select_time(this, opts.time_roi);
    end

    this.TimeInfo.Units = 'seconds';
    switch ty{1}
        case 'amp',  this.DataInfo.Units = 'uV';
        case 'temp', this.DataInfo.Units = 'deg C';
        otherwise,   this.DataInfo.Units = 'V';
    end

    this.DataInfo.Channels = channels.(chan_(ty)); 
    output.(upper(ty{1})) = this;
end

if nargout == 0
    assignin('caller','intanData',output)
    clear
end

return

function new = select_time(old,time_range)
% Enter a timeseries and a range.
% Returns a new timeseries according to range.
% Range should be a 1x2 matrix with [start end]
if isstruct(old), new = old([]);
else new = timeseries;
end

sel = old.time >= min(time_range) & old.time <= max(time_range); 

new.name = old.name;
new.time = old.time(sel);
new.data = old.data(sel,:);

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

function console_progress_bar( frac )

w = 35;
lhs = repmat('#',1,floor(frac * w));

if frac == 1, printInfo('[%s] ',lhs), return, end

rhs = repmat(' ',1,w-numel(lhs)-1);
ff = floor(10*mod(frac*w,1));

printInfo('[%s%d%s] ',lhs,ff,rhs)