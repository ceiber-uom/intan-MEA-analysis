function data = readKS(filename, varargin)
% data = tools.readKS([data], filename, varargin)
% 
% Read spike data as exported by KiloSort and append it to the input
% data structure (or generate a new data structure, compatible with tools
% that work with the output of readIntan). Adapted from
% selectSortedClusters code 
% 
% Options: 
% -quality [4]    - select cells with quality value <= Q (1 = best)
% -min-count [10] - select only cells with at least N spikes
% -field [SPIKE]  - append spike data to what INTAN wave type?
%                   .SPIKE is the default for this data. 
% 
% v0.1 - 4 June 2023 - Calvin Eiber <c.eiber@ieee.org>

named = @(n) strncmpi(varargin,n,length(n));
get_ = @(v) varargin{find(named(v))+1};

%% Get filename
if nargin == 0, filename = '?'; end
if isstruct(filename), data = filename; 
  if nargin > 2 && exist(varargin{1},'file'), filename = varargin{2};
  else ks_folder = fileparts(data(1).filename);
       filename = [ks_folder '/ks_sorted/ksrasters/ksrasters.mat'];
  end
end

% IF the filename isn't a good file, use a UI to pick a file. 
if ~exist(filename,'file')
  filename = strrep(filename,'.mat','*.mat');
  if isempty(dir(filename))
    [n,p] = uigetfile({filename;'ksrasters.mat';'*.mat'},[],filename);
    filename = [p n];
  else
    list = dir(filename);
    p_ = @(x) [x(1).folder filesep x(1).name]; % path expander
    [~,sel] = max([list.datenum]);
    filename = p_(list(sel)); 
  end
end
clear list sel n p p_     

% load the KS rasters
disp(['Reading ' filename])
spike = load(filename);

quality_threshold = 4;
if any(named('-q')), quality_threshold = get_('-q'); end

minimum_count = 10; 
if any(named('-m')), minimum_count = get_('-m'); end

%%

if ~exist('data','var')
    data = struct;
    data.filename = filename;
    data.notes  = struct;
    data.config = struct; 
end

data.notes.spikes = sprintf('Loaded from %s on %s with %s', ...
                             filename, datestr(now),'tools.readKS');

field = 'SPIKE';
if any(named('-field')), field = get_('-field'); end
if ~isfield(data,field), data.(field) = struct; 
elseif ~isstruct(data.(field)), field = [field '_spikes']; 
    if ~isfield(data,field), data.(field) = struct; end
end

% 3rd column is the rating of clusters
ok = spike.clusters(:,3) <= quality_threshold;

if any(~ok), fprintf('%d spikes removed %s %d (%d/%d remaining)\n', ...
                     sum(~ok), 'using -quality', quality_threshold, ...
                     sum(ok), numel(ok))
end

nS = sum(cellfun(@numel,spike.spike_times),2);

lc = (nS < minimum_count);
if any(lc(ok)), fprintf('%d spikes removed %s %d (%d/%d remaining)\n', ...
                     sum(lc(ok)), 'using -min-count', minimum_count, ...
                     sum(ok & ~lc), sum(ok))
end
ok(lc) = false; 


% The sorted spikes are chunked minute by minute, so unchunk. 

timing_info = spike.sort_params.stim_start_end / spike.fs;
timing_info = timing_info - timing_info(1);

spike.unit_id = cell(size(spike.amplitudes));

for ii = 1:numel(spike.unit_id)
    [kk,tt] = ind2sub(size(spike.unit_id), ii);
    spike.unit_id{ii} = uint16(0*spike.amplitudes{ii} + kk);
    spike.spike_times{ii} = spike.spike_times{ii} + timing_info(tt,1);
end

unit_channel_id = uint16([spike.sort_info.ch]); 

% Build output structure

data.(field).time    = cat(1,spike.spike_times{ok,:});
data.(field).channel = []; 
data.(field).unit    = cat(1,spike.unit_id{ok,:});
data.(field).channel = unit_channel_id(data.(field).unit)';
data.(field).shape   = cat(1,spike.amplitudes{ok,:});

data.(field).clusters = spike.clusters(ok);
data.(field).sort_info   = spike.sort_info(ok);
data.(field).template_info = spike.template_info(ok);

data.config.(['KS_' field]) = spike.sort_params;

data.config.shape_variables = {'amplitude'};

return