function data = readPlexon(filename, varargin)
% data = tools.readPlexon([data], filename, varargin)
% 
% Read tabular spike data as exported by Plexon and append it to the input
% data structure (or generate a new data structure, compatible with tools
% that work with the output of readIntan). 
% 
% The input table is required to have the "Timestamp", "Channel", and
% "Unit" fields. 
% 
% Options: 
% -field [SPIKE] - append spike data to what INTAN wave type?
%                  .SPIKE is the default for this data. 
% 
% -field [AMP]

named = @(n) strncmpi(varargin,n,length(n));
get_ = @(v) varargin{find(named(v))+1};

%% Get filename
if nargin == 0, filename = '?'; end
if isstruct(filename), data = filename; 
  if nargin > 2 && exist(varargin{1},'file')
       filename = varargin{2};
  else filename = strrep(data.filename,'.rhd','.csv');
  end
end

if ~exist(filename,'file')
  filename = strrep(filename,'.csv','*.csv');
  if isempty(dir(filename))
    [n,p] = uigetfile({'*.csv';'*.xls;*.xlsx'},[],filename);
    filename = [p n];
  else
    list = dir(filename);
    p_ = @(x) [x(1).folder filesep x(1).name]; % path expander
    [~,sel] = max([list.datenum]);
    filename = p_(list(sel)); 
  end
end
clear list sel n p p_     

disp(['Reading ' filename])
spikes = readtable(filename);

%%

if ~exist('data','var')
    data = struct;
    data.filename = filename;
    data.notes.spikes = sprintf('Loaded from %s on %s', filename, datestr(now));
    data.config = struct; 
end

field = 'SPIKE';
if any(named('-field')), field = get_('-field'); end
if ~isfield(data,field), data.(field) = struct; 
elseif ~isstruct(data.(field)), field = [field '_spikes']; 
    if ~isfield(data,field), data.(field) = struct; end
end

[shape_vars,svi] = setdiff(spikes.Properties.VariableNames, ...
                     {'Timestamp','Channel','Unit'});

data.(field).time    = spikes.Timestamp;
data.(field).channel = spikes.Channel;
data.(field).unit    = spikes.Unit;
data.(field).shape   = spikes{:,svi};

data.config.shape_variables = shape_vars;

return