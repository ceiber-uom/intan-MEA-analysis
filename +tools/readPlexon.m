function data = readPlexon(filename, varargin)
% data = tools.readPlexonCSV([data], filename, varargin)
% function header go here
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
    [n,p] = uigetfile('*.csv',[],filename);
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

AMP = 'AMP';
if any(named('-field')), AMP = get_('-field'); end
if ~isfield(data,AMP), data.(AMP) = struct; 
elseif ~isstruct(data.(AMP)), AMP = [AMP '_spikes']; 
    if ~isfield(data,AMP), data.(AMP) = struct; end
end

data.(AMP).time    = spikes.Timestamp;
data.(AMP).channel = spikes.Channel;
data.(AMP).unit    = spikes.Unit;
data.(AMP).shape   = spikes{:,4:end};

return