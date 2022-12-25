

function map = layout(data, varargin)
% channel_map = plots.layout(data, ... )
% implement standard channel mapping syntax
% 
% typical usage: 
% channel_map = plots.layout(data, varargin{:});
% 
% Options:
% -chan [c1 c2 c3 ...] % set channels 
% -map  [c1 c2; c3 c4; ...] % Set channels and arrangement of channels
% -rc   [rows cols]   % set subplot size (not needed if -map used)

named = @(s) strncmpi(s,varargin,numel(s));
get_ = @(v) varargin{find(named(v))+1};

if isfield(data,'data'),        nC = size(data.data,2);
elseif isfield(data,'channel'), nC = max(data.channel(:));
else error('unable to determine channels')
end

if nC == 64
     map = [ 1  2  3  4  5  6  7  8;
             9 10 11 12 13 14 15 16;
            17 18 19 20 21 22 23 24;
            25 26 27 28 29 30 31 32; 
            33 34 35 36 37 38 39 40;
            41 42 43 44 45 46 47 48;
            49 50 51 52 53 54 55 56;
            57 58 59 60 61 62 63 64]; 
else map = 1:nC;
end

if any(named('-map')), map = get_('-map');
elseif any(named('-ch')), map = get_('-ch');
end

if any(named('-rc')), sprc = get_('-rc'); 
  if prod(sprc) < numel(map)
    sprc(2) = ceil(numel(map)/sprc(1));
  end
elseif min(size(map)) == 1 && ~any(named('-map'))
  sprc = ceil(sqrt(numel(map))); 
  sprc(2) = ceil(numel(map) / sprc);
  map(end+1 : prod(sprc)) = 0; 
else sprc = []; 
end

if ~isempty(sprc)
  map = reshape(map,sprc(1),sprc(2)); 
end

map = map'; % change to subplot index convention (rows)
map(map > nC) = 0; % mask invalid values with 0
map(map < 0) = 0;
