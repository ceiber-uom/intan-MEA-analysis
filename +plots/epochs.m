

function epochs (data, varargin)
% function plots.epochs( data, ... )
%
% Options:
% -chan [c1 c2 c3 ...] % set channels 
% -map  [c1 c2; c3 c4; ...] % Set channels and arrangement of channels
% -sp   [rows cols]   % set subplot size (not needed if -map used)
%
% if any(named('-wave')) || ~isfield(...)
% end

named = @(s) strncmpi(s,varargin,numel(s));
get_ = @(v) varargin{find(named(v))+1};

if isfield(data,'config')
   this = @(d) plots.epochs(d, varargin{:});
   tools.forWaveTypes(data, this, varargin{:});
   return
end

channel_map = get_channelMap(get_, named, data);

if any(named('-roi')), t_roi = get('-roi'); 
    if numel(t_roi) == 2, error forach
    end
else t_roi = true(size(data.time)); 
end

for pp = 1:numel(channel_map)




end









return


function map = get_channelMap(get_, named, data)

nC = size(data.data,2);
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

if any(named('-sp')), spp = get_('-sp'); 
  if prod(spp) < numel(map)
    spp(2) = ceil(numel(map)/spp(1));
  end
elseif min(size(map)) == 1 && ~any(named('-map'))
  spp = ceil(sqrt(numel(map))); 
  spp(2) = ceil(numel(map) / spp);
else spp = []; 
end

if ~isempty(spp)
  map = reshape(map,spp(1),spp(2)); 
end

map = map'; % change to subplot index convention (rows)

map(map > nC) = 0; 
