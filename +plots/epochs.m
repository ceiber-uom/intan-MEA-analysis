

function epochs (data, varargin)
% function plots.epochs( data, ... )

% Options:
% -chan [c1 c2 c3 ...] % set channels 
% -map  [c1 c2; c3 c4; ...] % Set channels and arrangement of channels
% -sp   [rows cols]   % set subplot size (not needed if -map used)

% if any(named('-wave')) || ~isfield(...)
% end

named = @(s) strncmpi(s,varargin,numel(s));
get_ = @(v) varargin{find(named(v))+1};

channel_map = [ 1  2  3  4  5  6  7  8;
                9 10 11 12 13 14 15 16;
               17 18 19 20 21 22 23 24;
               25 26 27 28 29 30 31 32; 
               33 34 35 36 37 38 39 40;
               41 42 43 44 45 46 47 48;
               49 50 51 52 53 54 55 56;
               57 58 59 60 61 62 63 64]; 

if any(named('-map')), channel_map = get_('-map');
elseif any(named('-ch')), channel_map = get_('-ch');
  if any(named('-sp')), spp = get_('-sp'); 
    if prod(spp) < numel(channel_map)
      spp(2) = ceil(numel(channel_map)/spp(1));
    end
  elseif min(size(ch)) == 1
    spp = ceil(sqrt(numel(channel_map))); 
    spp(2) = ceil(numel(channel_map) / spp);
  else spp = []; 
  end
  if ~isempty(spp)
    channel_map = reshape(channel_map,spp(1),spp(2)); 
  end
end

channel_map = channel_map'; % change to subplot index convention (rows)

if any(named('-roi')), t_roi = get('-roi'); 
    if numel(t_roi) == 2, t_roi = []; 
    end
else t_roi = true(size(data.time)); 
end

for pp = 1:numel(channel_map)




end


