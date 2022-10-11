

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
    if numel(t_roi) == 2
        t_roi = (data.time >= min(t_roi) & ...
                 data.time <= max(t_roi));
    end
else t_roi = true(size(data.time)); 
end

subplot_nxy = num2cell(size(channel_map));

opts.ticks = (numel(channel_map) > 6) == any(named('-ti'));
%%
clf

time = data.time(t_roi); 
dsf = ceil(size(data.data,2)^0.5 * length(time)^0.1);
if any(named('-dt')), dsf = get('-dt'); end
time = time(1:dsf:end);

dy = median(max(max(data.data))); 
if any(named('-dy')), dy = dy * get_('-dy'); 
% else dy = dy;
end

C = lines(7); 

for pp = 1:numel(channel_map)

    if channel_map(pp) == 0, continue, end
    subplot(subplot_nxy{:}, pp), cc = channel_map(pp); 

    y = permute(data.data(t_roi,cc,:), [1 3 2]);
    y = y(1:dsf:end,:) + dy*(0:size(y,2)-1);

    plot(time, y,'Color',C(1,:)); % ,'Clipping','off')

    xlim(time([1 end])),
    tidyPlotForIllustrator
    p = get(gca,'Position'); 
    set(gca,'Position',p + [-1 -1 2 2].*p([3 4 3 4])/10)
    if ~opts.ticks, set(gca,'XTick',[],'YTick',[]); end
end

h = get(gcf,'Children');

if size(y,2) > 1 && dy > 0
     set(h,'YLim',[-0.75 size(y,2)-0.5]*dy)
else set(h,'YLim',[-1.01 1]*max(abs([h.YLim])));
end

linkaxes(h,'xy')
pos = cat(1,h.Position);

[~,blc] = min(pos*[1;1;0;0]);
set(h(blc),'XTickMode','auto');
set(h,'LineWidth',0.8)

%%
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
