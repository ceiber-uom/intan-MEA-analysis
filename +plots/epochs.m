

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

channel_map = plots.layout(data, varargin{:});

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

