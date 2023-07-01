
function epochs (data, varargin)
% function plots.epochs( data, ... )
%
% Basic plot of recorded signal waveform. If epoched data is supplied, one
% trace per epoch. 
% 
% Options:
% -chan [c1 c2 c3 ...] Set channels 
% -map  [c1 c2;  ... ] Set channels and arrangement of channels
% -rc   [rows cols]    Set subplot size (# rows / columns)
% -roi [t0 t1]         Set beginning and end of plotting window
% -dt  [x]             Set decimation factor for plotting. Otherwise, this
%                       is automatically determined based on the requested
%                       number of channels / samples. Integer only. 
% -dy [y]              Set vertical offset (as multiple of wave delta-z) 
% -ticks               Toggle default ticks behaviour (ticks enabled for
%                       6 or fewer channels, disabled otherwise)
% -single-figure       Plot all on one figure (old default)
% -pass [pass_ids]     Select passes to plot (for expoched data), 
% if any(named('-wave')) || ~isfield(...)
% end

named = @(s) strncmpi(s,varargin,numel(s));
get_ = @(v) varargin{find(named(v))+1};

if isfield(data,'config')
   this = @(d) plots.epochs(d, varargin{:});
   tools.forWaveType(data, this, varargin{:});
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
opts.one_figure = any(named('-single-f'));

%%
if numel(channel_map) > 1, clf, end

time = data.time(t_roi); 
dsf = ceil(size(data.wave,2)^0.5 * length(time)^0.1);
if any(named('-dt')), dsf = get('-dt'); 
  if dsf < 1, dsf = ceil(1/dsf); 
  else dsf = round(dsf); 
  end
end
time = time(1:dsf:end);

dy = median(max(max(data.wave))); 
if any(named('-dy')), dy = dy * get_('-dy'); 
% else dy = dy;
end

if any(named('-pass')), passes = get_('-pass');
else passes = true(size(data.wave,3),1);
end

C = lines(7); 

for pp = 1:numel(channel_map)

    cc = channel_map(pp);
    if cc == 0, continue, end
    if opts.one_figure, figure('Name',sprintf('channel %d', cc)), cla
    elseif numel(channel_map) > 1, subplot(subplot_nxy{:}, pp), 
    else cla reset
    end

    y = permute(data.wave(t_roi,cc,passes), [1 3 2]);
    y = y(1:dsf:end,:) + dy*(0:size(y,2)-1);

    plot(time, y,'Color',C(1,:)); % ,'Clipping','off')

    xlim(time([1 end])),
    plots.tidy
    p = get(gca,'Position'); 
    set(gca,'Position',p + [-1 -1 2 2].*p([3 4 3 4])/10)
    set(gca,'UserData',cc)
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

