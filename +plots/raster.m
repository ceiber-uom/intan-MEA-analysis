
function raster (data, varargin)
% function plots.raster( data, ... )
%
% Basic raster plot of recorded spike times. 
% If epoched spiketime data is supplied, one row per epoch and one axis
% per channel, otherwise plot all on a single axis.
% 
% Options:
% -chan [c1 c2 c3 ...] Set channels 
% -map  [c1 c2;  ... ] Set channels and arrangement of channels
% -unit [u1 u2 ... ]   Filter the displayed units 
% -rc   [rows cols]    Set subplot size (# rows / columns)
% -roi [t0 t1]         Set beginning and end of plotting window
% -ticks               Toggle default ticks behaviour (ticks enabled for
%                       6 or fewer channels, disabled otherwise)
% -pass [pass_ids]     Select passes to plot (for expoched data), 
% -no-epochs           Plot all spikes on a single axes, time vs. channel
% -expand-units        Expand out the different units (for visibility)
%                      enabled by default if epoch mode disabled. 
% 
% v0.1 - 23 December 2022 - Calvin Eiber <c.eiber@ieee.org>

named = @(s) strncmpi(s,varargin,numel(s));
get_ = @(v) varargin{find(named(v))+1};

if isfield(data,'config')
   this = @(d) plots.raster(d, varargin{:});
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

if any(named('-unit')), u_roi = get('-unit'); 
  if isnumeric(u_roi), u_roi = ismember(data.unit, u_roi); end
  if any(u_roi), t_roi = t_roi & u_roi; end
end

u_max = max(data.unit(:));

%%
if ~isfield(data,'epoch') || ~any(named('-no-e'))
    %% Basic raster plot
    
    y = data.channel(t_roi); 
    if ~any(named('-ex')), 
        y = y + data.unit(t_roi)/(u_max+1)-0.5; 
    end

    cla reset
    scatter(data.time(t_roi), y, [], data.unit(t_roi), '.')
    colormap([.5 .5 .5; lines(u_max)])
    caxis([0 u_max])
    plots.tidy

    xlabel('time, s'), ylabel('channel ID')
    ylabel(colorbar,'unit ID')
    return
end

%% One-per-channel raster plots

if any(named('-pass')), pass_ok = get_('-pass');
elseif ~isfield(data,'pass')
     data.pass = ones(size(data.time));
     pass_ok = true;
else pass_ok = true(size(data.data,3),1);
end
if ~islogical(pass_ok)
    pass_ok = ismember(1:max(data.pass), pass_ok); 
end

subplot_nxy = num2cell(size(channel_map));
opts.ticks = (numel(channel_map) > 6) == any(named('-ti'));

if numel(channel_map) > 1, clf, end
for pp = 1:numel(channel_map)

    if channel_map(pp) == 0, continue, end

    cc = channel_map(pp); 
    ok = (t_roi & data.channel == cc & pass_ok(data.pass)); 
    
    if ~any(ok), continue, end
    if numel(channel_map) > 1, subplot(subplot_nxy{:}, pp), 
    else cla reset
    end

    y = data.pass(ok); 
    if any(named('-ex')), y = y+data.unit(ok)/u_max-0.5; end 

    scatter(data.time(ok), y, [], data.unit(ok), '.')
    colormap([.5 .5 .5; lines(u_max)])
    caxis([0 u_max])
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

