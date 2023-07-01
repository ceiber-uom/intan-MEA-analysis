
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
% -per-unit            One panel per unit (default: one per channel)
% -isi                 View raster as time versus ISI (use plots.ISI)
% -shape [n]           Colour spikes by the Nth column of 'data.shape'
%                      (see tools.readPlexon and the source data file)
% 
% v0.1 - 4 June 2023 - Calvin Eiber <c.eiber@ieee.org>


named = @(s) strncmpi(s,varargin,numel(s));
get_ = @(v) varargin{find(named(v))+1};

if isfield(data,'config')
   this = @(d) plots.raster(d, varargin{:});
   tools.forWaveType(data, this, varargin{:});
   return
elseif numel(data) > 1 || any(named('--undo')) % e.g. the output from tools.simplify()
  % The --undo flag is added to force this behaviour in the case of plotting just one cell  
    data = tools.simplify(data,'-undo');
end

if any(named('-isi'))
    plots.ISI(data, varargin{:}, '-raster');
    return
end

if any(named('-per-u'))
    [~,~,data.channel] = unique([data.channel data.unit],'rows');
end

channel_map = plots.layout(data, varargin{:});

if any(named('-roi')), time_window = get('-roi'); 
else time_window = []; 
end

include = true(size(data.time)); 

if any(named('-unit')), u_roi = get('-unit'); 
  if isnumeric(u_roi), u_roi = ismember(data.unit, u_roi); end
  if any(u_roi), include = include & u_roi; end
end

u_max = double(max(data.unit(:)));

do_shape_var = [];
if any(named('-shape')), do_shape_var = get_('-shape'); end

%%
if ~isfield(data,'pass') || any(named('-no-e'))
    %% Basic raster plot without epoching


    if numel(time_window) >= 2
        include = include & (data.time >= min(include) & ...
                             data.time <= max(include));
    end

    
    include = include & ismember(data.channel,channel_map(:));

    y = double(data.channel(include)); 
    if ~any(named('-ex')), 
        unit_offset = double(data.unit(include))./(u_max+1)-0.5;
        y = y + unit_offset; 
    end

    u = data.unit(include);
    if ~isempty(do_shape_var), u = data.shape(include,do_shape_var); end

    cla reset
    scatter(data.time(include), y, [], u, '.')
    if isempty(do_shape_var), caxis([0 u_max]),
      if u_max > 7
           colormap([.5 .5 .5; turbo(u_max)]); 
      else colormap([.5 .5 .5; lines(u_max)]); 
      end
    end

    v_ = @(x) reshape(x,[],1);
    yl = double([min(v_(channel_map(channel_map>0))) max(channel_map(:))]);

    plots.tidy, ylim(yl + [-1.5 0.5])

    xlabel('time, s'), ylabel('channel ID')
    if isempty(do_shape_var), ylabel(colorbar,'unit ID'); end
    return
end

%% One-per-channel raster plots

if any(named('-pass')), pass_ok = get_('-pass');
elseif ~isfield(data,'pass')
     data.pass = ones(size(data.time));
     data.pass_begin = 0; 
     pass_ok = true;
else pass_ok = true(max(data.pass),1);
end
if ~islogical(pass_ok)
    pass_ok = ismember(1:max(data.pass), pass_ok); 
end

pass_ok = [false; pass_ok];
pass_t0 = [0; data.pass_begin];

subplot_nxy = num2cell(size(channel_map));
opts.ticks = (numel(channel_map) > 6) == any(named('-ti'));

if numel(channel_map) > 1, clf, end
for pp = 1:numel(channel_map)

    if channel_map(pp) == 0, continue, end

    cc = channel_map(pp); 
    ok = (include & data.channel == cc & pass_ok(data.pass+1)); 
    
    if numel(time_window) >= 2
        if ~any(ok), continue, end
        t = data.time(ok) - pass_t0(data.pass(ok)+1);
        ok(ok) = (t >= min(time_window) & t <= max(time_window));
    end

    if ~any(ok), continue, end
    if numel(channel_map) > 1, subplot(subplot_nxy{:}, pp), 
    else cla reset
    end

    y = double(data.pass(ok)); 
    if any(named('-ex')), y = y+double(data.unit(ok))/u_max-0.5; end 

    u = double(data.unit(ok));
    if ~isempty(do_shape_var), u = data.shape(ok,do_shape_var); end

    t = data.time(ok) - pass_t0(data.pass(ok)+1);

    

    scatter(t, y, [], u, '.')
    if isempty(do_shape_var), 
      caxis([0 u_max]), colormap([.5 .5 .5; lines(u_max)]); 
    end

    plots.tidy
    
    p = get(gca,'Position'); 
    set(gca,'Position',p + [-1 -1 2 2].*p([3 4 3 4])/10)
    set(gca,'UserData',cc)
    
    if ~opts.ticks, set(gca,'XTick',[],'YTick',[]); end
end

h = get(gcf,'Children');

set(h,'YLim',[0.5 max(abs([h.YLim]))])
linkaxes(h,'xy')
pos = cat(1,h.Position);

[~,blc] = min(pos*[1;1;0;0]);
set(h(blc),'XTickMode','auto');
set(h,'LineWidth',0.8)

%%
return

