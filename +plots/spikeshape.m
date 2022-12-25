
function spikeshape (data, varargin)
% function plots.spikeshape( data, ... )
%
% Basic plots of recorded spike shapes. 
% If epoched spiketime data is supplied, one row per epoch and one axis
% per channel, otherwise plot all on a single axis.
% 
% Options:
% -chan [c1 c2 c3 ...] Set channels 
% -map  [c1 c2;  ... ] Set channels and arrangement of channels
% -rc   [rows cols]    Set subplot size (# rows / columns)
% -roi [t0 t1]         Set beginning and end of plotting window
% -ticks               Toggle default ticks behaviour (ticks enabled for
%                       6 or fewer channels, disabled otherwise)
% -pass [pass_ids]     Select passes to plot (for expoched data), 
% -n [24]              Set number of example waves. 
%                      If 0, show average only +- SD
% -no-hash             Only show spikes with unit code > 0. 
% -mean                Show mean +- standard deviation 

named = @(s) strncmpi(s,varargin,numel(s));
get_ = @(v) varargin{find(named(v))+1};

if isfield(data,'config')
   this = @(d) plots.spikeshape(d, varargin{:});
   tools.forWaveType(data, this, varargin{:});
   return
end

if ~isfield(data,'shape')
    error('spike shape data not found'), 
end

channel_map = plots.layout(data, varargin{:});

if any(named('-time')), time = get('-time'); 
elseif size(data.shape,2) == 32
     time = ((1:32)-8) / 20000 ; % default time
else time = (1:size(data.shape,2)) / 20000 ; % default time
end
if numel(time) ~= size(data.shape,2)
    error('number of TIME samples (%d) must match %s (%d)\n%s', ...
           numel(time), 'number of samples for spike.shape.', ...
           size(data.shape,2), 'Use -time [t] to set.')
else time = reshape(time,1,[]);
end

if any(named('-roi')), include = get('-roi'); 
    if numel(include) == 2
        include = (data.time >= min(include) & ...
                   data.time <= max(include));
    end
else include = true(size(data.time)); 
end

u_max = max(data.unit(:));
C = lines(u_max); 

n_examples = 24;
if any(named('-n')), n_examples = get_('-n'); end
if any(named('-mean')), n_examples = 0; end

do_hash = ~any(named('-no-hash'));

%% One-per-channel spike shape plots

if any(named('-pass')), pass_ok = get_('-pass');
elseif ~isfield(data,'pass')
     data.pass = ones(size(data.time));
     pass_ok = true;
else pass_ok = true(max(data.pass),1);
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
    this_channel = (include & data.channel == cc & pass_ok(data.pass)); 
    if ~do_hash, this_channel = this_channel & data.unit > 0; end

    if ~any(this_channel), continue, end
    if numel(channel_map) > 1, subplot(subplot_nxy{:}, pp), 
    else cla reset
    end, hold on

    n_units = max(data.unit(this_channel));

    for u_id = 0:n_units

        ok = this_channel & data.unit == u_id;
        if u_id == 0, color = [.5 .5 .5];
        else color = C(u_id,:); 
        end

        if ~any(ok), continue, end

        avg_wave = mean(data.shape(ok,:),1);        
    
        if n_examples == 0
    
            sd = std(data.shape(ok,:),[],1);
            fill([time fliplr(time)], [avg_wave+sd fliplr(avg_wave-sd)], ...
                  color,'EdgeColor','none','FaceAlpha',0.3);
            plot(time, avg_wave, '-', 'Color', color,'LineWidth',1.2)
            continue
    
        elseif any(named('-rand'))
             [~,idx] = sort(rand(sum(ok),1));
             idx = idx(1:max(end,n_examples));
        else idx = unique(round(linspace(1,sum(ok),n_examples)));
        end
    
        ok_ids = find(ok); idx = ok_ids(idx); 

        plot(time, data.shape(idx,:),'Color', [color 0.5])
        plot(time, avg_wave, '-', 'Color', color,'LineWidth',1.2)    
    end

    plots.tidy    
    p = get(gca,'Position'); 
    set(gca,'Position',p + [-1 -1 2 2].*p([3 4 3 4])/10)
    set(gca,'UserData',cc)
    if ~opts.ticks, set(gca,'XTick',[],'YTick',[]); end

    


end

h = get(gcf,'Children');
% linkaxes(h,'xy')
pos = cat(1,h.Position);

[~,blc] = min(pos*[1;1;0;0]);
set(h(blc),'XTickMode','auto');
set(h,'LineWidth',0.8)

%%
return

