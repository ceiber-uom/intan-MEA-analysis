
function ISI (data, varargin)
% function plots.ISI( data, ... )
%
% Basic inter-spike-inteval histogram plot of recorded spike times. 
% 
% Options:
% -chan [c1 c2 c3 ...] Set channels 
% -map  [c1 c2;  ... ] Set channels and arrangement of channels
% -rc   [rows cols]    Set subplot size (# rows / columns)
% -roi  [t0 t1]        Set beginning and end of plotting window
% -unit [u1 u2 ... ]   Filter for the specified units 
% -pass [pass_ids]     Select passes to plot (for expoched data)
% -ignore-unit         Compute ISI for all spikes on the channel
% -no-hash             Only show spikes with unit code > 0. 
% -t [bin-times]       Bin times for ISI histogram (default 0:0.5 s, n=64)
% -labels              Add channel labels (default: corners only)
% -no-labels           Suppress channel labels for corners
% 
% v0.1 - 23 December 2022 - Calvin Eiber <c.eiber@ieee.org>

named = @(s) strncmpi(s,varargin,numel(s));
get_ = @(v) varargin{find(named(v))+1};

if isfield(data,'config')
   this = @(d) plots.ISI(d, varargin{:});
   tools.forWaveType(data, this, varargin{:});
   return
end

channel_map = plots.layout(data, varargin{:});

if any(named('-roi')), include = get('-roi'); 
    if numel(include) == 2
        include = (data.time >= min(include) & ...
                 data.time <= max(include));
    end
else include = true(size(data.time)); 
end

subplot_nxy = num2cell(size(channel_map));
opts.ticks = (numel(channel_map) > 6) == any(named('-ti'));
opts.hash  = ~any(named('-no-hash'));
opts.merge = any(named('-ignore-u'));
opts.log   = any(named('-log')); 

% Add channel labels
opts.label_all = any(named('-label'));
opts.label_none = any(named('-no-l'));

if any(named('-unit')), u_roi = get('-unit'); 
  if isnumeric(u_roi), u_roi = ismember(data.unit, u_roi); end
  if any(u_roi), include = include & u_roi; end
elseif opts.merge
  if opts.hash, data.unit(:) = 1; 
  else data.unit = sign(data.unit);
  end
end

if any(named('-pass')), pass_ok = get_('-pass');
elseif ~isfield(data,'pass')
     data.pass = ones(size(data.time));
     pass_ok = true;
else pass_ok = true(size(data.data,3),1);
end
if ~islogical(pass_ok), pass_ok = ismember(1:max(data.pass), pass_ok); end
if any(pass_ok), include = include & pass_ok; end

u_max = max(data.unit(:));
C = lines(u_max); 

bins = linspace(0, 0.5, 64); 
if any(named('-log')), bins = linspace(-3, 1, 64); end
if any(named('-t')), bins = get_('-t');
  if numel(bins) == 1, bins = linspace(0, bins, 64);
  elseif numel(bins) == 2, bins = linspace(0, bins(1), bins(2));
  end
end

if numel(channel_map) > 1, clf, end

for pp = 1:numel(channel_map)

    if channel_map(pp) == 0, continue, end
    
    cc = channel_map(pp);
    this_channel = (include & data.channel == cc & pass_ok(data.pass)); 
    if ~opts.hash, this_channel = this_channel & data.unit > 0; end

    if ~any(this_channel), continue, end
    if numel(channel_map) > 1, subplot(subplot_nxy{:}, pp), 
    else cla reset
    end, hold on

    n_units = max(data.unit(this_channel));

    yl = 1; 

    for u_id = 0:n_units

        ok = this_channel & data.unit == u_id;
        if u_id == 0, color = [.5 .5 .5];
        else color = C(u_id,:); 
        end

        if ~any(ok), continue, end

        isi = diff(data.time(ok)); 

        if opts.log, [y,x] = hist(log10(isi), bins); %#ok<HIST> 
        else         [y,x] = hist(isi, bins);        %#ok<HIST> 
        end

        style = {'FaceColor',color,'FaceAlpha',0.5,'EdgeColor','none'};
        if opts.merge, style{4} = 1; end

        bar(x,y,1,style{:})
        yl = max([yl y(2:end-1)]);
    end

    ylim([0 yl]), plots.tidy    
    p = get(gca,'Position'); 
    set(gca,'Position',p + [-1 -1 2 2].*p([3 4 3 4])/10)
    set(gca,'UserData',cc)
    if ~opts.ticks, set(gca,'XTick',[],'YTick',[]); end

    % Set labels
    if opts.label_none, continue, end
    [a,b] = find(channel_map == cc);
    do_label = opts.label_all || ...
         (a == 1 || a == size(channel_map,1)) && ...
         (b == 1 || b == size(channel_map,2));

    if do_label
        text(mean(xlim),mean(ylim),num2str(cc),'Color',[.4 .4 .4], ...
                                 'FontSize',14,'horiz','center')
    end
end

h = get(gcf,'Children');
% linkaxes(h,'xy')
pos = cat(1,h.Position);

[~,blc] = min(pos*[1;1;0;0]);
set(h(blc),'XTickMode','auto');
set(h,'LineWidth',0.8)

return

