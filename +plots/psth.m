
function psth (data, varargin)
% function plots.psth( data, ... )
%
% Peri-stimulus time histograms of recorded spike times. 
% 
% Options:
% -chan [c1 c2 c3 ...] Set channels 
% -map  [c1 c2;  ... ] Set channels and arrangement of channels
% -unit [u1 u2 ... ]   Filter the displayed units 
% -rc   [rows cols]    Set subplot size (# rows / columns)
% 
% -merge               Compute PSTH for all spikes on the channel
% -no-hash             Only show spikes with unit code > 0. 
% -time [bin-times]    Bin times for histogram (default: 50 ms bin width)
%                       if len(times) <= 3 they're expanded out 
% -time [bin-width]    Set bin width (default: 50 ms) in ms                  
% -labels              Add channel labels (default: corners only)
% -no-labels           Suppress channel labels for corners
% -ticks               Toggle default ticks behaviour (ticks enabled for
%                       6 or fewer channels, disabled otherwise)
% -per-unit            One panel per unit (default: one per channel)
% -pass [pass_ids]     Select passes to plot
% 
% v0.1 - 4 June 2023 - Calvin Eiber <c.eiber@ieee.org>

named = @(s) strncmpi(s,varargin,numel(s));
get_ = @(v) varargin{find(named(v))+1};

if isfield(data,'config')
   this = @(d) plots.psth(d, varargin{:},'--epochs',data.epochs);
   tools.forWaveType(data, this, varargin{:});
   return
end

if any(named('-per-u'))
    [~,~,data.channel] = unique([data.channel data.unit],'rows');
end

if any(named('--epochs')), epochs = get_('--epochs');
    epochs.fs = mean((epochs.finish-epochs.start)./epochs.duration);
else error guestimate_epoch_structure
end

opts.epochs = epochs;

opts.n_smooth = 3;
opts.bins = 50; % default
opts.merge = any(named('-merge')); 

if any(named('-n-s')), opts.n_smooth = get_('-n-s'); end
if any(named('-tim')), opts.bins = get_('-tim'); end

if numel(opts.bins) == 1, 
    bin_w = opts.bins;
    opts.bins = [-(.5 : -1e3*epochs.frame_size(1)/epochs.fs/bin_w) ...
                   .5 :  1e3*epochs.frame_size(2)/epochs.fs/bin_w ];
    opts.bins = sort(opts.bins*bin_w,'ascend');

elseif numel(opts.bins) == 2, 
    bin_w = 5; 
    opts.bins = [-(.5 : -1e3*min(opts.bins)/epochs.fs/bin_w) ...
                   .5 :  1e3*max(opts.bins)/epochs.fs/bin_w ];
    opts.bins = sort(opts.bins*bin_w,'ascend');

elseif numel(opts.bins) == 3, 
    bin_w = opts.bins(3); 
    opts.bins = [-(.5 : -1e3*opts.bins(1)/epochs.fs/bin_w) ...
                   .5 :  1e3*opts.bins(2)/epochs.fs/bin_w ];
    opts.bins = sort(opts.bins*bin_w,'ascend');
end

tools.forChannels(data, @psth_plot, varargin{:}, ...
                   '--opts', opts, '--subplot')

ax = get(gcf,'Children');
pos = cat(1,ax.Position);
[~,blc] = min(pos*[1;1;0;0]);

xlabel(ax(blc),'time (s)'), 
ylabel(ax(blc),'response (imp/s)')

nP = numel(epochs.start);

arrayfun(@(a)fix_vertical_offset(a,nP),ax);

return


function psth_plot(data, index, opts)
        
bins = opts.bins;
dx = mean(diff(bins))/2;

nP = numel(data.pass_begin);

if index(2) > 0, color = lines(index(2)); color = color(end,:);
else color = [.6 .6 .6];
end

style = {color,'FaceAlpha',0.5,'EdgeColor',color};
if opts.merge, style{3} = 1; end

for pp = 1:nP

    t = data.time(data.pass == pp) - data.pass_begin(pp);    
    y = hist(t,bins/1e3); y = y/(2*dx)*1e3;

    if ~any(y), continue, end

    if opts.n_smooth > 0, 
        y = conv(y,ones(1,opts.n_smooth)/opts.n_smooth,'same'); 
    end

    x = [bins-dx;bins+dx]; x = x([1 1:end end])/1e3;
    y = [y;y]; y = [0; y(:); 0]; %#ok<AGROW> 

    fill(x,y,style{:},'userdata',[index pp])
end

return


function fix_vertical_offset(ax,n_passes)

h = findobj(ax,'type','patch');

indices = cat(1,h.UserData); 
pass_id = double(indices(:,3));

dy = arrayfun(@(x) max(x.YData), h); 
dy = 1.05 * max(dy); 

for ii = 1:numel(h)
    h(ii).YData = h(ii).YData - h(ii).YData(1) + dy*(pass_id(ii)-1);
end

ax.YLim = [-0.5 1.1*n_passes]*dy;

h = findobj(ax,'type','text');
if isempty(h), return, end
h(1).Position(2) = ax.YLim(2);

return