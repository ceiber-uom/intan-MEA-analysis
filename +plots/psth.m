
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
% -no-bar              Disable stimulus bar on PSTH plots
% -ticks               Toggle default ticks behaviour (ticks enabled for
%                       6 or fewer channels, disabled otherwise)
% -per-unit            One panel per unit (default: one per channel)
% -per-pass            One row per pass (default if epochs.condition_id not
%                       set)
% -sem                 Show SEM around average PSTH
%                       (only relevent if -per-pass not enabled)
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

opts.n_smooth = 3; % n-point moving average smooth for psth display
opts.bins = 50;    % << set default bin size in ms
opts.merge = any(named('-merge')); 
opts.trial_avg = isfield(epochs,'condition_id') && ~any(named('-per-p')); 
opts.stim_bar = ~any(named('-no-b'));
opts.show_sem = any(named('-sem'));

if any(named('-n-s')), opts.n_smooth = get_('-n-s'); end
if any(named('-tim')), opts.bins = get_('-tim'); end

if any(named('-roi')), time_roi = get_('-roi');
  if numel(time_roi) == 1, time_roi = [0 time_roi]; end
else time_roi = epochs.frame_size/epochs.fs;
end

if numel(opts.bins) == 1, 
    bin_w = opts.bins;
    opts.bins = [-(.5 : -1e3*time_roi(1)/bin_w) ...
                   .5 :  1e3*time_roi(2)/bin_w ];
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

oa = findobj(0,'type','axes'); % axes before plotting

tools.forChannels(data, @psth_plot, varargin{:}, ...
                   '--opts', opts, '--subplot')

ax = setdiff(find(0,'type','axes'), oa);
% ax = get(gcf,'Children');

if numel(unique([ax.Parent])) > 1, 
     blc = 1:numel(ax); 
else pos = cat(1,ax.Position);
  [~,blc] = min(pos*[1;1;0;0]);
end

xlabel(ax(blc),'time (s)'),         % TODO - test for not -single-fig
ylabel(ax(blc),'response (imp/s)')

if opts.trial_avg, 
     nY = numel(unique(epochs.condition_id));
else nY = numel(epochs.start);
end

arrayfun(@(a) fix_vertical_offset(a,nY,opts),ax);

if isfield(epochs,'block_id'), 
    error TODO_implement_block_logic
end

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
if opts.trial_avg
    [conds,~,cond_id] = unique(opts.epochs.condition_id);
    nC = numel(conds); 
    y_sum = zeros(nP,1) * bins;
end

for pp = 1:nP

    t = data.time(data.pass == pp) - data.pass_begin(pp);    
    y = hist(t,bins/1e3); y = y/(2*dx)*1e3;

    if ~any(y), continue, end

    if opts.trial_avg, y_sum(pp,:) = y; continue, end

    if opts.n_smooth > 0, 
        y = conv(y,ones(1,opts.n_smooth)/opts.n_smooth,'same'); 
    end

    x = [bins-dx;bins+dx]; x = x([1 1:end end])/1e3;
    y = [y;y]; y = [0; y(:); 0]; %#ok<AGROW>
    fill(x,y,style{:},'userdata',[index pp])

    if opts.stim_bar
      dur = opts.epochs.duration(pp); 
      rectangle('position',[0 -max(y)*0.03 dur max(y)*0.03], ...
                'facecolor',[.3 .3 .3],'userdata',[index pp], ...
                'clipping','off','edgecolor',[.3 .3 .3])
    end
end

if ~opts.trial_avg, return, end

for cc = 1:nC

    sel = (cond_id == cc);

    y = mean(y_sum(sel,:),1); 
    if opts.n_smooth > 0, 
        y = conv(y,ones(1,opts.n_smooth)/opts.n_smooth,'same'); 
    end

    if opts.show_sem

        se = std(y_sum(sel,:),[],1) ./ sqrt(sum(sel));
        x = bins/1e3;

        if opts.n_smooth > 0, 
          se = conv(se,ones(1,opts.n_smooth)/opts.n_smooth,'same'); 
        end
        
        style{5} = 'none';
        fill([x fliplr(x)],[y-se fliplr(y+se)], ...
               style{:},'userdata',[index cc])
        plot(x,y,'-','color',color,'linewidth',1.5,'userdata',[index cc])



    else
        x = [bins-dx;bins+dx]; x = x([1 1:end end])/1e3;
        y = [y;y]; y = [0; y(:); 0]; %#ok<AGROW>
        fill(x,y,style{:},'userdata',[index cc])
    end

    if opts.stim_bar
      dur = opts.epochs.duration(find(sel,1)); 
      rectangle('position',[0 -max(y)*0.03 dur max(y)*0.03], ...
               'facecolor',[.3 .3 .3],'userdata',[index cc], ...
               'clipping','off','edgecolor',[.3 .3 .3])
    end
end

return


function fix_vertical_offset(ax,n_passes,opts) %#ok<INUSD> 

graphics_object = {'patch','line','rectangle'};
dy = []; 

for g = 1:3
    
  h = findobj(ax,'type',graphics_object{g});    
  if isempty(h), continue, end
        
  indices = cat(1,h.UserData); 
  pass_id = double(indices(:,3));    

  if g < 3 && isempty(dy)
    dy = arrayfun(@(x) max(x.YData), h); 
    dy = 1.05 * max(dy); 
    ax.YLim = [-0.5 1.1*n_passes]*dy;
  end
  if g == 3
    rect_h = cat(1,h.Position);
    rect_h = max(rect_h(:,4)); 
  end
           
  for ii = 1:numel(h)
    if g < 3, % patch or line
      h(ii).YData = h(ii).YData + dy*(pass_id(ii)-1);
    else % rectangle
      h(ii).Position(2) = dy*(pass_id(ii)-1) - rect_h;
      h(ii).Position(4) = rect_h;
    end
  end
end

h = findobj(ax,'type','text');
if isempty(h), return, end
h(1).Position(2) = ax.YLim(2);

return