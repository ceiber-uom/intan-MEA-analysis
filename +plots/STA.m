function results = STA(data, varargin)
% results = plots.STA( spike_data, [trigger], ... )
% 
% Run a spike-triggered average analysis to look at correlations
% between a selected spike  (the 'trigger' spike) and the other
% channels and units in the supplied data. The following data 
% structure is returned for each channel / unit : 
% 
%  .channel_unit   : [channel unit] of analysis
%  .is_trigger     : was this channel/unit the trigger channel/unit?
%  .n_spikes       : how many spikes were in the source data?
%  .baseline       : average firing rate in impulses / second
%  .mean_sta       : average firing rate around trigger spike
%  .time           : time vector for mean_sta. Spikes which occur on the 
%                    observation channel just after a spike in the trigger 
%                    channel contribute to STA density at positive times. 
% 
% By default, a heatmap plot is generated with channels on the Y axis and
% lags on the X axis which shows the cross-correlation (STA) between the
% trigger channel and each other channel in imp/s. The channels are ordered
% roughly by shape similarity; a tick on the Y axis shows the trigger
% channel
% 
% This heatmap is interactive: click to get a measurement of the STA for a
% particular time and channel (output to the MATLAB console). 
% 
% Options: 
%  -pdf              : Generate a PDF of this figure for each channel/unit
%  -trig [chan unit] : set trigger channel/unit 
%                      (default: closest to average)
%  -trig-max         : trigger off the channel/unit with the most spikes
%  -win [t0:t1]      : set spike correlation view window. 
%                      The default window is  -20 : 20 ms (n=101). 
%                      The window can also be set by:
%                      -win [t0] (symm), -win [t0 t1], or -win [t0 t1 n] 
%  -no-plot          : suppress output plot generation
%  -replot           : re-plot the output results (replace data) possibly
%                      with different options 
% 
% Visualisation options: 
%  -raw             : show output without scaling to +- baseline
%  -hist            : show stack-of-histograms output 
%                     (only legible if the number of channels is small)
%  -no-pca          : show output image in channel order
%                     (default: order by similarity)
%  -pc [1]          : use PCA #n as measure of similarity
%  -units-s         : show plot in X-axis units of s (default: ms)
%  -label           : toggle channel labels
% 
% (inherited from tools.forChannels)
% -chan [c1 c2 ... ] : Select channels to analyse
% -unit [u1 u2 ... ] : Filter for the specified units 
% -roi [start end]  : set time analysis window (default: whole wave)
% -merge-units       : Combined analysis of all spikes on each channel
% -pass [pass_ids]     Filter for epoch IDs to analyse
% -no-hash             Only analyse spikes with unit code > 0. 
% 
% v0.1 - 31 December 2022 - Calvin Eiber <c.eiber@ieee.org>

if nargin == 0, try data = evalin('caller','data'); end, end %#ok<TRYNC> 

named = @(s) strncmpi(s,varargin,numel(s));
get_ = @(v) varargin{find(named(v))+1};

if any(named('-replot')) || isfield(data,'is_trigger')
    make_summary_graphic(data, varargin{:});
    return
end

if any(named('-pdf')), run_PDF_script(data, varargin{:}); return, end

disp(datestr(now)), disp('Running spike-triggered average analysis')

if any(named('-trig-max')), chan_unit = []; 
elseif any(named('-trig')), chan_unit = get_('-trig');
elseif nargin > 1 && isnumeric(varargin{1}), chan_unit = varargin{1}; 
else chan_unit = []; 
end

tfc_opt = {'--ordered'}; % For "get spike counts"
if any(named('-merge')) && ~any(named('-merge-sta')), 
    tfc_opt = [{'-merge'} tfc_opt]; 
end

if numel(chan_unit) <= 1 && ~isstruct(chan_unit) % auto detect
    if ~isempty(chan_unit), tfc_opt = [tfc_opt {'-chan'} chan_unit]; end
    [chan_unit, spike_count] = tools.forChannels(data, ...
                                   @(d,varargin) length(d.time), ...
                                                 tfc_opt{:}); 
    if any(named('-trig-max')), [~,idx] = max(spike_count); 
    else [~,idx] = min(abs(spike_count-mean(spike_count))); 
    end
    chan_unit = chan_unit(idx,:); 
    fprintf('Selected Trigger: c%d.u%d (%d spikes)\n', ...
                            chan_unit, spike_count(idx))
end

if isstruct(chan_unit), T = chan_unit;
  if ~isfield(T,'channel_unit'), 
    T.channel_unit = [T.channel(1) T.unit(1)]; 
  end
else
  [~,T] = tools.forChannels(data, @(d,varargin) d, tfc_opt{1}, ...
                            '-chan',chan_unit(1),'-unit',chan_unit(2));
  T.channel_unit = chan_unit; 
end

printInfo(); 
[~, results] = tools.forChannels(data, @run_sta, T, varargin{:}); 
disp('Done! ')

seq = cat(1,results.channel_unit);
[~,seq] = sort(seq(:,1) + (seq(:,2)/max(seq(:,2)+1)));
results = results(seq);

make_summary_graphic(results, varargin{:})

return

%% per-channel / per-unit analysis
function stats = run_sta(data, index, trig, varargin)

printInfo('STA analysis (c%d.u%d) ~ trigger: c%d.u%d ', index, trig.channel_unit);

named = @(s) strncmpi(s,varargin,numel(s));
get_ = @(v) varargin{find(named(v))+1};

window = 20;
if any(named('-win')),  window = get_('-win');  end
if ~any(named('units-s')), window = window / 1e3; end

if numel(window) <= 3
  switch(numel(window)),
    case 1, window = {-window(1), window(1), 102};
    case 2, window = { window(1), window(2), 102};
    case 3, window = { window(1), window(2), window(3)};
  end
  window = linspace(window{:});
  time = conv(window,[1 1]/2,'valid');
else
    time = window; 
    window = conv(window,[1 1]/2,'valid');
    window = [2*window(1)-window(2) window ...
              2*window(end) - window(end-1)];
end

recording_time = [0 inf]; 
if any(named('-roi')), recording_time = get_('-roi'); end % match .forChannels
recording_time = estimate_ROI(data, recording_time);
dt = mean(diff(time)); 

sr_hist = [];
nT = numel(trig.time);

for tt = 1:nT

    % positive number for delta: channel spike happens after trigger spike
    delta = data.time-trig.time(tt);
    sel = (delta >= window(1) & delta <= window(end)); 
    dy = hist(delta(sel),time); %#ok<HIST> 

    if isempty(sr_hist), sr_hist = zeros(size(dy)); end
    sr_hist = sr_hist + dy./nT./dt;  
    
    % see if trimming data.time speeds things up?
    % data.time(1:find(sel,1)-1) = []; 

    % if epoch data is supplied, a 'control' histogram should also be constructed 
    % by shifting spikes randomly among passes
end

% sr_hist is in units of impulses / second * (n_windows/n_windows)
% the expected hist base rate is (n_spikes/overall window)*total window width
baseline = numel(data.unit) / range(recording_time); 

%% Construct output data structure

stats = struct;
stats.channel_unit = index;
stats.is_trigger = all(index == trig.channel_unit);
stats.n_spikes   = numel(data.time);
stats.baseline   = baseline;
stats.mean_sta   = sr_hist;
stats.time       = time;


%% make empirical ROI window based on observed spike-times 
function time_ROI = estimate_ROI(data, time_ROI)

if nargin == 1, time_ROI = [0 inf]; end
if all(isfinite(time_ROI)), return, end

%% make ROI window based on observed spike-times 
derived_ROI = [min(data.time) max(data.time)]; 
padding = abs(time_ROI - derived_ROI);
padding(~isfinite(padding)) = []; 
if isempty(padding), padding = 0; 
else padding = mean(padding); 
end

if ~isfinite(time_ROI(1)), time_ROI(1) = derived_ROI(1) - padding; end
if ~isfinite(time_ROI(2)), time_ROI(2) = derived_ROI(2) + padding; end


%% Make summary visualisation
function make_summary_graphic(results, varargin)

named = @(s) strncmpi(s,varargin,numel(s));
get_ = @(v) varargin{find(named(v))+1};

if any(named('-no-plot')), return, end

clf
smooth = @(x,n) conv2(x,ones(1,n)/n,'same');
C = lines(7);

x = results(1).time;
if ~any(named('-units-s')), x = x*1e3; 
     time_unit = 'ms';
else time_unit = 's';
end


if any(named('-hist'))
  x = x([1 1:end end]);
  for ii = 1:numel(results)

    if results(ii).is_trigger && ~any(named('-keep-autocorr'))
      results(ii).mean_sta(results(ii).time == 0) = 0;
    end

    y = smooth(results(ii).mean_sta,5);
    b = results(ii).baseline / max(y(:)) * 1.2 + ii;
    y = y/max(y(:))*1.2 + ii; 

    if results(ii).is_trigger, c = C(3,:);
    else c = [.7 .7 .7];
    end

    fill(x, [ii y ii], c, 'EdgeColor','none','FaceAlpha',0.8),  hold on
    plot(x([1 end]), [b b           ],'-','Color',[c 0.3])
  end
else

  img = cat(1,results.mean_sta);
  if ~any(named('-raw'))
    img = img ./ cat(1,results.baseline);
  end

  if ~any(named('-no-pca'))
       [~,weight] = pca(img);

       pc_id = 1;
       if any(named('-pc')), pc_id = get_('-pc'); end
       [~,seq] = sort(weight(:,pc_id));
       img = img(seq,:);
  else seq = 1:size(img,1);
  end

  nnz = mean(img~=0,2);
  nnz_threshold = 0.25; 
  ok = (nnz>=nnz_threshold) | [results.is_trigger]';

  img = img(ok,:);
  seq = seq(ok);

    % img = smooth(img,5);

  h = imagesc(x, 1:size(img,1), img,'userdata',seq);
  caxis([0 quantile(img(:),0.99)])
  hold on

  if any(named('-label'))
    xt = xlim * [-0.02 -0.04; 1.02 1.04];
    cu = cat(1,results(seq).channel_unit);
    for ii = 1:numel(seq)
      text(xt(mod(ii,2)+1), ii, sprintf('c%d.%d',cu(ii,:)),'FontSize',6)
    end
  end

  if ~any(named('-no-in'))
    set(h,'ButtonDownFcn',@(a,b)write_channel(a,b,results(seq),time_unit));
  end

end

axis xy tight, plots.tidy
xlabel(['time, ' time_unit]), 
set(gca,'YTick',[]), 
ylabel('channels')

sel = [results.is_trigger];

if any(sel), 
     tcu = results(sel).channel_unit;
     title(sprintf('STA with trigger = c%d.u%d', tcu)) 
     set(gca,'YTick',find(sel(seq)),'YTickLabel','')
else title('STA')
end



%%

return

%% UI function
function write_channel(hobj,event,results,ms)

ip = round(event.IntersectionPoint(2)); 
[~,ix] = min(abs(hobj.XData-event.IntersectionPoint(1)));

fprintf('c%d.%d @ %0.2f %s: %0.4f\n', results(ip).channel_unit, ...
         hobj.XData(ix), ms, hobj.CData(ip,ix))

return

%% Script: Loop over this analysis to generate PDF for each channel/unit
function run_PDF_script(data, varargin)

named = @(s) strncmpi(s,varargin,numel(s));
varargin(named('-pdf')) = []; % prevent recursion

[~,T] = tools.forChannels(data, @(d,varargin) d, varargin{:},'--ordered');

plots.PDF_tools('setup'); 

for ii = 1:numel(T)
    plots.STA(data, '-trigger', T(ii), varargin{:});
    text(xlim*[.2;.8],ylim*[-0.02;1.02], ...
         sprintf('%d trigger spikes',numel(T(ii).time)),...
         'FontSize',7,'Color',[.4 .4 .4])
    plots.PDF_tools(gcf, 'page-%04d.ps', ii)
end

plots.PDF_tools('compile','spike-triggered averages (%d).pdf')





