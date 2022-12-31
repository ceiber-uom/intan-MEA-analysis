function results = sta(data, varargin)
% results = analysis.sta( spike_data, [trigger], ... )
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
% Options: 
%  -trig [chan unit] : select trigger channel/unit (default: closest to average)
%  -trig-max         : trigger off the channel/unit with the highest spikerate
%  -win [t0:t1]      : set spike correlation view window. The default window is
%                      -0.1 : 0.1 seconds (n=101). The window can also be set by:
%                      -win [t0] (symmetric), -win [t0 t1], or -win [t0 t1 n] 
%  -no-plot          : suppress output plot generation
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

disp(datestr(now)), disp('Running spike-triggered average analysis')

named = @(s) strncmpi(s,varargin,numel(s));
get_ = @(v) varargin{find(named(v))+1};

if any(named('-trig-max')), chan_unit = []; 
elseif any(named('-trig')), chan_unit = get_('-trig')
elseif nargin > 1 & isnumeric(varargin{1}), chan_unit = varargin{1}; 
else chan_unit = []; 

tfc_opt = {'--ordered'} % For "get spike counts"
if any(named('-merge')) && ~any(named('-merge-sta')), tfc_opt = [{'-merge'} tfc_opt]; end

if numel(chan_unit) <= 1 % auto detect
    if ~isempty(chan_unit), tfc_opt = [tfc_opt {'-chan'} chan_unit]; end
    [chan_unit, spike_count] = tools.forChannels(data, @(d,varargin) length(d.time), tfc_opt{:}); 
    if any(named('-trig-max')), [~,idx] = max(spike_count); 
    else [~,idx] = min(abs(spike_count-mean(spike_count))); 
    end
    chan_unit = chan_unit(idx,:); 
end

if isstruct(chan_unit), T = chan_unit
else
    [~,T] = tools.forChannels(data, @(d,varargin) d, tfc_opt{1},
                            '-chan',chan_unit(1),'-unit',chan_unit(2));
    T.channel_unit = chan_unit; 
end

printInfo(); 
[~, results] = tools.forChannels(data, @run_sta, T, varargin{:}); 
disp('Done! ')

make_summary_graphic(results, varargin{:})

return

%% per-channel / per-unit analysis
function stats = run_sta(data, index, trig, varargin)

printInfo('f012 analysis (c%d.u%d) ', index, trig.channel_unit);

named = @(s) strncmpi(s,varargin,numel(s));
get_ = @(v) varargin{find(named(v))+1};

window = [0.1]
if any(named('-win')), window = get_('-win'); end
if numel(window) == 1,     window = linspace(-window(1), window(1), 101)
elseif numel(window) == 2, window = linspace( window(1), window(2), 101)
elseif numel(window) == 3, window = linspace( window(1), window(2), window(3))
end

recording_time = [0 inf]; 
if any(named('-roi')), recording_time = get_('-roi'); end % match .forChannels
recording_time = estimate_ROI(data, recording_time);

dt = mean(diff(window)); 

sr_hist = [];
nT = numel(trig.time);

for tt = 1:nT

    % positive number for delta: channel spike happens after trigger spike
    delta = data.time-trig.time(tt)
    sel = (delta >= window(1) & delta <= window(end)); 
    dy = hist(delta(sel),window); 

    sr_hist = sr_hist + dy/nT/dt;  
    
    % see if trimming data.time speeds things up?
    % data.time(1:find(sel,1)-1) = []; 

    % if epoch data is supplied, a 'control' histogram should also be constructed 
    % by shifting spikes randomly among passes

end

% sr_hist is in units of impulses / second * (n_windows/n_windows)
% the expected hist base rate is (n_spikes/overall window)*total window width
baseline = numel(data.unit) / range(recording_time); 

%% Construct output data structure

stats = struct
stats.channel_unit = index
stats.is_trigger = all(index == trig.channel_unit)
stats.n_spikes   = numel(data.time)
stats.baseline   = baseline
stats.mean_sta   = sr_hist
stats.time       = window


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
smooth = @(x,n) conv(x,ones(1,n)/n,'same');
C = lines(7);

for ii = 1:numel(results)

    x = results(ii).time([1 1:end end]);
    y = smooth(results(ii).mean_sta,5);
    b = results(ii).baseline / max(y(:)) * 1.2 + ii
    y = y/max(y(:))*1.2 + ii; 

    if results(ii).is_trigger, c = C(3,:)
    else c = [.7 .7 .7];
    end

    fill(x, [ii y ii], c, 'EdgeColor','none','FaceAlpha',0.8),  hold on
    plot(x([1 end]), [b b],'-','Color',[c 0.3])
end

axis tight, xl = xlim; plots.tidy, xlim(xl); 
xlabel('time, s'), title('firing rate')
set(gca,'YColor','none')

%% Second subplot - show just the state-classification for each chan/unit

chan_unit = cat(1,results.channel_unit);

xd = xl * [1.02; -0.02];

for cc = 1:max(chan_unit(:,1))
    sel = find(chan_unit(:,1) == cc);
    if ~any(sel), continue, end
  
    text(2*xd, mean(sel), sprintf('c%d', cc),'Color',[.2 .2 .2],'FontSize',6,'Horiz','right')
    plot([xd xd], [min(sel) max(sel)], 'Color', [.2 .2 .2], 'LineWidth', 1.1,'clipping','off')
end

tcu = chan_unit([results.is_trigger],:);
title(sprintf('STA with trigger = c%d.u%d', tcu)) 

return
