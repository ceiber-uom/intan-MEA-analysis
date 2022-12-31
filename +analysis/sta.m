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
%  -pdf              : Generate a PDF of this analysis for each channel/unit
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

if numel(chan_unit) <= 1 % auto detect
    if ~isempty(chan_unit), tfc_opt = [tfc_opt {'-chan'} chan_unit]; end
    [chan_unit, spike_count] = tools.forChannels(data, ...
                                   @(d,varargin) length(d.time), ...
                                                 tfc_opt{:}); 
    if any(named('-trig-max')), [~,idx] = max(spike_count); 
    else [~,idx] = min(abs(spike_count-mean(spike_count))); 
    end
    chan_unit = chan_unit(idx,:); 
    fprintf('Selected Trigger: c%d.u%d (%d spikes)\n', chan_unit, spike_count(idx))
end

if isstruct(chan_unit), T = chan_unit;
else
    [~,T] = tools.forChannels(data, @(d,varargin) d, tfc_opt{1}, ...
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

printInfo('STA analysis (c%d.u%d) ~ trigger: c%d.u%d ', index, trig.channel_unit);

named = @(s) strncmpi(s,varargin,numel(s));
get_ = @(v) varargin{find(named(v))+1};

window = 0.1;
if any(named('-win')),  window = get_('-win'); end

if numel(window) <= 3
  switch(numel(window)),
    case 1, window = {-window(1), window(1), 102};
    case 2, window = { window(1), window(2), 102};
    case 3, window = { window(1), window(2), window(3)};
  end
  window = linspace(window{:});
  window = conv(window,[1 1]/2,'valid');
end

recording_time = [0 inf]; 
if any(named('-roi')), recording_time = get_('-roi'); end % match .forChannels
recording_time = estimate_ROI(data, recording_time);

dt = mean(diff(window)); 

sr_hist = [];
nT = numel(trig.time);

for tt = 1:nT

    % positive number for delta: channel spike happens after trigger spike
    delta = data.time-trig.time(tt);
    sel = (delta >= window(1) & delta <= window(end)); 
    dy = hist(delta(sel),window); %#ok<HIST> 

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
stats.time       = window;


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
% get_ = @(v) varargin{find(named(v))+1};

if any(named('-no-plot')), return, end

clf
smooth = @(x,n) conv2(x,ones(1,n)/n,'same');
C = lines(7);


if any(named('-hist'))
    
    for ii = 1:numel(results)
    
        x = results(ii).time([1 1:end end]);
        y = smooth(results(ii).mean_sta,5);
        b = results(ii).baseline / max(y(:)) * 1.2 + ii;
        y = y/max(y(:))*1.2 + ii; 
    
        if results(ii).is_trigger, c = C(3,:);
        else c = [.7 .7 .7];
        end
    
        fill(x, [ii y ii], c, 'EdgeColor','none','FaceAlpha',0.8),  hold on
        plot(x([1 end]), [b b],'-','Color',[c 0.3])
    end
else

    x = results(1).time;
    img = cat(1,results.mean_sta);
    if ~any(named('-raw'))
      img = img ./ cat(1,results.baseline);
    end

    if ~any(named('-no-pca'))
         [~,weight] = pca(img);
         [~,seq] = sort(weight(:,1));
         img = img(seq,:);
    else seq = 1:size(img,1);
    end

    nnz = mean(img~=0,2);
    nnz_threshold = 0.25; 

    img = img(nnz>=nnz_threshold,:);
    seq = seq(nnz>=nnz_threshold);

    % img = smooth(img,5);

    imagesc(x, 1:size(img,1), img,'userdata',seq)
    caxis([0 quantile(img(:),0.99)])
    hold on

end

axis xy tight, xl = xlim; plots.tidy
xlabel('time, s'), set(gca,'YTick',[]), ylabel('channels')

tcu = results([results.is_trigger]).channel_unit;
title(sprintf('STA with trigger = c%d.u%d', tcu)) 

%%

return


%% Script: Loop over this analysis to generate PDF for each channel/unit
function run_PDF_script(data, varargin)

named = @(s) strncmpi(s,varargin,numel(s));
varargin(named('-pdf')) = []; % prevent recursion

[~,T] = tools.forChannels(data, @(d,varargin) d, varargin{:},'--ordered');

plots.PDF_tools('setup'); 

for ii = 1:numel(T)
    analysis.sta(data, '-trigger', T(ii), varargin{:})
    plots.PDF_tools(gcf, 'page-%04d.ps', ii)
end

plots.PDF_tools('compile','spike-triggered averages (%d).pdf')





