
function results = states(data, varargin)
% results = analysis.states( spike_data, ... )
% 
% For each channel and unit in the supplied data, this analysis looks at
%  the time-varying firing rate of that unit and attempts to identify
%  periods of high and low firing rate by fitting a simple gaussian
%  mixture model. A results structure with the following fields is 
%  returned for each channel / unit : 
% 
%  .channel_unit   : [channel unit] of analysis
%  .n_classes      : how many classes were fit to this unit?
%  .bin_size       : what was the time resolution of the firing rate? 
%                   (default: pick resolution to give mean 1 spike/bin)
%  .wave  : information about firing rate and state: 
%  ..time          : time vector 
%  ..rate          : firing rates at each moment in time
%  ..class         : classification of firing rate (1 = low, 2 = high)
%  ..start_indices : index of times where classification changes
%  ..duration      : duration of each state occurance (in sec)
% 
%  .mixture_model  : Gaussian fit object relating object to class ID 
%  .mixture_params : amplitude, center, and width (a,b,c) for each Gaussian
%                    [3 x n_classes matrix] 
%  .spikerate      : stats about the spikerate in each state.
%                    mean, median, standard deviation, and quartiles.
%                    each statistic is a [1 x n_classes] matrix; IQR is
%                    returned as a [2 x n_classes]: 25%, 75% points. 
%  .duration       : stats about the duration of occurances of each state. 
%                    the format is the same as .spikerate
%  .time_fraction  : what fraction of time is spent in each state?
% 
% Options: 
%  -pdf              : generate a PDF documenting the analysis
%  -roi [start end]  : set time analysis window (default: whole wave)
%  -bins [width]     : set fixed bin width for spikerate 
%                     (default: auto, 1 spike/bin on average)
%  -nc [2]           : set number of classes to attempt to fit
%  -min-dur [2 s]    : set minimum duration of a state occurance (in sec)
%  -smooth [4 s]     : set spikerate smoothing window width in sec
%                     (default: 2*min-dur seems to work well)
%                      this spikerate smoothing window is applied up to 4
%                      times in order to get a smooth enough histogram of
%                      instantaneous spike-rates to fit the Gaussian model.
%  -smooth-kern []   : explicitly set the smoothing kernel for spike-rates
%  -plot             : generate a plot for each channel / unit
%  -no-plot          : suppress generation of the summary plots showing the
%                      combined results 
% (inherited from tools.forChannels)
% -chan [c1 c2 ... ] : Select channels to analyse
% -unit [u1 u2 ... ] : Filter for the specified units 
% -merge-units       : Combined analysis of all spikes on each channel
% -pass [pass_ids]     Filter for epoch IDs to analyse
% -no-hash             Only analyse spikes with unit code > 0. 
% 
% v0.1 - 28 December 2022 - Calvin Eiber <c.eiber@ieee.org>

named = @(s) strncmpi(s,varargin,numel(s));

if nargin == 0, try data = evalin('caller','data'); end, end %#ok<TRYNC> 

disp(datestr(now)), disp('Running spike-rate multi-state analysis')

do_PDF = any(named('-pdf'));
plots.PDF_tools('setup',do_PDF);

printInfo(); 
[~, results] = tools.forChannels(data, @run_msa, varargin{:},'--ordered'); 
disp('Done! ')

make_summary_graphic(results, varargin)
plots.PDF_tools(gcf, do_PDF, 'page-%04d-%04d.ps', 0, 1)
plots.PDF_tools('compile',do_PDF,'states-analysis (%d).pdf')


return

%% per-channel / per-unit analysis
function stats = run_msa(data, index, varargin)

printInfo('multi-state analysis (c%d.u%d) ', index);

named = @(s) strncmpi(s,varargin,numel(s));
get_ = @(v) varargin{find(named(v))+1};

%% Build variables to analyse

time = [0 inf]; 
if any(named('-roi')), time = get_('-roi'); end % match .forChannels
time = estimate_ROI(data, time);

if any(named('-bin')), bin_size = get_('-bin');
else bin_size = diff(time) / numel(data.time);
end

n_classes = 2;
if any(named('-nc')), n_classes = get_('-nc'); end

minimum_dur = 2; % s
if any(named('-min')), minimum_dur = get_('-min'); end
minimum_bins = minimum_dur / bin_size;

%%

time = time(1):bin_size:time(2); 
spike_rate = hist(data.time,time); %#ok<HIST> 
spike_rate = spike_rate / bin_size;

if any(named('-smooth-k')), sk = get_('-smooth-k');
elseif any(named('-sm')), sk = ones(1,ceil(get_('-sm')/bin_size));  
else sk = ones(1,ceil(2*minimum_bins)); % empirically derived default
end, sk = sk / sum(sk); 

% apply conv filter up to 4 times
for iter = 1:4
    if isempty(sk), break, end
    spike_rate = conv(spike_rate,sk,'same');
    if numel(unique(spike_rate)) < 3*(n_classes+1), break, end
end

n_hist = min(32, numel(unique(spike_rate)));
if any(named('-nh')), n_hist = get_('-nh'); end

[rate_count,inst_rates] = hist(spike_rate, n_hist); %#ok<HIST> 

v_ = @(x) reshape(x,[],1); 

%% Fit Gaussian Mixture model

if numel(rate_count) < 3*n_classes, n_classes = 1; end 
if numel(rate_count) < 3, stats = []; return, end

lb = [0; 0; 0] * ones(1,n_classes); 
ub = [2*max(rate_count); max(inst_rates); max(inst_rates)] * ones(1,n_classes);

gmm = fit(v_(inst_rates), v_(rate_count), ... 
          sprintf('gauss%d', n_classes),  ... 
          'Lower', v_(lb), 'Upper', v_(ub));
gmm_pars = reshape(coeffvalues(gmm),3,[]); 

% gmm = fitgmdist( v_(spike_rate), n_classes );
g1c = @(p,x) p(1)*exp(-((x-p(2))./p(3)).^2) ; % p(1,2,3) -> a,b,c

est_state = zeros(size(spike_rate)); 
est_prob  = zeros(size(spike_rate)); 

for cc = 1:n_classes
    this_prob = g1c( gmm_pars(:,cc), spike_rate ); 
    sel = (this_prob > est_prob); 
    est_state(sel) = cc; 
    est_prob(sel) = this_prob(sel); 
end

clear this_prob cc sel est_prob lb ub sk

% for a "schmidtt trigger" style approach: 
% inflection_points = gmm_pars(2,:) + [-1; 1]*gmm_pars(3,:); 
% Not 100% obvious how to set this up for n > 2 classes, 

%% Minimum duration in state approach

change_idx = find([0 est_state] ~= [est_state 0]);
duration = diff(change_idx)+1; change_idx(end) = []; 

state = est_state;

[~,seq] = sort(duration,'ascend');

for ss = 1:numel(seq), ii = seq(ss);
  if duration(ii) >= minimum_bins, continue, end
  if duration(ii) == 0, continue, end

  idx = change_idx(ii) + (0:duration(ii));
  % if idx(1) == 1, continue, end
  if ii == 1, continue, end

  previous = state(idx(1)-1);
  state(idx) = previous;
  
  duration(ii-1) = duration(ii-1) + duration(ii); 
  if idx(end) >= numel(est_state), continue, end
  if state(idx(end)+1) == previous
    duration(ii-1) = duration(ii-1) + duration(ii+1);
    duration(ii+1) = 0;
  end
end

state(numel(est_state)+1 : end) = []; % trim ? 

change_idx = find([0 state] ~= [state 0]);
duration = diff(change_idx)+1; change_idx(end) = []; 
new_state = state(change_idx); 
duration = duration * bin_size;

the = @(f,varargin) arrayfun(f,1:n_classes,varargin{:});
msr = the(@(c) median(spike_rate(state == c)));

[~,seq] = sort(msr); 
state = seq(state); 
gmm_pars = gmm_pars(:,seq);

%% Format output

stats = struct; 

stats.channel_unit = index; 
stats.n_classes = n_classes;
stats.bin_size = bin_size; 

stats.wave.time          = time; 
stats.wave.rate          = spike_rate;
stats.wave.class         = state; 
stats.wave.start_indices = change_idx;
stats.wave.duration      = duration;

stats.mixture_model      = gmm; 
stats.mixture_params     = gmm_pars;

for var = {'spikerate','duration'}
  switch var{1}
    case 'spikerate', y = @(c) spike_rate(state == c);   
    case 'duration',  y = @(c) duration(new_state == c);
  end

  stats.(var{1}).median =  the(@(c) median(y(c)));
  stats.(var{1}).mean   =  the(@(c)   mean(y(c),'omitnan'));
  stats.(var{1}).iqr    = [the(@(c) quantile(y(c),0.25));
                           the(@(c) quantile(y(c),0.75))];
  stats.(var{1}).stdev  =  the(@(c)      std(y(c),'omitnan')); 

end

stats.time_fraction   = the(@(c) mean(state == c));

%% If requested, generate plot

do_PDF = any(named('-pdf'));
if ~any(named('-plot')) && ~do_PDF, return, end

%%

% sbw = gausswin(ceil(3/bin_size))
smooth = @(x,n) conv(x,ones(1,n)/n,'same');

clf
subplot(1,4,[1 3])
plot(time, smooth(spike_rate,15),'LineWidth', 1.2,'Color','k')
plots.tidy, yl = ylim; 

hold on
% scatter(time, 0*time + 0.95*max(ylim), [], naive_class, '.')
try scatter(time, 0*time + 0.5*min(ylim), [], state, '.'), end %#ok<TRYNC> 
if n_classes > 1
   caxis([1 n_classes]), colormap(gca, lines(n_classes));
else caxis([0 1]), colormap([.5 .5 .5; .5 .5 .5])
end

subplot(1,4,4)
barh(inst_rates, rate_count, 1, 'FaceColor', [.2 .2 .2], ... 
                                'FaceAlpha', 0.4, ...
                                'EdgeColor','none')
plots.tidy, ylim(yl), set(gca,'YColor','none'), hold on
xlim(xlim)

irx = linspace(inst_rates(1), inst_rates(end), 101);
C = lines(n_classes); 

for ii = 1:n_classes
    plot(g1c(gmm_pars(:,ii),irx), irx, 'Color', C(ii,:),'LineWidth',1.2)
end

suptitle(sprintf('channel %d unit %d', stats.channel_unit))
pause(0.01)

do_PDF = any(named('-pdf'));
plots.PDF_tools(gcf, do_PDF, 'page-%04d-%04d.ps', index)

return

%% Results summary graphics
function make_summary_graphic(results, varargin)

named = @(s) strncmpi(s,varargin,numel(s));
get_ = @(v) varargin{find(named(v))+1};

if any(named('-no-plot')), return, end

R = [results.wave];
[R.index] = deal(0);
for ii = 1:numel(R)
    R(ii).index = ii+0*R(ii).time; 
end

%% First subplot - show firing rates for each channel / unit 

clf
subplot(1,2,1)
smooth = @(x,n) conv(x,ones(1,n)/n,'same');
C = lines(max([results.n_classes])); 

for ii = 1:numel(results)

    x = R(ii).time([1 1:end end]);
    y = smooth(R(ii).rate,15); 
    y = y/max(y(:))*1.2 + ii; 

    for cc = 1:results(ii).n_classes
        py = y;
        py(R(ii).class ~= cc) = ii; 
        fill(x, [ii py ii], C(cc,:),'EdgeColor','none','FaceAlpha',0.5 )
        hold on
    end
end

axis tight, xl = xlim; plots.tidy, xlim(xl); 
xlabel('time, s'), title('firing rate')
set(gca,'YColor','none')

%% Second subplot - show just the state-classification for each chan/unit

subplot(1,2,2), cla
scatter([R.time],[R.index],[],[R.class],'.')

axis tight, plots.tidy, xlim(xl), hold on
xlabel('time, s'), title('state')
set(gca,'YColor','none')

colormap( (lines(max([R.class]))+1)/2 );
chan_unit = cat(1,results.channel_unit);

xd = xl * [1.02; -0.02];

for sp = 1:2
  subplot(1,2,sp)
  for cc = 1:max(chan_unit(:,1))
    sel = find(chan_unit(:,1) == cc);
    if ~any(sel), continue, end
  
    text(2*xd, mean(sel), sprintf('c%d', cc),'Color',[.2 .2 .2],'FontSize',6,'Horiz','right')
    plot([xd xd], [min(sel) max(sel)], 'Color', [.2 .2 .2], 'LineWidth', 1.1,'clipping','off')
  end
end

% suptitle(sprintf('channel %d unit %d', chan_unit))









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