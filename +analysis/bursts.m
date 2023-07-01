
function results = bursts(data, varargin)
% results = analysis.bursts( spike_data, ... )
% 
% For each channel and unit in the supplied data, this analysis looks at
%  the bursts (closely spaced spikes on the basis of ISI histograms). 
%  Results with the following data are returned for each channel / unit: 
% 
%  .channel_unit   : [channel unit] of analysis
%  .isi_threshold  : If two spikes occurred within this many seconds of one
%                    another, they were considered a burst for the purpose
%                    of this analysis
%  .mixture_model  : Gaussian fit object for burst vs non-burst (if auto) 
% 
%  .bursts : struct array of information about intra-burst firing rates: 
%  ..size          : data for bursts of X spikes (or 0 for non-bursting) 
%  ..count         : how many bursts were observed in the record? 
%  ..intra_isi     : stats about the ISI between spikes in the burst
%                    mean, median, standard deviation, and quartiles.
%  ..post_isi      : stats about the ISI between the last spike in the 
%                    burst and the next spike (mean, median, etc). 
%  .spike : struct with classification information for each spike
%  ..isi           : ISI for each spike (last spike gets median isi)
%  ..burst_id      : what burst is each spike a member of? (0 = non-burst)
% 
% This code can also be used to reproduce the plot for a results set: 
% >> analysis.multiState( results(1) ).
% 
% Spike-times may be shifted when re-plotting data 
% 
% Options: 
%  -pdf              : generate a PDF documenting the analysis
%  -cut []           : set explicit isi_threshold for burst detection
%  -lim [0.02 0.1]   : Set the upper bound of burst ISI center (s), 
%                              lower bound of non-burst ISIs (s)
%                      (defaults may adjust from these values if needed)
%  -replot           : Make plot from previously analysed result (default
%                       if result is the first input argument) 
% 
%  -min-dur [2 s]    : set minimum duration of a state occurance (in sec)
%  -smooth [4 s]     : set spikerate smoothing window width in sec
%                     (default: 2*min-dur seems to work well)
%                      this spikerate smoothing window is applied up to 4
%                      times in order to get a smooth enough histogram of
%                      instantaneous spike-rates to fit the Gaussian model.
%  -smooth-kern []   : explicitly set the smoothing kernel for spike-rates
%  -no-plot          : suppress generation plot for each channel / unit
% 
% (inherited from tools.forChannels)
%  -chan [c1 c2 ... ] : Select channels to analyse
%  -unit [u1 u2 ... ] : Filter for the specified units 
%  -roi [start end]   : set time analysis window (default: whole wave)
%  -merge-units       : Combined analysis of all spikes on each channel
%  -pass [pass_ids]   : Filter for epoch IDs to analyse
%  -no-hash           : Only analyse spikes with unit code > 0. 
% 
% v0.1 - 28 December 2022 - Calvin Eiber <c.eiber@ieee.org>

named = @(s) strncmpi(s,varargin,numel(s));

if nargin == 0, try data = evalin('caller','data'); end %#ok<TRYNC> 
elseif isfield(data,'isi_threshold') || any(named('-replot'))
    clf, make_singleChannel_plot(data(1)), return
end

do_PDF = any(named('-pdf'));
plots.PDF_tools('setup',do_PDF);


disp(datestr(now)), disp('Running spike burst analysis')

printInfo(); 
[~, results] = tools.forChannels(data, @run_burst_analysis, ...
                                        varargin{:}, '--ordered'); 
disp('Done! ')

% make_summary_graphic(results, varargin)

plots.PDF_tools('compile',do_PDF,'bursts-analysis (%d).pdf')

return



%% per-channel / per-unit analysis
function stats = run_burst_analysis(data, index, varargin)

printInfo('spike burst analysis (c%d.u%d) ', index);

named = @(s) strncmpi(s,varargin,numel(s));
get_ = @(v) varargin{find(named(v))+1};

%% Build variables to analyse
isi = diff(data.time);
isi(end+1) = median(isi);

n_spk = numel(data.time);
[density,log_x] = hist(log10(isi),sqrt(numel(isi))); %#ok<HIST> 
d_max = max(density);

if numel(density) < 6, stats = []; return, end

b_lims = [0.02 0.1]; % upper bound of burst ISI center (s)
                     % lower bound of non-burst ISIs (s)
b_lims(2) = min(b_lims(2), quantile(isi, 0.4));
b_lims(1) = min(b_lims(1), b_lims(2)/2); 

if any(named('-lim')), b_lims = get_('-lim'); end
if any(b_lims < 0), b_lims = 10.^b_lims; end
% test for multimodality with Hartigan's Dip Test

if any(named('-cut')), 
    burst_threshold = get_('-cut');
    gmm = []; 
else
    %%
    qi = range(quantile(log10(isi),([1 9])/10)); % isi quantiles
    lb = [d_max/50 -inf          qi/20 d_max/50 log10(b_lims(2)) qi/20 ];
    ub = [n_spk*2  log10(b_lims(1)) qi n_spk*2  inf              qi    ];
    gmm = fit(log_x', density', 'Gauss2', 'Upper', ub, 'Lower', lb);    
    gmm_pars = reshape(coeffvalues(gmm),3,[]); 
    
    % gmm = fitgmdist( v_(spike_rate), n_classes );
    g1c = @(p,x) p(1)*exp(-((x-p(2))./p(3)).^2) ; % p(1,2,3) -> a,b,c

    bursty = g1c(gmm_pars(:,1), log_x) > g1c(gmm_pars(:,2), log_x);
    idx = find(bursty ~= bursty(1), 1);

    if ~any(idx) % not a bursty cell
        % disp('no bursts detected'),  printInfo(); 
        stats = []; return
    end

    burst_threshold = 10.^mean(log_x(idx-[0 1])); 

    clear qi ub lb gmm_pars bursty idx
end

%% Given a burst threshold, identify bursts

bb = 0; 
burst_id = zeros(size(isi));
burst_nk = []; 

for tt = 1:numel(isi)-1
  if isi(tt) > burst_threshold, continue, end
  if burst_id(tt) == 0, bb = bb+1;  burst_nk(bb) = 1; end %#ok<AGROW> 
  burst_id(tt+[0 1]) = bb;
  burst_nk(bb) = burst_nk(bb) + 1; %#ok<AGROW> 
end

[bsf,bsx] = hist(burst_nk,1:max(burst_nk)); %#ok<HIST> 
bsf(1) = sum(burst_id == 0); 

spike_nib = burst_id; 
spike_nib(spike_nib>0) = burst_nk(spike_nib(spike_nib>0));

B = []; 

for bb = [0 2:numel(bsx)] % statistics about each burst size category

    this = struct;
    this.size = bb; % if range, give range
    this.count = bsf(max(1,bb));
    if this.count == 0, continue, end

    sel = (spike_nib == bb);
    if bb < 2
         sel = sel | [sel(2:end); 1];
    else sel = sel & [sel(2:end); 0];
    end

    ss{1} = sel; % intra_isi
    ss{2} = (spike_nib == bb) & ~sel; % post_isi
    var = {'intra_isi','post_isi'};

    for ii = 1:numel(var)
        this.(var{ii}).median =   median(isi(ss{ii}));
        this.(var{ii}).mean   =     mean(isi(ss{ii}));
        this.(var{ii}).iqr    = quantile(isi(ss{ii}),[.25 .75]);
        this.(var{ii}).stdev  =      std(isi(ss{ii}));
    end

    if isempty(B), B = this;
    else B(end+1) = this; %#ok<AGROW> 
    end
end

%% Format output

stats = struct; 

stats.channel_unit  = index; 
stats.isi_threshold = burst_threshold;
stats.bursts        = B; 
stats.mixture_model = gmm; 

stats.spike.isi      = isi; 
stats.spike.burst_id = burst_id;

%% Visualisation of ISI stats

if any(named('-no-plot')), return, end

make_singleChannel_plot(stats, data, spike_nib)

do_PDF = any(named('-pdf'));
plots.PDF_tools(gcf, do_PDF, 'page-%04d-%04d.ps', index)

return



function make_singleChannel_plot(stats, data, spike_nib)

clf
style = {'LineWidth',1.2};
font  = {'FontSize',6,'Color'};

bt = stats.isi_threshold; 

if nargin == 1
    data.time = cumsum(stats.spike.isi([end 1:end-1])); 
end

%%

C = lines(7); 
subplot(2,3,[1 2])
semilogy(data.time, stats.spike.isi,'.','Color',[.7 .7 .7]), hold on
plot(data.time([1 end]), [1 1]*bt,'-',style{:})
text(max(xlim),bt,sprintf('%0.5f',bt),font{:},C(2,:))

plots.tidy, yl = ylim; xl = xlim;
xlabel('time (s)'), ylabel('Burst ISI (s)')

% Recompute this
[density,log_x] = hist(log10(stats.spike.isi), ...
                        sqrt(numel(stats.spike.isi))); %#ok<HIST> 

subplot(2,3,3), cla
barh(log_x, density,1,'FaceColor',[.7 .7 .7],'EdgeColor','none')

hold on
if ~isempty(stats.mixture_model)
  plot(stats.mixture_model(log_x), log_x, style{:})
end

plot(xlim, [1 1] .* log10(bt), style{:},'Color',C(2,:))
plots.tidy, ylim(log10(yl)), set(gca,'YColor','none')
xlabel('spike count')

% title(hartigan's dip test p = %g)

if nargin < 3, 
    spike_nib = stats.spike.burst_id; 
    burst_nk = arrayfun(@(x) sum(spike_nib == x), 1:max(spike_nib));
    spike_nib(spike_nib>0) = burst_nk(spike_nib(spike_nib>0));
end

subplot(2,3,[4 5])
plot(data.time, spike_nib, '.', 'Color', [.7 .7 .7])
axis tight, plots.tidy, xlabel('time (s)'), ylabel('Burst n\_pikes')
bl = ylim; xlim(xl)


B = stats.bursts; 

for ii = 1:numel(B)
    text(max(xlim),B(ii).size,num2str(B(ii).count),font{:},[.4 .4 .4])
end

subplot(2,3,6)
si = [B.intra_isi]; delta = abs(cat(1,si.iqr) - cat(1,si.median));
errorbar([si.median], [B.size], [], [], delta(:,1), delta(:,2), ...
         '-o',style{:},'Color',C(2,:),'MarkerFaceColor',C(2,:));
hold on
si = [B.post_isi]; delta = abs(cat(1,si.iqr) - cat(1,si.median));
errorbar([si.median], [B.size], [], [], delta(:,1), delta(:,2), ...
         '-s',style{:},'Color',C(1,:),'MarkerFaceColor','w');

if numel(B) > 1
    text(B(2).intra_isi.median,1,'intra-burst',font{:},C(2,:),'horiz','center')
    text(B(2).post_isi.median,1,'post-burst',font{:},C(1,:),'horiz','center')
end
set(gca,'XScale','log'), plots.tidy, set(gca,'YColor','none')
plot([1 1]*bt,ylim,'-','Color',[0 0 0 0.3])
xlabel('ISI (s)'), ylim(bl);

suptitle(sprintf('channel %d unit %d', stats.channel_unit))
pause(0.01)

return



%% Visualuse summary results (not yet implemented)
function make_summary_graphic(results, varargin)

named = @(s) strncmpi(s,varargin,numel(s));
get_ = @(v) varargin{find(named(v))+1};

error TODO_make_graphic