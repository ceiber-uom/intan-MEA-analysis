

function [epochs,opts] = segmentEpochs( trig, varargin )
% [epochs,opts] = segmentEpochs( trigger, [data], ... )
% 
% Options:
% 
% -opts [opts] : use specified options structure
% -min 0.4 : lower trigger level for stimulus offset
% -max 0.6 : upper trigger level for stimulus onset
% -before [0.5] : what fraction of the inter-stimulus interval is
%                 considered 'before' the next stimulus
%                 (as opposed to 'after' the previous)
% -jitter [0.01] : If variation in start times is greater than X, halt
%                  (error may indicate bad min/max values)
% -assess [0.90] : for determining input stimulus properties, use the
%                  middle X % of the recorded data
% 
% V0.1 - 27 September 2022 - Calvin Eiber <c.eiber@ieee.org>
% 

named = @(n) strncmpi(varargin,n,length(n));
get_ = @(v) varargin{find(named(v))+1};


if isstruct( trig )
  data = trig; 
  if isfield(trig,'trigger'), trig = trig.trigger; noun = '.trigger';
  elseif isfield(trig,'ADC'), trig = trig.ADC;     noun = '.ADC';
  end
else noun = ''; 
end

if isstruct( trig )
     assert( isfield(trig,'Data') && isfield(trig,'Time'), ...
         'expected data%s as a timeseries object', noun)
else assert( isa(trig,'timeseries'), ...
         'expected data%s as a timeseries object', noun)
end

%% Parse options

opts.threshold_min = 0.4;
opts.threshold_max = 0.6;
opts.before_fraction = 0.5; % for determining frame size
opts.allowable_jitter = 0.01; % coefficient of variation 
opts.assess_fraction = 0.9; % for getting properties

if any(named('-op')),  opts = get_('-op'); end
if any(named('-min')), opts.threshold_min = get_('-min'); end
if any(named('-max')), opts.threshold_max = get_('-max'); end
if any(named('-be')),  opts.before_fraction = get_('-be'); end
if any(named('-ji')),  opts.allowable_jitter = get_('-ji'); end
if any(named('-as')),  opts.assess_fraction = get_('-as'); end

%%

epochs.start     = [];
epochs.finish       = [];
epochs.frame_size = [];
epochs.duration  = [];
epochs.frequency = [];
epochs.amplitude = [];
epochs.phase     = []; 
epochs.average  = [];
epochs.units = struct; 

dd = sort([opts.threshold_min opts.threshold_max]); 

%%
stim_active = false; 
for tt = 1:numel(trig.Data)
  if stim_active
    if trig.Data(tt) > dd(1), continue, end
    epochs.finish(end+1,1) = tt;
    stim_active = false;
  else
    if trig.Data(tt) < dd(2), continue, end
    epochs.start(end+1,1) = tt;
    stim_active = true;
  end
end

if stim_active
  warning('mea:epochs:cutoff', ...
          'Recording ended while a stimulus was active. Data may be cut off.')
  epochs.finish(end+1) = tt; 
end
assert(numel(epochs.start) == numel(epochs.finish), 'start/end number of edges mismatch')

delta_t = diff(epochs.start); 
observed_jitter = std(delta_t) / mean(delta_t); 
if observed_jitter > opts.allowable_jitter

    clf
    plot(trig), hold on
    plot(xlim,[1 1]*dd(2),'--')
    plot(xlim,[1 1]*dd(1),'--')
    plot(trig.Time(epochs.start), ...
         dd(2) * ones(size(epochs.start)),'k^','markerFacecolor','k')

    plot(trig.Time(epochs.finish), ...
         dd(1) * ones(size(epochs.finish)),'kv','markerFacecolor','k')
    error('The observed variation in start timestamps (%0.2f%%) %s (%0.1f%%). %s', ...
          100*observed_jitter, 'exceeds the allowable jitter', ...
          100*opts.allowable_jitter, 'If this is OK, please increase this parameter')
end

%% Determine pre/post frame 

n_sam_active = mean(epochs.finish - epochs.start);
n_sam_gap = mean(diff(epochs.start)) - n_sam_active;
n_sam_pre = floor(opts.before_fraction * n_sam_gap);
n_sam_post = ceil((1-opts.before_fraction) * n_sam_gap);

epochs.frame_size = round([-n_sam_pre n_sam_active + n_sam_post]);

%% Determine stim properties

for ii = 1:numel(epochs.start)

    t = trig.Time(epochs.start(ii):epochs.finish(ii));
    y = trig.Data(epochs.start(ii):epochs.finish(ii));

    t = t - t(1); 
    epochs.duration(ii,1) = t(end);

    if opts.assess_fraction < 1
        n_cut = floor(numel(t) * (1-opts.assess_fraction)/2);
        t(1 : n_cut) = []; t(end-n_cut+1 : end) = []; 
        y(1 : n_cut) = []; y(end-n_cut+1 : end) = []; 
    end
    
    est = fit(t,y,'fourier1'); 

    epochs.frequency(ii,1) = est.w/2/pi;
    epochs.amplitude(ii,1) = sqrt( est.a1^2 + est.b1^2 );
    epochs.phase(ii,1)     = atan2( est.b1, est.a1 );
    epochs.average(ii,1)  = est.a0;
    
    if false
      %%
      cla %#ok<UNRCH> 
      plot(t,y), hold on
      plot(t, epochs.average(ii) + epochs.amplitude(ii) * ...
                    real(exp(2i*pi*epochs.frequency(ii)*t - ... 
                                1i*epochs.phase(ii))))
    end

    %%
    continue
    
    % I had a crack at trying to estimate the stimulus temporal frequency
    % from the spectrum directly ... invalid if the tf is not a clean
    % integer multiple of the sample rate.

    if isempty(fs), fs = 1./mean(diff(t)); end %#ok<UNRCH> 
    yy = fft(y); yy(ceil(end/2):end) = []; 
    hz = fs*(0:numel(yy)-1)/numel(t);
   [pk_amp,pk_idx] = max(abs(yy(2:end)));

end

epochs.units.start  = 'samples of ADC.Time';
epochs.units.finish = 'samples of ADC.Time';
epochs.units.frame_size = 'samples';
epochs.units.duration = trig.TimeInfo.Units;
epochs.units.frequency = 'Hz (if duration in s)';
epochs.units.amplitude = [trig.Name ' ' trig.DataInfo.Units];
epochs.units.phase = 'radians';
epochs.units.average = [trig.Name ' ' trig.DataInfo.Units];

%% 

if any(named('-d')), data = get_('-d'); end
if ~exist('data','var'), return, end

error todo_apply_epoch_segmentation_to_data


