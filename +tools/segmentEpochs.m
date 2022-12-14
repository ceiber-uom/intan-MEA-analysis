
function [epochs,opts] = segmentEpochs( data, varargin )
% [epochs,opts] = tools.segmentEpochs( [data], [trigger], ... )
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

verbose = ~any(named('-v')); 

if nargin == 0
  if evalin('caller','exist("intanData","var")')
    data = evalin('caller','intanData');
    fprintf('Please call tools.segmentEpochs( intanData )')
  else error('Not enough input arguments.')
  end
end

if any(named('-epochs')) % apply epoch filter to existing data
    data = apply_segmentation(data, get_('-epochs'));
end

noun = ''; 
if isa(data,'timeseries'), trig = data; data = []; 
elseif ~isstruct(data) 
    error('Expected an intanData object (as returned by tools.readIntan)')
else 
  if isfield(data,'trigger'), trig = data.trigger; noun = '.trigger';
  elseif isfield(data,'ADC'), trig = data.ADC;     noun = '.ADC';
  elseif isfield(data,'Data') && isfield(data,'Time') 
                              trig = data; data = []; 
  end
end

if verbose
    fprintf('segmenting epochs based on photocell trigger in data%s ... \n', noun)
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
for tt = 1:length(trig.time)
  if stim_active
    if trig.data(tt) > dd(1), continue, end
    epochs.finish(end+1,1) = tt;
    stim_active = false;
  else
    if trig.data(tt) < dd(2), continue, end
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

v_ = @(x) reshape(x,[],1); % vertical vector

for ii = 1:numel(epochs.start)

    t = trig.time(epochs.start(ii):epochs.finish(ii));
    y = trig.data(epochs.start(ii):epochs.finish(ii));

    t = t - t(1); 
    epochs.duration(ii,1) = t(end);

    if opts.assess_fraction < 1
        n_cut = floor(numel(t) * (1-opts.assess_fraction)/2);
        t(1 : n_cut) = []; t(end-n_cut+1 : end) = []; 
        y(1 : n_cut) = []; y(end-n_cut+1 : end) = []; 
    end
    
    est = fit(v_(t),v_(y),'fourier1'); 

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
epochs.units.amplitude = [trig.name ' ' trig.DataInfo.Units];
epochs.units.phase = 'radians';
epochs.units.average = [trig.name ' ' trig.DataInfo.Units];

%% 
if any(named('-d')), data = get_('-d'); 
elseif isempty(data) && nargin > 1 && isstruct(varargin{1})
    data = varargin{1};
end
if ~exist('data','var'), return, end

data.epochs = epochs;
data.config.epochs = opts; 
data = apply_segmentation(data, epochs);

% shift output argument sequence
opts = epochs; epochs = data;

return


function data = apply_segmentation(data, epochs)

data_fields = fieldnames(data); 
all_upper = @(s) all(ismember(s,['A':'Z' '_'])); 
data_fields = data_fields(cellfun(all_upper,data_fields));

for channel_type = reshape(data_fields, 1, [])
  %%
  source = data.(channel_type{1});
  is_ts = ~isstruct(source);

  if is_ts, this = timeseries;
       this.name = source.name;
       this.TimeInfo = source.TimeInfo;
       this.DataInfo = source.DataInfo;
  else this = source;
       this.time = [];
       this.data = []; 
       if size(source.data,3) > 1
         error('this data appears to have already been segmented into epochs')
       end
  end

  for tt = 1:numel(epochs.start)
    sel = epochs.start(tt) + (epochs.frame_size(1):epochs.frame_size(2));
    ok = (sel > 0 & sel < numel(source.time));
    sel(~ok) = 1; 

    if isempty(this.time)
        this.time = (sel - epochs.start(tt)) * mean(diff(source.time));
    end

    this.data = cat(3, this.data, source.data(sel,:));
    this.data(~ok,:,end) = nan;
  end
  %%
  data.(channel_type{1}) = this; 
end

if isfield(data,'spikes')
    error('TODO: implement epoch for spiketimes')
end

return
