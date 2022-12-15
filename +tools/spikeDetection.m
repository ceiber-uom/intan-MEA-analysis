
function [data,opts] = spikeDetection(data, varargin)
% function data = spikeDetection( data, varargin )
% 
% TODO reimplement this from my old Ph.D. code.
% 

if isfield(data,'config')
    this = @(d) tools.spikeDetection(d, varargin{:});
   [spikes, info] = tools.forWaveTypes(data, this, varargin{:});
    data.config.SpikeDetection = info;
    error here_append_spikes
    return
end

named = @(s) strncmpi(s,varargin,numel(s));
get_ = @(v) varargin{find(named(v))+1};

opts = parse_options(varargin{:}); 

if opts.detect_symmetric, epochs = abs(epochs); end

epochs = data.data; % [time x chan x ... ]

dt = mean(diff(data.time));

if isempty(opts.channel_threshold)
    error here_get_channel_RMS
end

epoch_size = size(epochs);
assert(epoch_size(1) > 1, ...
      'epochs should be a (time x channels x ... ) matrix')
epochs = reshape(epochs, size(epochs,1), []); 

if any(named('-ch')), channel_list = get_('-ch');
else channel_list = 1:size(epochs,2);
end

if numel(opts.channel_thresholds) == 1
    opts.channel_thresholds = repmat(opts.channel_thresholds, ...
                                     epoch_size(1), 1);
end

%% Core detection loop

spike_indices = cell(epoch_size(2:end)); 

tt = tic;

for dd = 1:numel(spike_indices) 

  cc = mod(dd-1,epoch_size(1))+1;
  if ~any(channel_list == cc), continue, end

  spk_cutoff = opts.channel_thresholds(cc,:); 
  spk_min_dt = opts.minimum_seperation / dt; 

  spike_indices{dd} = tools.peakseek(epochs(:,dd), spk_min_dt, spk_cutoff);

  if toc(tt)>5, tt = tic;  
    fprintf('Detecting Spikes, %d/%d\n ', dd, numel(spike_indices))
  end
end

%% Format output 

spike.times  = cellfun(@(x) x*dt, spike_indices, 'unif', 0);
spike.counts = cellfun(@numel, spike_indices);

spike.waveform = cell(size(spike_indices));

spk_window = reshape(opts.spike_wave_indices,1,[]); 

tt = tic;
if ~any(named('-no-w'))
  for dd = 1:numel(spike_indices) 
    if isempty(spike_indices{dd}), continue, end

    idx = max(min(1, spk_window + reshape(spike_indices{dd},1,[])), ...
                   size(epochs,1));
  
    spike.waveform{dd} = epochs(idx,dd); 
  
    if toc(tt)>5, tt = tic; 
      fprintf('Gathering Waveforms, %d/%d\n ', dd, numel(spike_indices))
    end
  end
end

return





function opts = get_detect_options(varargin) 

named = @(s) strncmpi(s,varargin,numel(s));
get_ = @(v) varargin{find(named(v))+1};

opts = struct;

if any(named('-th'))
    opts.threshold_algorithm = 'as specified';
    opts.channel_thresholds = get_('-th');
% elseif ...
else
    opts.threshold_algorithm = 'RMS';     
    opts.channel_thresholds = []; % auto from RMS algorithm
    opts.relative_RMS_threshold = 3.5; % default threshold for threshold_RMS
    if any(named('-rms')), opts.relative_RMS_threshold = get_('-rms'); end
end

opts.minimum_seperation = 0.001; % in units of data.time
opts.detect_symmetric = any(named('-abs')) || any(named('-sym')); 

if ~any(named('-no-w'))
  opts.spike_wave_indices = ((1:32)-8); % TODO - check this vs Plexon defaults
end

if any(named('-o')), oi = get_('-o'); 
  for f = fieldnames(oi)',
    if isfield(opts, f{1}), opts.(f{1}) = oi.(f{1}); end
  end
end

return


%% From below is OLD code <thesis-era>, my professional opinion is that I
% may have done a bunch of reinventing the wheel in a way that was not
% useful - basically that code structure was almost a self-contained
% load-run-and-save system that is way clearer as a script. 

% function [epochs, obj] = preprocesser_Wavelets(epochs, obj )
% 
%  if ~isfield(obj.opts,'CWT_wavelet'), obj.opts.CWT_wavelet = 'Bior 1.3'; end
%  if ~isfield(obj.opts,'CWT_scales'), obj.opts.CWT_scales = unique(round(logspace(0,1.5,20))); end
%  if ~isfield(obj.opts,'CWT_reject'), obj.opts.CWT_reject = 0; end
%  if ~isfield(obj.opts,'CWT_reduce'), obj.opts.CWT_reduce = @(x) geomean(abs(x)); end
%  disp('Applying Continuous Wavelet Transform')
%  tic
%  cwr = obj.opts.CWT_reject;
%  for c = 1:size(epochs,3)
%  for r = 1:size(epochs,2)
%  cwt_epoch = cwt(epochs(:,r,c),obj.opts.CWT_scales, obj.opts.CWT_wavelet);
%  if cwr > 0, cwt_epoch(abs(cwt_epoch) < cwr) = 0; end
%  epochs(:,r,c) = obj.opts.CWT_reduce(cwt_epoch);
%  end
%  disp(['CWT: chan ' num2str(c)]);
%  end
%  toc

function [AB_stats, opts] = spikeTimeDistStats(data, varargin )

opts = struct; 
opts.stats_roi = [0 inf];
opts.target_spikerate = 10; % imp/s

if any(named('-tail')), opts.AB_tail = get_('-tail'); end
if any(named('-transform')), opts.AB_transform = get_('-transform'); end
if any(named('-roi')), opts.AB_roi = get_('-roi'); end
if any(named('-r')), opts.target_spikerate = get_('-r'); 
elseif any(named('-tar')), opts.target_spikerate = get_('-tar');
end

if numel(opts.stats_roi) == 1, 
    opts.stats_roi = opts.stats_roi * [1 inf]; 
end

delta_t = mean(diff(data.time)); 

data_size = size(data.data); 
if numel(data_size) < 3, data_size(3) = 1; end

target_spikes = opts.target_spikerate * delta_t * data_size(1) * ...
                                             prod(data_size(3:end)); 

fprintf('Estimating thresholds to achieve %0.0f spikes/channel... \n', target_spikes)






t1 = tic;

roi = data.time >= opts.AB_roi(1) & data.time < opts.AB_roi(2); 

fake_spikes = rand(ceil([1 targetSpikes]))*size(data.data,1);
fake_spikes = fake_spikes(roi(ceil(fake_spikes)));

clear ed
obj.opts.algorithm = @threshold_ABS;



P = zeros(size(epochs,3),1);
T = zeros(size(epochs,3),1);

for ci = 1:size(epochs,3)

  obj.opts.threshold(ci) = 1;
  if max(max(epochs(:,:,ci))) == 0, continue, end
  wave = epochs(:,:,ci);
  wave = wave(:);

  indices = obj.opts.algorithm(wave, ci, obj.data, obj.opts);

  ROI_indices = mod(indices,size(epochs,1));
  inROI = (ROI_indices > min(roi) & ROI_indices < max(roi));

  while sum(inROI) < length(ref)/size(epochs,3)
    obj.opts.threshold(ci) = 0.1 * obj.opts.threshold(ci);
    indices = obj.opts.algorithm(wave, ci, obj.data, obj.opts);
    ROI_indices = mod(indices,size(epochs,1));
    inROI = (ROI_indices > min(roi) & ROI_indices < max(roi));
  end

  indices = indices(inROI);
  [peaks, sortIndex] = sort(abs(wave(indices)), 'descend');
  T(ci) = (peaks(round(targetSpikes))+peaks(round(targetSpikes) + 1))/2;

  indices = indices(sortIndex);
  indices = indices(1:round(targetSpikes));
  ROI_indices = mod(indices,size(epochs,1));

  % ref + median(ROI_indices)
  [~,P(ci)] = ansaribradley(ref, ROI_indices, 'Tail', obj.opts.AnsariBradlyTail);

  if max(obj.opts.display_autoThreshold == ci)
    % autoThresholdPlot(obj, peaks, T(ci), targetSpikes, size(epochs), ci, P(ci));
  end
  if toc(t1) > 5 || P(ci) < 0.05,
    fprintf('... C%d: p = %.3f at %d µV\n', ci, P(ci), T(ci))
    t1 = tic;
  end
end

T = T ./ obj.data.rmsValues';
obj.opts.algorithm = @threshold_RMS;
obj.opts.autoThresholdLevels = T;
obj.opts.autoThresholdPvals = P;

P = obj.opts.AnsariBradlyXform(P);
P(~isfinite(P)) = 0;
P(isnan(P)) = 0;
obj.opts.threshold = nansum(T.*P)/nansum(P);
disp(['p-Weighted threshold = ' num2str(obj.opts.threshold) 'xRMS'])

return