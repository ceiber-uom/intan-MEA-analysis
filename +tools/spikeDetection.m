
% function data = spikeDetection( data, varargin )

% TODO implement this from my old Ph.D. code.

function data = spikeDetection(data)





% classdef SpikeDetect

%  properties opts data, end
%  properties(Hidden) indices rawEpochs filename filepath, end
% 
%  methods

 % Constructor
 function obj = SpikeDetect(varargin)
     obj.opts.algorithm = @threshold_RMS;
     obj.opts.spikeWindow = ((1:32)-8);  % 
     obj.opts.threshold = 3.5;
     obj.opts.minInterval = 0.001;
     obj.opts.preprocesing = {@preprocesser_CommonMode};
     obj.opts.saveWaves = 'raw';
     obj.opts.commonMode = true;
     obj.opts.display = false;
     obj = applyVargs2Opts (obj, varargin);
 end

function obj = Detect(obj)
    if isempty(obj.data), obj = obj.Load(); end
    if ~isfield(obj.data.settings, 'samplingRate'),
     obj.data.settings.samplingRate = obj.data.settings.General.samplingRate;
    end

 if isempty(obj.rawEpochs), epochs = obj.data.epochs;
      obj.rawEpochs = epochs;
 else epochs = obj.rawEpochs;
 end
 obj.data = rmfield(obj.data,'epochs');
 % Execute extensible preprocessing steps
 for pi = 1:length(obj.opts.preprocesing)
 [epochs, obj] = obj.opts.preprocesing{pi}(epochs, obj);
 end
 t1 = tic;
 obj.indices = cell(size(epochs,2), size(epochs,3));
 if isfield(obj.opts, 'Channels'), C_range = obj.opts.Channels;
 else C_range = 1:size(epochs,3);
 end

 for ci = C_range
 if max(max(epochs(:,:,ci))) == 0, continue, end
 for ei = 1:size(epochs, 2)
 obj.indices{ei,ci} = obj.opts.algorithm(epochs(:,ei,ci), ci, ...
 obj.data, ...
 obj.opts);
 end
 if toc>15, disp(['Formatting Channel ' num2str(ci)]), tic, end
 end
 obj.data.epochs = epochs(:,:,C_range);
 end

 function obj = Save(obj, varargin)
 obj.opts.spikeWindow = obj.opts.spikeWindow(:); % convert to column
 [spikeCounts,spikeTimes,spikeShapes] = Output(obj, varargin);

 % get variables to save
 params = obj.data.params;
 rmsValues = obj.data.rmsValues;
 settings = obj.data.settings;
 tags = obj.data.tags;

 settings.SpikeDetect = cleanStruct(obj.opts);

 disp(['Saving ' strrep(obj.filename,'-epochs','-spikes') '...' ]);
 save([obj.filepath strrep(obj.filename,'-epochs','-spikes')], ...
 'spikeShapes', 'spikeCounts', 'spikeTimes', ...
 'params', 'settings', 'rmsValues', 'tags');
 disp(['Results saved to file ' strrep(obj.filename,'-epochs','-spikes')]);
 end
 
 function t = times(obj, varargin)
 % --- 14 lines of code hidden ---
 % Code is same as “Output” except only spikeTimes are returned
 end

 function s = shapes(obj, varargin)
 % --- 25 lines of code hidden ---
 % Code is same as “Output” except only spikeShapes are returned
 end

 function [spikeCounts,spikeTimes,spikeShapes] = Output(obj, varargin)

 obj.opts.spikeWindow = obj.opts.spikeWindow(:); % convert to column vector
 spikeTimes = cell(size(obj.data.epochs,2), size(obj.data.epochs,3));
 spikeCounts = int32(zeros(size(obj.data.epochs,2), size(obj.data.epochs,3)));
 spikeShapes = cell(size(obj.data.epochs,2), size(obj.data.epochs,3));

 if strncmpi(obj.opts.saveWaves, 'raw', 3), epochs = obj.rawEpochs;
 else epochs = obj.data.epochs;
 end

 tic
 for ci = 1:size(obj.data.epochs, 3)
 for ei = 1:size(obj.data.epochs, 2)
 if ~isempty(obj.indices{ei,ci})
 temp = obj.opts.spikeWindow * ones(size(obj.indices{ei,ci})) + ...
 ones(size(obj.opts.spikeWindow)) * (obj.indices{ei,ci});
 temp(temp < 1) = 1;
 temp(temp > size(obj.data.epochs,1)) = size(obj.data.epochs,1);

 spikeCounts(ei, ci) = length(obj.indices{ei,ci});
 spikeTimes{ei,ci} = obj.indices{ei,ci} / obj.data.settings.samplingRate;
 spikeShapes{ei, ci} = reshape(epochs(temp,ei,ci),size(temp));

 end
 end
 if toc>15, disp(['Formatting Channel ' num2str(ci)]), tic, end
 end
 end
 function obj = Load(obj, varargin)

 if length(varargin) >= 1 && exist(varargin{1},'file')

 [obj.filename, obj.filepath] = strtok(fliplr(varargin{1}),'\');
 obj.filename = fliplr(obj.filename);
 obj.filepath = fliplr(obj.filepath);
 else
 [obj.filename, obj.filepath] = uigetfile({'*epochs.mat'});
 end
 disp(['Loading ' obj.filename])
 obj.data = load([obj.filepath obj.filename]);
 end
 
 % property access methods
 function re = raw(obj), re = obj.rawEpochs; end
 function fn = name(obj), fn = obj.filename; end
 function fp = path(obj), fp = obj.filepath; end
 function obj = rename(obj, fullname)
   [fn, fp] = strtok(fliplr(fullname),'\');
   obj.filename = fliplr(fn);
   obj.filepath = fliplr(fp);
 end

% Threshold crossing detector variants:
function indices = threshold_RMS(wave, chan, data, opts)
 if length(opts.threshold) >= chan, T = opts.threshold(chan);
 else T = opts.threshold;
 end
 T = T * data.rmsValues(chan);
 W = opts.minInterval * data.settings.samplingRate;
 indices = peakseek(abs(wave), W, T);

    
    
 function indices = threshold_ABS(wave, chan, data, opts)
 if length(opts.threshold) >= chan, T = opts.threshold(chan);
 else T = opts.threshold;
 end
 W = opts.minInterval * data.settings.General.samplingRate;
 indices = peakseek(abs(wave), W, T);

%% Preprocessers (passed as function handles):
% function [epochs, obj] = preprocess_ ... ( epochs, obj )

    


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



function [epochs, obj] = preprocesser_AnsariBradly(epochs, obj )



    
if ~isfield(obj.opts, 'AnsariBradlyTail'), 
     obj.opts.AnsariBradlyTail = 'Right'; end
end
if ~isfield(obj.opts, 'AnsariBradlyXform'),
     obj.opts.AnsariBradlyXform = @(p) sqrt(-log10(p));
end

 if ~isfield(obj.opts, 'targetSpikerate'), obj.opts.targetSpikerate = 10; end % spikes / s
 targetSpikes = obj.opts.targetSpikerate*numel(epochs(:,:,1))/obj.data.settings.samplingRate;

 if max(obj.opts.display_autoThreshold), figure('Color','w'), end
 disp(['Calculating Thresholds to achieve ' num2str(round(targetSpikes)) ' spikes/channel.'])
 t1 = tic;

 ed = obj.data.settings.SaveEpochs.epochDuration;
 if ~isfield(obj.opts, 'referenceROI'), obj.opts.referenceROI = [ed-.085 ed]; end
 roi = obj.opts.referenceROI/ed*size(epochs,1); % window of interest

 ref = rand(round([1 targetSpikes]))*size(epochs,1);
 ref = ref(ref > min(roi) & ref < max(roi));
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
 autoThresholdPlot(obj, peaks, T(ci), targetSpikes, size(epochs), ci, P(ci));
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



