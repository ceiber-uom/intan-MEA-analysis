
function [data,info] = removeCommonMode(data, varargin)
% data = tools.removeCommonMode(data) 
% Common-mode removal algorithm - remove the (weighted) average signal
%   from all channels in the input signal: 
% - FAST mode (default) : the first PCA component is removed from the input
% - LINEST mode (-lin)  : the average signal is computed and scaled for
%                         each channel to minimise the mean abs deviation.
% 
% If the input data is an intanData structure, then by default this 
% operation is applied to the 'AMP' channels; other channels can be
% specified. For instance:
%  data = tools.removeCommonMode(data, 'AMP','AUX') removes common-mode
%  from the AMP and the AUX channels (seperately). 
% 
% see also https://doi.org/10.26190/unsworks/18556 (appendix 3)
% 
% v2.0 - 10 October 2022 - Calvin Eiber <c.eiber@ieee.org>

named = @(s) strncmpi(s,varargin,numel(s));

if isfield(data,'config')
   this = @(d) tools.removeCommonMode(d, varargin{:});
  [data, info] = tools.forWaveTypes(data, this, varargin{:});
   data.config.commonMode = info;
   return
end

wave_size = size(data.data);
wave = double(data.data);
wave_type = class(data.data);

% convert 'wave' to a channels x time matrix (flatten epochs)
% save shape info to convert back later

if length(wave_size) > 2 % epochs mode
     wave = reshape( permute(wave,[2 1 3]), wave_size(2), []);
elseif wave_size(1) == length(data.time)
     wave = wave';
else wave_size = []; 
end

if any(named('lin')) || any(named('-lin'))

    % this takes about 10x longer than 
    
    commonMode = mean(wave,1);
    weights = zeros(size(wave,1), 1);
    residual = zeros(size(wave,1), 1);

    for chan = 1:size(wave,1)
      if max(wave(chan,:)) == min(wave(chan,:)), continue, end
      fitFcn = @(A) mean(abs(wave(chan,:) - A*commonMode));
      [weights(chan), residual(chan)] = fminsearch(fitFcn, 1);
      printInfo('Channel %d processed.', chan)
    end

    rcm_method = 'fminsearch (linfit)';
else

    [weights,commonMode] = pca(wave','centered',false,'NumComponents',1);  
    commonMode = commonMode .* max(abs(weights));
    weights = weights .* (1./max(abs(weights)));
    commonMode = commonMode'; 

    rcm_method = 'pca';
end

wave = wave - weights * commonMode;
commonMode = cast(commonMode,wave_type);

if isempty(wave_size), data.data = cast(wave, wave_type); 
else 
    wave_size = num2cell(wave_size);
    wave = reshape( wave', wave_size{:} );
    data.data = cast(wave, wave_type); 
end

info.profile = commonMode;
info.weights = weights;

if exist('residual','var'), info.residual = residual; end

info.method = rcm_method; 
