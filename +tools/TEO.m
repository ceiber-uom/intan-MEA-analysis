
function [data,opts] = TEO(data, varargin)
% [epochs,opts] = tools.TEO( [data], [trigger], ... )
% 
% Implement the Teager Energy Operator (TEO). The TEO is a function that 
%   accepts a time-series Y(t) and returns a new timeseries:  
%   Y(t)' = Y(t)*Y(t) - Y(t-w)*Y(t+w)
% 
% The TEO is a nonlinear operation which emphasises high-frequency
%   transients in the data, and so is good for spike-detection. I found
%   that the TEO can improve spike detection for low SNR data and provide
%   estimates of the output SNR (based on the expected peak height spectra 
%   of spikes vs. noise detections), as detailed in chapter <5> of <thesis>
% 
% This implementation is adapted from the implementation used in 
%   <my thesis>, and originates from <ref>. It has been adapted to 
%   natively process the data returned by tools.readIntan or 
%   tools.segmentEpochs.  
% 
% Options: 
% -smooth [1]   : set wave smoothing factor (default: no smoothing). 
% -kernel [vec] : set smoothing kernel (default: pascal kernel). 
% -width [1]    : set window width for TEO
% -raw          : if set, return in units-squared - default is to return 
%                 sqrt(|y'(t)|) which is in the same units as y. 
% 
% Note: smoothing is applied /after/ TEO computation. 

if isstruct(data) && isfield(data,'config')
    this = @(d) tools.TeagerEnergyOperator(d, varargin{:});
   [data, info] = tools.forWaveType(data, this, varargin{:});
    data.config.TEO = info;
    return
end

named = @(s) strncmpi(s,varargin,numel(s));
get_ = @(v) varargin{find(named(v))+1};

sf = 1; % smoothing factor (obj.opts.TEO_smoothing)
dt = 1; % window width (obj.opts.TEO_size)

if any(named('-sm')), sf = get_('-sm'); end
if any(named('-ke')), window = get_('-ke'); % (obj.opts.TEO_kernel)
else window = diag(rot90(pascal(sf))); window = window/sum(window);
end
if any(named('-wi')), dt = get_('-wi'); end

% remind again what pascal window is? 

disp('Applying Teager energy operator')

teo_ = @(y) y.^2 - circshift(y,[dt 0 0]) .* circshift(y,[-dt 0 0]);

no_unpack = isnumeric(data);

if no_unpack, wave = teo_(data); else wave = teo_(data.data); end
if length(window) > 1, wave = convn(wave,window,'same');      end
if ~any(named('-raw')), wave = sqrt(abs(double(wave)));       end
if no_unpack, data = wave; else data.data = wave;             end

opts.TEO_width = dt;
opts.TEO_smoothing = numel(window);
if any(named('-ke')), opts.TEO_custom_kernel = window; end
if opts.TEO_smoothing == 1, opts.TEO_smoothing = false; end
opts.return_raw_TEO = any(named('-raw'));


return

tic %#ok<UNRCH> 

% Old code for extra-huge files: 
if size(wave,3) >= 1600 % for extra-large files, break into chunks
  wave = double(wave);
  for cc = 1:size(wave,3)
    wave(:,:,cc) = teo_(wave(:,:,cc));
    if length(window) > 1 
        wave(:,:,cc) = convn(wave(:,:,cc),window,'same'); 
    end
    wave(:,:,cc) = sqrt(abs(double(wave(:,:,cc))));
    if toc > 5, disp(['TEO: chan ' num2str(cc)]), tic; end
  end
end
