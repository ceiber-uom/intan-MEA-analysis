
function data = TeagerEnergyOperator(data, varargin)

if isfield(data,'config')
    this = @(d) tools.TeagerEnergyOperator(d, varargin{:});
   [data, info] = tools.forWaveTypes(data, this, varargin{:});
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

disp('Applying Teager energy operator')

teo_ = @(y) y.^2 - circshift(y,[dt 0 0]) .* circshift(y,[-dt 0 0]);

wave = teo_(data.data);
if length(window) > 1, wave = convn(wave,window,'same'); end
wave = sqrt(abs(double(wave)));

data.data = wave; 
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
