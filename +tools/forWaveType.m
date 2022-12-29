

function [data, varargout] = forWaveType( data, fun, varargin )
% [data, ...] = tools.forWaveTypes( data, function, [wave_types], ... )
% 
% Apply FUNCTION to a specified wave_type in DATA. This was a common code
%   pattern to enable top-level function calls on the intanData object.
%   By default, if the user did not specify one of the following, FUNCTION
%   is applied to data.AMP (if data.SPIKE is not present) or data.SPIKE 
% 
% types = .AMP, .AUX, .VOLT, .ADC, .DI, .DO, .TEMP, (.SPIKE)
% 
% Example Usage: 
% if isfield(data,'config')
%    this = @(d) tools.removeCommonMode(d, varargin{:});
%   [data, info] = tools.forWaveTypes(data, this, varargin{:});
%    data.config.commonMode = info;
%    return
% end
% 
% v0.1 - 10 October 2022 - CE

named = @(s) strncmpi(s,varargin,numel(s));
types = {'SPIKE','AMP','AUX','VOLT','ADC','DI','DO','TEMP'};

% implement the [types] in the input args
if nargin > 2 && iscell(varargin{1}), types = varargin{1}; end

varargout = cell(1,nargout-1);
if ~isfield(data,'config'), return, end

out = varargout; 
varargout = cellfun(@(x) struct, varargout,'unif',0);

do_type = cellfun(@(t) any(named(t)), types); 
if ~any(do_type), 
    do_type(find(1)) = true;
    do_type(1) = true; 
end % by default 

for ty = types(do_type)
  if ~isfield(data,ty{1}), continue, end

  if nargout == 0, fun(data.(ty{1}))
  else [data.(ty{1}),out{:}] = fun(data.(ty{1}));
  end
  
  if sum(do_type) == 1, varargout = out;
  else
    for nn = 1:numel(out)
      varargout{nn}.(ty{1}) = out{nn};
    end
  end
end
