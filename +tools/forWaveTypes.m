

function [data, varargout] = forWaveTypes( data, fun, varargin )
% [data, ...] = tools.forWaveTypes( data, function, [wave_types], ... )
% Apply FUNCTION to each specified wave_type in DATA. This was a common
% code pattern to enable top-level function calls on the intanData object.
% 
% types = {'AMP','AUX','VOLT','ADC','DI','DO','TEMP'}; (if not specified)
% 
% v0.1 - 10 October 2022 - CE

named = @(s) strncmpi(s,varargin,numel(s));
types = {'AMP','AUX','VOLT','ADC','DI','DO','TEMP'};

% implement the [types] in the input args
if nargin > 2 && iscell(varargin{1}), types = varargin{1}; end

varargout = cell(1,nargout-1);
if ~isfield(data,'config'), return, end

out = varargout; 
varargout = cellfun(@(x) struct, varargout,'unif',0);

do_type = cellfun(@(t) any(named(t)), types); 
if ~any(do_type), do_type(1) = true; end

for ty = types(do_type)
  if ~isfield(data,ty{1}), continue, end
  [data.(ty{1}),out{:}] = fun(data.(ty{1}));

  if sum(do_type) == 1, varargout = out;
  else
    for nn = 1:numel(out)
      varargout{nn}.(ty{1}) = out{nn};
    end
  end
end

  


