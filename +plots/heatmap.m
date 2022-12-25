
function heatmap (data, varargin)
% function plots.heatmap( data, [Y], ... )
%
% Basic heatmap of recorded spike counts (or other channel data)
% 
% Options:
% -y [y1 y2 ... ]      Set data to visualuse 
% -chan [c1 c2 c3 ...] Set channels
% -map  [c1 c2;  ... ] Set channels and arrangement of channels
% -rc   [rows cols]    Set subplot size (# rows / columns)
% -roi  [t0 t1]        Set beginning and end of plotting window
% -pass [pass_ids]     Select passes to plot (for expoched data)
% -unit [u1 u2 ... ]   Filter the displayed units 
% -labels              Add channel labels (default: corners only)
% -no-labels           Suppress channel labels for corners
% -scale               Toggle colorscale (default: enabled if subplot 111)
% -count               Plot spike count (default: mean spikerate)
%                      (only if Y data is not otherwise specified)
% 
% v0.1 - 25 December 2022 - Calvin Eiber <c.eiber@ieee.org>



named = @(s) strncmpi(s,varargin,numel(s));
get_ = @(v) varargin{find(named(v))+1};

if isfield(data,'config')
   this = @(d) plots.heatmap(d, varargin{:});
   tools.forWaveType(data, this, varargin{:});
   return
end

channel_map = plots.layout(data, varargin{:});
nC = max(channel_map(:));

if any(named('-y')), y_data = get_('-y');
elseif nargin > 1 && isnumeric(varargin{1}), y_data = varargin{1}; 
else
    %% Select spikes to include by applying filters
    
    if any(named('-roi')), include = get('-roi'); 
        if numel(include) == 2
            time_ROI = [min(include) max(include)]; 
            include = (data.time >= time_ROI(1) & ...
                       data.time <= time_ROI(2));
        end
    else include = true(size(data.time)); 
         time_ROI = [0 inf]; 
    end
    
    if any(named('-unit')), unit_ok = get('-unit'); 
      if isnumeric(unit_ok), unit_ok = ismember(data.unit, unit_ok); end
      if any(unit_ok), include = include & unit_ok; end
    end
    
    if any(named('-pass')), pass_ok = get_('-pass');
    elseif ~isfield(data,'pass')
         data.pass = ones(size(data.time));
         pass_ok = true;
    else pass_ok = true(size(data.data,3),1);
    end
    if ~islogical(pass_ok)
        pass_ok = ismember(1:max(data.pass), pass_ok); 
    end
    if any(pass_ok), include = include & pass_ok; end

    dci = data.channel(include);

    y_data = zeros(nC,1);
    for ii = 1:numel(dci)
        y_data(dci(ii)) = y_data(dci(ii)) + 1;
    end
    
    if ~any(named('-count'))
      if any(~isfinite(time_ROI))
        %% make ROI window based on observed spike-times 
        derived_ROI = [min(data.time(include)) ...
                       max(data.time(include))]; 
        padding = abs(time_ROI - derived_ROI);
        padding(~isfinite(padding)) = []; 
        if isempty(padding), padding = 0; 
        else padding = mean(padding); 
        end

        if ~isfinite(time_ROI(1))
            time_ROI(1) = derived_ROI(1) - padding; 
        end
        if ~isfinite(time_ROI(2))
            time_ROI(2) = derived_ROI(2) + padding; 
        end
      end

      y_data = y_data / diff(time_ROI); % to spikes/second
    end
end

if numel(y_data) ~= nC
  error implement_procedure_for_sparse_data
end

%%

img = y_data(channel_map); 

cla reset, imagesc(img')
axis image ij off

p = get(gca,'Position'); 
do_colorbar = (prod(p(3:4)) > 0.5) ~= any(named('-scale')); 
if do_colorbar, colorbar; end

% Add channel labels
label_all = any(named('-label'));
label_none = any(named('-no-l'));

for cc = 1:nC
  if label_none, break, end
  [a,b] = find(channel_map == cc);
  ok = label_all || ...
         (a == 1 || a == size(channel_map,1)) && ...
         (b == 1 || b == size(channel_map,2));
  if ~ok, continue, end
  text(a,b,num2str(cc),'Color',[.9 .9 .9],'FontSize',14,'horiz','center')
end

return
