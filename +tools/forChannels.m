
function [list, varargout] = forChannels (data, fun, varargin)
% [channel, ...] = tools.forChannels( data, function, ... )
% 
% Apply FUNCTION for each channel (and optionally unit, for sorted spikes) 
% in DATA. This is a common code pattern which has been factored out here.
% This function includes an optional call to tools.forWaveTypes and so can
% be called on the top-level intan data object as returned by
% tools.readIntan or tools.readPlexon. 
% 
% NOTE at this time this code is only implemented for spike data, not for
%      wave data. 
% 
% FUNCTION should have the following syntax: 
%  >> FUNCTION(data, [channel unit], varargin) 
% 
% This function provides the following options as a shared syntax: 
% 
% (inherited from plots.layout): 
% -chan [c1 c2 c3 ...] Set channels 
% -map  [c1 c2;  ... ] Set channels and arrangement of channels
% -rc   [rows cols]    Set subplot size (# rows / columns)
% 
% (implemented in this function):
% -roi  [t0 t1]        Set beginning and end of plotting window
% -unit [u1 u2 ... ]   Filter for the specified units 
% -merge-units         Ignore unit codes  
% -pass [pass_ids]     Select passes to plot (for expoched data)
% -no-hash             Only show spikes with unit code > 0. 
% 
% The following meta-directives are supported:
% --subplot            Generate a subplot for each call of FUNCTION. 
% -enable-empty-call   If set, FUNCTION is called even if there is no data
% 
% Exposed options for SUBPLOT mode:
% -ticks               Toggle Ticks display for subplots
% -no-tidy             Supress plots.tidy and basic subplot format
% -labels              Label each subplot (default: corners only)
% -no-labels           Suppress labeling of subplots, even corners. 
% --label-color [rgb]  Set label color, default [.4 .4 .4] (grey)
% 
% v0.1 - 23 December 2022 - Calvin Eiber <c.eiber@ieee.org>

named = @(s) strncmpi(s,varargin,numel(s));
get_ = @(v) varargin{find(named(v))+1};

if isfield(data,'config')
   varargout = cell(1,nargout-1);
   this = @(d) tools.forChannels(d, fun, varargin{:});
   if nargout == 0, tools.forWaveType(data, this, varargin{:},'--get-types');
   else [list,varargout{:}] = tools.forWaveType(data, this, varargin{:},'--get-types');
   end
   if nargout > 0, list = list.(types{1}); end %#ok<USENS> 
   return
end

channel_map = plots.layout(data, varargin{:});
subplot_nxy = fliplr(num2cell(size(channel_map)));

opts.ticks = (numel(channel_map) > 6) == any(named('-tick'));
opts.hash  = ~any(named('-no-hash'));
opts.merge = any(named('-ignore-u')) || any(named('-merge'));

% Add channel labels
opts.label_all    =  any(named('-label'));
opts.label_none   =  any(named('-no-l'));
opts.skip_empty   = ~any(named('-enable-empty'));
opts.do_subplot   =  any(named('--subplot')); 
opts.dynamic_YLIM =  any(named('--set-y'));

opts.label_color = [.4 .4 .4];
if any(named('--label-c')), opts.label_color = get_('--label-c'); end

if any(named('--opts')), function_opts = get_('--opts'); 
else function_opts = []; 
end

list = []; 
varargout = cell(1,nargout-1);
if opts.dynamic_YLIM, varargout = {[]}; end

%% Implement syntax to control which spikes get let through

if any(named('-roi')), include = get_('-roi'); 
    if numel(include) == 2
        include = (data.time >= min(include) & ...
                   data.time <= max(include));
    end
else include = true(size(data.time)); 
end

if any(named('-unit')), u_roi = get_('-unit'); 
  if isnumeric(u_roi), u_roi = ismember(data.unit, u_roi); end
  if any(u_roi), include = include & u_roi; end
end

if any(named('-pass')), pass_ok = get_('-pass');
elseif ~isfield(data,'pass')
     data.pass = ones(size(data.time));
     pass_ok = true;
elseif isfield(data,'data'), % for wave data
     pass_ok = true(size(data.wave,3),1);
else pass_ok = true(size(data.time)); % for spike data
end
if ~islogical(pass_ok), pass_ok = ismember(1:max(data.pass), pass_ok); end
if any(pass_ok), include = include & pass_ok; end

if numel(channel_map) > 1 && opts.do_subplot, clf, end

for pp = 1:numel(channel_map)

    if channel_map(pp) == 0, continue, end
    
    cc = channel_map(pp);
    this_channel = (include & data.channel == cc & pass_ok(data.pass)); 
    if ~opts.hash, this_channel = this_channel & data.unit > 0; end

    if ~any(this_channel) && opts.skip_empty, continue, end
    if opts.do_subplot
      if numel(channel_map) > 1, subplot(subplot_nxy{:}, pp), 
      else cla reset
      end, hold on
    end

    n_units = max(data.unit(this_channel));
    if opts.merge, n_units = 0; end

    for u_id = 0:n_units

      if opts.merge, ok = this_channel;
      else ok = this_channel & data.unit == u_id;
      end
        
      if ~any(ok) && opts.skip_empty, continue, end
      this = data; 

      % Select spikes to include
      for f = fieldnames(data)'
        if numel(ok) ~= size(this.(f{1})), continue, end
        this.(f{1}) = this.(f{1})(ok,:,:,:,:,:,:,:);
      end

      out = cell(size(varargout));

      % Invoke the function to be applied
      if isempty(function_opts)
           [out{:}] = fun(this, [cc u_id], varargin{:});
      else [out{:}] = fun(this, [cc u_id], function_opts);
      end

      list = [list; cc u_id]; %#ok<AGROW> 
      for nn = 1:numel(varargout)
        if isempty(out{nn}), list(end,:) = []; break, end
        if isstruct(out{nn}) && ~isstruct(varargout{nn})
               varargout{nn} = out{nn}; 
         else varargout{nn} = [varargout{nn}; out{nn}]; 
        end
      end
    end

    %% Optional subplot formatting
    if ~opts.do_subplot, continue, end
    
    if opts.dynamic_YLIM
        yl = max(varargout{1}(list(:,1) == cc,1));
        ylim([0 yl])
    end
    plots.tidy
    p = get(gca,'Position'); 
    set(gca,'Position',p + [-1 -1 2 2].*p([3 4 3 4])/10)
    set(gca,'UserData',cc)
    if ~opts.ticks, set(gca,'XTick',[],'YTick',[]); end

    % Set labels
    if opts.label_none, continue, end
    [a,b] = find(channel_map == cc);
    do_label = opts.label_all || ...
         (a == 1 || a == size(channel_map,1)) && ...
         (b == 1 || b == size(channel_map,2));

    if do_label
        text(mean(xlim),mean(ylim),num2str(cc), ...
                'Color',opts.label_color, ...
                'FontSize',14,'horiz','center')
    end
end

if opts.do_subplot
  %% Final subplot formatting  
  h = get(gcf,'Children');
  if any(named('--link')), linkaxes(h,get_('--link')); end
  pos = cat(1,h.Position);
    
  [~,blc] = min(pos*[1;1;0;0]);
  set(h(blc),'XTickMode','auto');
  set(h,'LineWidth',0.8)
end

if nargout == 0, clear, end
return

