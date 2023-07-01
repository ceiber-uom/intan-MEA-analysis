
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
% -single-figure, -1f  Old default option, create one figure with subplots
%                      instead of one figure per cell (new default)
% -ticks               Toggle Ticks display for subplots
% -no-tidy             Supress plots.tidy and basic subplot format
% -labels              Label each subplot (default: corners only)
% -no-labels           Suppress labeling of subplots, even corners. 
% --label-color [rgb]  Set label color, default [.4 .4 .4] (grey)
% 
% v0.1 - 23 December 2022 - Calvin Eiber <c.eiber@ieee.org>

named = @(s) strncmpi(s,varargin,numel(s)); % varargin parser
get_ = @(v) varargin{find(named(v))+1};

if isfield(data,'config') % is this a top-level INTAN object? if so call forWaveType. 
   varargout = cell(1,nargout-1);
   this = @(d) tools.forChannels(d, fun, varargin{:}); % this function
   if nargout == 0, tools.forWaveType(data, this, varargin{:},'--get-types');
   else [list,varargout{:}] = tools.forWaveType(data, this, varargin{:},'--get-types');
   end
   if nargout > 0, list = list.(types{1}); end %#ok<USENS> 
   return
elseif numel(data) > 1 || any(named('--undo')) % e.g. the output from tools.simplify()
  % The --undo flag is added to force this behaviour in the case of plotting just one cell  
    data = tools.simplify(data,'-undo');
end

channel_map = plots.layout(data, varargin{:});
subplot_nxy = fliplr(num2cell(size(channel_map)));

% control over how the data is split up and what gets passed through
opts.hash  = ~any(named('-no-hash'));
opts.merge = any(named('-ignore-u')) || any(named('-merge'));
opts.skip_empty   = ~any(named('-enable-empty'));

% Control over single/multi figure and figure labels
opts.do_subplot   = any(named('--subplot')); 
opts.ticks        =(numel(channel_map) > 6) == any(named('-tick'));
opts.one_figure   = any(named('-single-f')) || any(named('-1f'));
opts.label_all    = any(named('-label'));
opts.label_none   = any(named('-no-l'));
opts.dynamic_YLIM = any(named('--set-y'));

if opts.one_figure, 
    opts.ticks = true; 
    opts.label_none = ~opts.label_all; 
end

opts.label_color = [.4 .4 .4];
if any(named('--label-c')), opts.label_color = get_('--label-c'); end

if any(named('--opts')), function_opts = get_('--opts'); 
else function_opts = []; 
end

list = []; 
varargout = cell(1,nargout-1); % because arg_out[1] is list
if opts.dynamic_YLIM, varargout = {[]}; end 

%% Implement syntax to control which spikes get let through

% select passes or all passes 
if any(named('-pass')), pass_ok = get_('-pass');
elseif ~isfield(data,'pass')
     data.pass = ones(size(data.time));
     pass_ok = 1;
elseif isfield(data,'data'), % for wave data
     pass_ok = true(size(data.wave,3),1);
else pass_ok = 1:max(data.pass); % for spike data
end
if ~islogical(pass_ok), pass_ok = ismember(data.pass, pass_ok); end

if any(pass_ok) && ~isfield(data,'data'), include = pass_ok;
else include = true(size(data.time));
end

% set time window ROI
if any(named('-roi')), roi_window = get_('-roi'); 
  if numel(roi_window) == 1, roi_window = [0 roi_window];
  elseif numel(roi_window) ~= 2, 
      error('-roi must be a one- or two-element vector')
  end
  if isfield(data,'pass_begin') % spike times need to be adjusted
    t0 = [0; data.pass_begin];
    t_spk = data.time - t0(data.pass+1);

    include = include & (t_spk >= min(roi_window) & ...
                         t_spk <= max(roi_window));

  else
    include = include & (data.time >= min(roi_window) & ...
                         data.time <= max(roi_window));
  end
end

% set unit ROI (this makes more sense with a ch.unit where unit is e.g.
% 0,1,2), KiloSort does it differently with more like a cell_id across 
% channels.  

if any(named('-unit')), u_roi = get_('-unit'); 
  if isnumeric(u_roi), u_roi = ismember(data.unit, u_roi); end
  if any(u_roi), include = include & u_roi; end
end

if numel(channel_map) > 1 && opts.do_subplot, clf, end

%% Core loop

for pp = 1:numel(channel_map)

    if channel_map(pp) == 0, continue, end
    
    cc = channel_map(pp);
    this_channel = (include & data.channel == cc); 
    if ~opts.hash, this_channel = this_channel & data.unit > 0; end

    if ~any(this_channel) && opts.skip_empty, continue, end
    if opts.do_subplot
      if ~opts.one_figure, % will make figure inside of unit loop
      elseif numel(channel_map) > 1, subplot(subplot_nxy{:}, pp), 
      else cla reset
      end, hold on
    end

    n_units = max(data.unit(this_channel));
    if opts.merge, n_units = 0; end

    for u_id = 0:n_units

      if opts.merge, ok = this_channel;
      else ok = this_channel & data.unit == u_id;
      end

      if ~opts.one_figure,
        figure('Name',sprintf('channel %d, unit %d', cc, u_id))
      end
        
      if ~any(ok) && opts.skip_empty, continue, end
      this = data; 

      % Select spikes to include
      for f = fieldnames(data)'  %' filter variables to this spike
        if numel(ok) ~= size(this.(f{1})), continue, end
        this.(f{1}) = this.(f{1})(ok,:,:,:,:,:,:,:);
      end

      % Select per-pass variables to include
      if isfield(this,'index') && size(this.index,2) >= 2

          idx = all(this.index == [cc u_id],2);
    
          % Select pass variables to include
          for f = fieldnames(data)' %' filter to this pass
            if numel(idx) ~= size(this.(f{1})), continue, end
            this.(f{1}) = this.(f{1})(idx,:,:,:,:,:,:,:);
          end
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

if opts.do_subplot && opts.one_figure
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

