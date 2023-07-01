
function responseCurve (data, varargin)
% function plots.responseCurve( data, ... )
%
% Response curve analysis of recorded spikes. For simple stimuli this code 
% works almost automatically to figure out what way of plotting makes the
% most sense, for more complicated stimuli (i.e. parameter grid protocols)
% you might need to give it a bit of help. 
% 
% by default this will plot data.SPIKES.spike_count (in the post-stimulus
% interval i.e. t>0) and data.SPIKES.base_count (in the pre-stimulus
% interval i.e. t<0). If these fields don't exist yet they are computed 
% 
% Options:
% -chan [c1 c2 c3 ...] Set channels 
% -map  [c1 c2;  ... ] Set channels and arrangement of channels
% -unit [u1 u2 ... ]   Filter the displayed units 
% -rc   [rows cols]    Set subplot size (# rows / columns)
% -x [variable]        Set X axis variable (default: from epochs)
% -y [variable]        Set Y axis variable (default: spike count)
% 
% -merge               Compute response curve for all spikes on channel
% -no-hash             Only show spikes with unit code > 0. 
% -roi [start end]     Set analysis window ROI (in s)
% -labels              Add channel labels (default: corners only)
% -no-labels           Suppress channel labels for corners
% -ticks               Toggle default ticks behaviour (ticks enabled for
%                       6 or fewer channels, disabled otherwise)
% -per-unit            One panel per unit (default: one per channel)
% -no-average          Do not attempt averaging over conditions
% -pass [pass_ids]     Select passes to plot
% 
% Dependencies: tools.forChannels
% 
% v0.1 - 4 June 2023 - Calvin Eiber <c.eiber@ieee.org>

named = @(s) strncmpi(s,varargin,numel(s));
get_ = @(v) varargin{find(named(v))+1};

if isfield(data,'config')
   this = @(d) plots.responseCurve(d, varargin{:},'--epochs',data.epochs);
   tools.forWaveType(data, this, varargin{:});
   return
end

if any(named('-per-u'))
    [~,~,data.channel] = unique([data.channel data.unit],'rows');
end

if any(named('--epochs')), epochs = get_('--epochs');
    epochs.fs = mean((epochs.finish-epochs.start)./epochs.duration);
else error guestimate_epoch_structure
end

opts.epochs = epochs;
opts.log_x = false;

if any(named('-log-x')), opts.log_x = true; end

if any(named('-x'))
     opts.x_var = get_('-x');
     opts.x_label = strrep(get_('-x'),'_',' ');
else % attempt to auto-detect
    opts = guess_axis_structure(epochs, opts);
end

if any(named('-y')), opts.y_var = get_('-y'); 
     opts.y_label = strrep(get_('-y'),'_',' ');
     opts.y_style = {'o','s','v','x'};
else
    %% Select spikes to include by applying filters
    opts.y_label = 'spike count (imp)';
    opts.y_var = {'spike_count','base_count'};
    opts.y_style = {'o','-(0.3)'};
    
    if ~isfield(data,'spike_count')
      [chan_unit, y_data, y_base] = tools.forChannels(data, @pass_spike_counts, ...
                                                 '--ordered', varargin{:}); 
      data.spike_count = y_data;
      data.base_count = mean(y_base,2);
      data.index = chan_unit;
    end
end

if any(named('-sty')), opts.y_style = get_('-sty'); end

tools.forChannels(data, @plot_curve, varargin{:}, ...
                   '--opts', opts, '--subplot')

h = get(gcf,'Children');
pos = cat(1,h.Position);
[~,blc] = min(pos*[1;1;0;0]);

xlabel(h(blc),opts.x_label); ylabel(h(blc),opts.y_label);

set(h,'YLim',[min([h.YLim]) max([h.YLim])])

arrayfun(@(x) clean_up_plot(x,opts), h) % fix axis labelling

if isfield(epochs,'condition_id') && ~any(named('-no-av')) && ... 
        numel(unique(epochs.condition_id)) < numel(epochs.condition_id)

    arrayfun(@(x) average_across_replicates(x,opts), h) % fix axis labelling

end

return

function plot_curve(data, index, opts)


if ~iscell(opts.y_var), opts.y_var = {opts.y_var}; end

if index(2) > 0, color = lines(index(2)); color = color(end,:);
else color = [.6 .6 .6];
end


for ii = 1:numel(opts.y_var)
    
    x = opts.epochs.(opts.x_var);
    y = data.(opts.y_var{ii});
       
    mkr = opts.y_style{mod(ii-1,numel(opts.y_style))+1};
    style = {'Color',color,'lineWidth',1.1,'markerSize',5};  

    if any(mkr == '(')

        alpha = str2double(regexp(mkr,'0\.\d+','match'));
        if ~isempty(alpha), style{2}  = [color alpha(1)]; end
        % in theory could put other stuff here as well. 

        mkr(find(mkr == '(',1):end) = []; 
    end

    if numel(y) == 1 && any(mkr == '-'), 
      x = [min(x) max(x)];
      y = [y y];
    end


    plot(x,y,mkr,style{:},'userdata',[index ii])
    hold on
end

if opts.log_x, set(gca,'XScale','log'); end

return

function clean_up_plot(h, opts)

t = findobj(h,'type','text');
if isempty(t), return, end
if opts.log_x, t(1).Position(1) = sqrt(prod(h.XLim)); end
t(1).Position(2) = max(h.YLim);
t(1).VerticalAlignment = 'top';

return

function [n_spike,n_spk_pre] = pass_spike_counts(data, ~, varargin)

nP = numel(data.pass_begin);

n_spike = zeros(1,nP);
n_spk_pre = zeros(1,nP);

for pp = 1:nP
    t = data.time(data.pass == pp) - data.pass_begin(pp);    
    n_spike(pp) = sum( t>0 );
    n_spk_pre(pp) = sum(t <= 0);
end

return

function opts = guess_axis_structure(epochs, opts)
%% Determine best-spaced X-axis variable
consistency = [];
selection   = [];

vars = {'duration','frequency','amplitude','phase','average'};

if isfield(epochs,'block_id'), 
    error TODO_implement_block_logic
end

for ii = 1:numel(vars)

    x = sort(epochs.(vars{ii}));
    if any(x > 1e6), continue, end % bad value = do not use. 

    dx = diff(x); 
    dx(dx < min(x)/100) = []; % remove replicates

    consistency(end+1,:) = (mean(dx) / std(dx)) * sqrt(numel(dx));
    selection(end+1,:) = [ii 0];
        
    if any(x < 0), continue, end % do not consider log-transform

    x = log10(x);
    dx = diff(x); 
    dx(dx < 0.02) = []; 

    consistency(end+1,:) = mean(dx) / std(dx) * sqrt(numel(dx));
    selection(end+1,:) = [ii 1];

end

[~,sel] = max(consistency); 

opts.x_var = vars{selection(sel,1)};
opts.x_label = sprintf('%s (%s)', strrep(opts.x_var,'_',' '), ... 
                                   epochs.units.(opts.x_var));
opts.log_x = (selection(sel,1) ~= 0);


function average_across_replicates(h, opts)

epochs = opts.epochs;
cond = unique(epochs.condition_id);

x = h.Children(end).XData;
y = h.Children(end).YData;

x_avg = zeros(size(cond));
y_avg = zeros(size(cond));
y_sem = zeros(size(cond));

for cc = 1:numel(cond)

    sel = (epochs.condition_id == cond(cc));
    y_avg(cc) = mean(y(sel));
    y_sem(cc) = std(y(sel)) / sqrt(sum(sel));
    x_avg(cc) = mean(x(sel));

end

% guess covariate, fewest 
covariate = []; 
covariate_label = {}; 
vars = {'duration','frequency','amplitude','phase','average'};
cov = @(x) std(x) / mean(x);

for ii = 1:numel(vars)
  if strcmp(opts.x_var,vars{ii}), continue, end

  civ = arrayfun(@(c) cov(epochs.(vars{ii})(epochs.condition_id == c)), ...
                       cond);

  % disp(mean(civ))
  if mean(civ) > 0.01, continue, end % too much variance within condition ID

  if numel(unique(epochs.(vars{ii}))) > 2*numel(cond), continue, end

  covariate = [covariate epochs.(vars{ii})]; 
  covariate_label = [covariate_label vars(ii)];
end

[covariate_levels,~,cov_id] = unique(covariate,'rows');

%% Make errorbars

if numel(covariate_levels) > 3
     color = turbo(size(covariate_levels,1)) .* [1 .9 1];
else color = lines(7); 
end

for ii = 1:numel(covariate_levels)

    sel = ismember(cond, epochs.condition_id(cov_id == ii));
    if mod(ii,2) == 1, mfc = 'w';
    else mfc = (color(ii,:) + 0.3)/1.3;
    end

    info = struct;
    for lbl = 1:numel(covariate_label)
        info.(covariate_label{lbl}) = covariate_levels(ii,lbl); 
    end

    errorbar(x_avg(sel),y_avg(sel),y_sem(sel),'o', ...
        'LineWidth',1.1,'Color',color(ii,:),'MArkerFaceColor',mfc, ...
        'PArent',h,'userdata',info)

end

h.Children(end).Visible = 'off';

return