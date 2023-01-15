
function culm (data, varargin)
% function plots.culm( data, [Y], [cell_ids], ... )
%
% Cumulative density of recorded spike counts (or other channel data)
% 
% Options:
% -y [y1 y2 ... ]      Set data to visualuse 
% -chan [c1 c2 c3 ...] Set channels
% -rc   [rows cols]    Set subplot size (# rows / columns)
% -roi  [t0 t1]        Set beginning and end of plotting window
% -pass [pass_ids]     Select passes to plot (for expoched data)
% -unit [u1 u2 ... ]   Filter the displayed units 
% -scale               Toggle colorscale (default: enabled if subplot 111)
% -count               Plot spike count (default: mean spikerate)
%                      (only if Y data is not otherwise specified)
% 
% v0.1 - 25 December 2022 - Calvin Eiber <c.eiber@ieee.org>

named = @(s) strncmpi(s,varargin,numel(s));
get_ = @(v) varargin{find(named(v))+1};

if isfield(data,'config')
   this = @(d) plots.culm(d, varargin{:});
   tools.forWaveType(data, this, varargin{:});
   return
end

y_label = '';

chan_unit = [];
if any(named('-cel')), chan_unit = get_('-cel'); end
if any(named('-y')), y_data = get_('-y');
elseif nargin > 1 && isnumeric(varargin{1}), y_data = varargin{1}; 
    if nargin > 2 && isnumeric(varargin{2}), chan_unit = varargin{2}; end
else
    %% Select spikes to include by applying filters
    [chan_unit, y_data] = tools.forChannels(data, ...
                                       @(d,varargin) length(d.time), ...
                                               '--ordered', varargin{:}); 
    y_label = 'spike count (imp)';
end

if isempty(chan_unit), 
    [chan_unit,~] = tools.forChannels(data, @(d,varargin) 1, '--ordered');
end

if ~any(named('--no-t'))
  if size(y_data,2) > size(y_data,1), y_data = y_data'; end
  if size(chan_unit,2) > size(chan_unit,1), chan_unit = chan_unit'; end
end

if size(y_data,1) ~= size(chan_unit,1)
    error('the size of Y (%d entries) does not match %s (%d entries)', ...
           size(y_data,1),'the number of cell IDs',size(chan_unit,1))
end

%%

n_init = numel(get(gca,'children')); 
hold_mode = get(gca,'NextPlot');
if strcmp(hold_mode,'replace'), cla reset, n_init = 0; end

for ii = 1:size(y_data,2)
    [y,seq] = sort(y_data(:,ii));
    h = plot( y, linspace(0,100,numel(y)),'linewidth',1.2);
    hold on

    if any(named('-no-i')), continue, end
    set(h,'ButtonDownFcn',@(a,b)tooltip(a,b,y,chan_unit(seq,:)))
end

if any(named('-tidy')) == (n_init == 0), return, end
plots.tidy, xlabel(y_label)
set(gca,'YTickMode','manual', ...
        'YTickLabel',strcat(get(gca,'YTickLabel'),'%'));
set(gca,'NextPlot',hold_mode)


return


function tooltip(hobj,eventdata,y,ids)
%%
delete(findobj(hobj.Parent,'UserData','tooltip')) % remove previous
ii = round(numel(y)*eventdata.IntersectionPoint(2)/100); % add new

x = 100*ii/numel(y);

hold_mode = get(gca,'NextPlot');

hold on
plot(y(ii),x,'o','MarkerFaceColor',hobj.Color,'Color',hobj.Color,'UserData','tooltip','HitTest','off')

label = sprintf('.%g', ids(ii,:));

text(y(ii),x-1,sprintf('%0.1f%%: %s', x, label(2:end)), ...
    'FontSize',8,'Color',hobj.Color,'VerticalAlignment','top','userData','tooltip','HitTest','off')
set(gca,'NextPlot',hold_mode)


return


