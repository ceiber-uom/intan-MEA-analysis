
function ISI (data, varargin)
% function plots.ISI( data, ... )
%
% Basic inter-spike-inteval histogram plot of recorded spike times. 
% 
% Options:
% -chan [c1 c2 c3 ...] Set channels 
% -map  [c1 c2;  ... ] Set channels and arrangement of channels
% -rc   [rows cols]    Set subplot size (# rows / columns)
% -roi  [t0 t1]        Set beginning and end of plotting window
% -unit [u1 u2 ... ]   Filter for the specified units 
% -pass [pass_ids]     Select passes to plot (for expoched data)
% -merge               Compute ISI for all spikes on the channel
% -no-hash             Only show spikes with unit code > 0. 
% -time [bin-times]    Bin times for ISI histogram (default 0:0.5 s, n=64)
%                       if len(times) <= 3 they're expanded using linspace.
% -labels              Add channel labels (default: corners only)
% -no-labels           Suppress channel labels for corners
% -log                 Plot ISI on log scale (default if -raster set)
% -raster              Plot ISI-spiketime raster instead of histogram
% 
% v0.1 - 23 December 2022 - Calvin Eiber <c.eiber@ieee.org>

named = @(s) strncmpi(s,varargin,numel(s));
get_ = @(v) varargin{find(named(v))+1};

opts.log   = any(named('-log')); 
opts.raster = any(named('-ras'));
opts.bins = linspace(0, 0.5, 64);

if opts.log, opts.bins = linspace(-3, 1, 64); end
if any(named('-tim')), opts.bins = get_('-tim');
  if numel(opts.bins) == 1, opts.bins = linspace(0, opts.bins, 64);
  elseif numel(opts.bins) == 2, 
    opts.bins = linspace(0, opts.bins(1), opts.bins(2));
  elseif numel(opts.bins) == 3, 
    opts.bins = linspace(opts.bins(1), opts.bins(2), opts.bins(3));
  end
end

opts.merge = any(named('-merge')); 

tools.forChannels(data, @isi_plot, varargin{:}, ...
                   '--opts', opts, '--subplot', '--set-y')

h = get(gcf,'Children');
pos = cat(1,h.Position);
[~,blc] = min(pos*[1;1;0;0]);

if opts.raster, xlabel(h(blc),'time (s)'), ylabel(h(blc),'ISI (s)')
elseif opts.log, xlabel(h(blc),'log_{10} ISI (s)')
else xlabel(h(blc),'ISI (s)')
end


function yl = isi_plot(data, index, opts)
        
isi = diff(data.time); 
bins = opts.bins;

if index(2) > 0, color = lines(index(2)); color = color(end,:);
else color = [.6 .6 .6];
end

if ~opts.raster
  if opts.log, [y,x] = hist(log10(isi), bins); %#ok<HIST> 
  else         [y,x] = hist(isi, bins);        %#ok<HIST> 
  end
end

style = {'FaceColor',color,'FaceAlpha',0.5,'EdgeColor','none'};
if opts.merge, style{4} = 1; end

if opts.raster
  isi(end+1) = median(isi);
  plot(data.time,isi,'.','Color',color)
  if ~opts.log, set(gca,'YScale','log'), end
  yl = max(ylim);
else 
  bar(x,y,1,style{:})
  yl = max(y(2:end-1));
end

return

