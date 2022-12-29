
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
% -t [bin-times]       Bin times for ISI histogram (default 0:0.5 s, n=64)
% -labels              Add channel labels (default: corners only)
% -no-labels           Suppress channel labels for corners
% 
% v0.1 - 23 December 2022 - Calvin Eiber <c.eiber@ieee.org>

named = @(s) strncmpi(s,varargin,numel(s));
get_ = @(v) varargin{find(named(v))+1};

opts.log   = any(named('-log')); 
opts.bins = linspace(0, 0.5, 64);

if opts.log, opts.bins = linspace(-3, 1, 64); end
if any(named('-t')), opts.bins = get_('-t');
  if numel(opts.bins) == 1, opts.bins = linspace(0, opts.bins, 64);
  elseif numel(opts.bins) == 2, 
    opts.bins = linspace(0, opts.bins(1), opts.bins(2));
  end
end

opts.merge = any(named('-merge')); 

tools.forChannels(data, @isi_plot, varargin{:}, ...
                   '--opts', opts, '--subplot', '--set-y')

h = get(gcf,'Children');
pos = cat(1,h.Position);
[~,blc] = min(pos*[1;1;0;0]);

if opts.log, xlabel(h(blc),'log_{10} ISI (s)')
else xlabel(h(blc),'ISI (s)')
end



function yl = isi_plot(data, index, opts)
        
isi = diff(data.time); 
bins = opts.bins;

if index(2) > 0, color = lines(index(2)); color = color(end,:);
else color = [.6 .6 .6];
end

if opts.log, [y,x] = hist(log10(isi), bins); %#ok<HIST> 
else         [y,x] = hist(isi, bins);        %#ok<HIST> 
end

style = {'FaceColor',color,'FaceAlpha',0.5,'EdgeColor','none'};
if opts.merge, style{4} = 1; end

bar(x,y,1,style{:})
yl = max(y(2:end-1));

return

