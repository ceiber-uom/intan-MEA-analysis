
% This assumes that data has already been loaded and is this file:
assert(contains(epochs.filename,'01_SQ1100NF565040423'))

n_cells = numel(epochs.SPIKE.clusters); 

do_PDF = true;
plots.PDF_tools('setup',do_PDF);

for cell_id = 1:n_cells % for each cell_id

    clf    
    p = get(gca,'position'); % default single axes

    % this is the code which actually generates the plot
    plots.response_curve(epochs,'-per-unit','-chan',cell_id, '-x', ...
                        'average','-roi',[-0.05 0.2]); %  see just this unit
    
    % modify the axis position and x/y limits
    set(gca,'position',p), axis tight, xlim([-0.1 10]), ylim(ylim);
    xlabel('stimulus intensity') % was 'average'
    
    h = get(gca,'children');
    nD = length(findobj(h,'Type','ErrorBar')); % number of durations

    for dur = 1:nD % make into small-replicates axes by handle manipulation
    
        ax = copyobj(h(1).Parent,gcf); % original axes with all plots/data
        ax.Position(2) = ax.Position(2) + (dur-1)*ax.Position(4)/(nD-0.5);
        ax.Position(4) = ax.Position(4) / (nD+1);
        
        cull = setdiff(1:numel(h),[dur numel(h)-1]); % this one and baseline
        delete(ax.Children(cull))
    
        if dur > 1 % remove axis ticks and labels
            set(ax,'XTick',[],'YTickLAbel','')
            ylabel(ax,''), xlabel(ax,'')
        end
    
        ax.Children(2).Color = [0 0 0 0.3]; % make feint grey
        ax.Children(2).XData(1) = 0; % make go to zero
    
        info = ax.Children(1).UserData; % get the stimulus info
        % this is put there by the 
        text(max(xlim),ylim*[.1;.9], ...
             sprintf('duration = %d ms', info.duration*1e3), ...
             'FontSize',10,'Color',ax.Children(1).Color, ...
             'Horiz','right','parent',ax)
    end
    
    delete(h(1).Parent)
    
    suptitle(sprintf('[%d/%d] cell %d : channel %d, quality %d', ...
                     cell_id, n_cells, ...
                     data.SPIKE.sort_info(cell_id).clus, ...
                     data.SPIKE.sort_info(cell_id).ch, ...
                     data.SPIKE.sort_info(cell_id).quality))
    
     pause(0.5)

     plots.PDF_tools(gcf, do_PDF, 'page-%04d.ps', cell_id)

end

plots.PDF_tools('compile',do_PDF,'my-analysed-data (%d).pdf')

