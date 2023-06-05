

n_cells = numel(epochs.SPIKE.clusters); 

do_PDF = true;
plots.PDF_tools('setup',do_PDF);

for cell_id = 1:n_cells

    clf    
    p = get(gca,'position'); % default single axes

    % this is the code which actually generates the plot
    plots.response_curve(epochs,'-per-unit','-chan',cell_id, '-x', ...
                        'average','-roi',[-0.05 0.2]); %  see just this unit
    
    set(gca,'position',p), axis tight, xlim([-0.1 10])
    
    xlabel('stimulus intensity')
    
    
    h = get(gca,'children');
    
    for dur = 1:7
    
        ax = copyobj(h(1).Parent,gcf);
        ax.Position(4) = ax.Position(4) / 8;
        ax.Position(2) = ax.Position(2) + (dur-1)*ax.Position(4)*(8/6.5);
    
        cull = setdiff(1:10,[dur 9]);
    
        delete(ax.Children(cull))
    
        if dur > 1
            set(ax,'XTick',[],'YTickLAbel','')
            ylabel(ax,''), xlabel(ax,'')
        end
    
        ax.Children(2).Color = [0 0 0 0.3];
        ax.Children(2).XData(1) = 0;
    
        info = ax.Children(1).UserData;
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

