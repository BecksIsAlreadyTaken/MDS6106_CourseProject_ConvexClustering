function [f] = convergence_plot(ng,path,plot_title,showfig)
    f = figure;
    if showfig == false
        set(f,'visible','off');
    end
    hold on;
    plot(1:size(ng,1),log10(ng),'LineWidth',1);
    hold off;
    ylim([-1 inf]);
    xlim([1 inf]);
    xlabel('Number of iterations');
    title(plot_title);
    saveas(gcf,path);
end