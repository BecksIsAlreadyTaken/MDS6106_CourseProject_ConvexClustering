function [f] = generate_plot(result,label,path,xlimit,ylimit,plot_title,showfig)
    f = figure;  
    if showfig == false
        set(f,'visible','off');
    end
    hold on;
    for j = 1:size(label,1)
        mask = result(:,1) == label(j);
        col = j.*ones(1,size(result(mask,:),1));
        scatter(result(mask,2),result(mask,3),[],col,'filled');
    end
    hold off;
    title(plot_title);
    xlim(xlimit);
    ylim(ylimit);
    saveas(gcf,path);
end