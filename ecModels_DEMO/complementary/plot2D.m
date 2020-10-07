function plot2D(dataX,dataY,titleStr,xStr,yStr,cumDist)
    %figure   
    if ~cumDist
        plot(dataX,dataY,'LineWidth',5)
        set(gca,'FontSize',14)
    else
        [f, x] = ecdf(dataX);
        plot(x,f,'LineWidth',5)
        set(gca, 'XScale', 'log','FontSize',14)
        xlim([0 1000])
        ylim([0 1])
        %axis(axisLimits)
    end
    title(titleStr,'FontSize',18)
    xlabel(xStr,'FontSize',18)
    ylabel(yStr,'FontSize',18)
end