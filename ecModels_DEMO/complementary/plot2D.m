function plot2D(dataX,dataY,titleStr,xStr,yStr,cumDist)
    figure   
    if ~cumDist
        plot(dataX,dataY,'LineWidth',5)
    else
        [f, x] = ecdf(dataX);
        plot(x,f,'LineWidth',5)
        set(gca, 'XScale', 'log')
        xlim([0 1000])
        ylim([0 1])
        %axis(axisLimits)
    end
    title(titleStr)
    xlabel(xStr)
    ylabel(yStr)
end