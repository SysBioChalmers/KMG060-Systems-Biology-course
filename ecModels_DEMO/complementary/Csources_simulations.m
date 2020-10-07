function [flux_dist, MRE_tot,conditions] = Csources_simulations(ecModel_batch)
% Csources_simulations
% 
% Ivan Domenzain. created:       2017-09-21
% Ivan Domenzain. Last modified: 2020-10-06

data    = readtable('../data/growthRates_data_carbonSources.txt','delimiter','\t','ReadVariableNames',true);
media   = {'Min' 'MAA' 'YEP'};
conditions = [];
figure
axis square
%Create a table for storing all flux distributions
flux_dist = table(ecModel_batch.rxns,ecModel_batch.rxnNames,'VariableNames',{'rxns' 'rxnNames'});
%Identify position of the growth pseudoreaction in model structure
gR_pos = find(strcmpi(ecModel_batch.rxnNames,'growth'));
%Set growth as objective to maximize
ecModel_batch.c(:) = 0;
ecModel_batch.c(gR_pos) = 1;
MRE_tot = [];
markers = {'o' 's' 'd'};
colors  = {'black' 'red' 'blue'};
legendStr = {};
for i=1:length(media)
    %Identify all entries for the i-th media type
    mediaIdxs = find(strcmpi(data.media,media{i}));
    relErr = [];
    gRates_sim = zeros(length(mediaIdxs),1);
    for j=1:length(mediaIdxs)
        model = ecModel_batch;
        %create string for the corresponding carbon source uptake reaction
        %name, then set media conditions with unbounded cSource uptake
        c_source = [data.c_source{mediaIdxs(j)} ' exchange (reversible)'];
        model    = changeMedia_batch(model,c_source,data.media{mediaIdxs(j)});
        %Get a flux distribution
        solution      = solveLP(model);
        gRates_sim(j) = solution.x(gR_pos);
        %Get relative error in gRate prediction
        res    = round(abs(data.gRate(mediaIdxs(j))-gRates_sim(j))*100/data.gRate(mediaIdxs(j)),2);
        relErr = [relErr;res];
        str    = [data.media{mediaIdxs(j)} '_' data.short_name{mediaIdxs(j)}];
        conditions = [conditions;str];
        %Append flux distribution to results table
        eval(['flux_dist.' str ' = solution.x;'])
        %Place c source str in the corresponding coordinates
        
    end
    %hold on
    plot(data.gRate(mediaIdxs),gRates_sim,markers{i},'MarkerSize',15,'MarkerFaceColor',colors{i},'MarkerEdgeColor',colors{i})
    text(1.06*data.gRate(mediaIdxs),1.06*gRates_sim,data.short_name(mediaIdxs),'FontSize',14)
    hold on
    MRE_tot = [MRE_tot;relErr];
    MRE = mean(relErr);
    legendStr{i} = [media{i} ' / MRE =' num2str(MRE) '%'];
end
%Calculate overall mean relative error
MRE_tot = mean(MRE_tot);
set(gca,'FontSize',14)
%title('Max growth rate on different carbon sources','FontSize',30,'FontWeight','bold')
ylabel('\mu_{max} predicted [h^{-1}]','FontSize',18,'FontWeight','bold');
xlabel('\mu_{max} experimental [h^{-1}]','FontSize',18,'FontWeight','bold');
xlim([0 0.6])
ylim([0 0.6])
x1 = linspace(0,1,1000);
%legend(legendStr,'FontSize',16,'Location','southeast')
legend(legendStr,'FontSize',16,'Location','southeast')
plot(x1,x1)
hold off

end
