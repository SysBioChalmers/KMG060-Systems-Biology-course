%ecModels DEMO
%
%  - Enzyme-constrained model exploration
%  - The role of enzyme constraints in growth simulations
%  - Growth on different carbon sources
%  - Crabtree effect simulation
%  
%  based on the latest
%  version of ecYeastGEM (GitHub: https://github.com/SysBioChalmers/ecModels/ecYeastGEM/model)
%
% Ivan Domenzain.	Last modified 2020-10-06

%% Enzyme-constrained model exploration
%Load ecModel_batch with total protein pool constraint
load('models/ecYeastGEM_batch.mat');
%Load original model (yeastGEM)
load('models/yeastGEM.mat');
%Get the total number of promiscuous enzymes isoenzymes and rxns with Kcats
cd complementary
[isoEnzymes,Promiscuous,Complexes,RxnWithKcat] =  rxnCounter(ecModel_batch,model);

%% Visualizing solution vectors as cumulative distributions
%set minimal mineral medium and allow a unit uptake of glucose [mmol/gDwh]
tempModel  = minimal_Y6(model,glucIN,-1); 
tempEcModel = changeMedia_batch(ecModel_batch,'D-glucose exchange (reversible)','Min',1); 
%take a look to the minimal_Y6 and changeMedia_batch functions. They're used
%allow the uptake of metabolites that are present in a given media
%formulation and close the uptake for those that are not present

%solve linear programing problem using the "solveLP" function from
%RAVEN
solution   = solveLP(model);
solutionEC = solveLP(tempEcModel);
%This will show a cumulative distribution for the absolute values of the
%fluxes in the simulation
FluxDist = solution.x(abs(solution.x)>1E-7);
plot2D(abs(FluxDist),[],'FBA simulation (yeastGEM)','Fluxes value [mmol/gDCW h]','Cumulative distribution',true)
%Look at the new flux cumulative distribution
FluxDist = solutionEC.x(abs(solutionEC.x)>1E-7);
plot2D(abs(FluxDist),[],'FBA simulation (ecYeastGEM)','Fluxes value [mmol/gDCW h]','Cumulative distribution',true)

Growth rate vs Glucose uptake rate
%Find relevant reaction positions first
growthIndex   = strcmpi(model.rxnNames,'growth');
growthIndexEC = strcmpi(ecModel_batch.rxnNames,'growth');
glucIN        = model.rxns{strcmpi(model.rxnNames,'D-glucose exchange')};
glucIndex     = strcmpi(model.rxnNames,'D-glucose exchange');
c_source      = 'D-glucose exchange (reversible)'; 
glucIndexEC   = find(strcmpi(ecModel_batch.rxnNames,c_source));
%% Plot growth rate vs GUR
gRates    = 0; %Initialize vectors of results
GURates   = 0;
gRatesEC  = 0;
GURatesEC = 0;
%Perform a growth maximization for both GEM and ecModel at each level of
%glucose uptake rate (GUR)
for i=1:18
    %Set new lb for glucose uptake at every iteration
    GUR = i;
    %Set minimal media conditions (just allow uptake of glucose and
    %essential minerals)
    tempModel      = minimal_Y6(model,glucIN,-GUR); 
    [tempEcModel,~] = changeMedia_batch(ecModel_batch,'D-glucose exchange (reversible)','Min',GUR); 
    %solve linear programing problem
    solution       = solveLP(tempModel);
    solutionEC     = solveLP(tempEcModel);
    %If the simulations are feasible then save the results
    if ~isempty(solution.x)
        gRates  = [gRates; solution.x(growthIndex)];
        GURates = [GURates;-solution.x(glucIndex)];
    end
    if ~isempty(solutionEC.x)
        gRatesEC  = [gRatesEC; solutionEC.x(growthIndexEC)];
        GURatesEC = [GURatesEC;solutionEC.x(glucIndexEC)];
    end
end
%Plot results
plot2D(GURates,gRates,'','GUR [mmol/gDw h]','Growth rate [g/gDCW h]',false)
hold on
plot2D(GURatesEC,gRatesEC,'','GUR [mmol/gDw h]','Growth rate [g/gDCW h]',false)
legend({'yeastGEM' 'ecYeastGEM'},'Location','best','FontSize',16)
hold off

%% Growth on diverse environments
%Simulate max. growth on diverse environments using different media composition
%and different carbon sources. No numerical constraints are set by the
%Csources_simulations function
[flux_dist, MRE,conditions] = Csources_simulations(ecModel_batch);
hold off

%% Simulate Crabtree-effect
% Run 50 chemostat simulations in the range 0 growth - max. growth, store
% exchange fluxes and total protein demandsfrom each simulation and plot 
% them against dilution rate

%Get the indices or IDs for some relevant exchange fluxes
Csource_id  = ecModel_batch.rxns(glucIndexEC);
oxyIndex    = find(strcmpi(ecModel_batch.rxnNames,'oxygen exchange (reversible)'));
CO2Index    = find(strcmpi(ecModel_batch.rxnNames,'carbon dioxide exchange'));
ethIndex    = find(strcmpi(ecModel_batch.rxnNames,'ethanol exchange'));
aceIndex    = find(strcmpi(ecModel_batch.rxnNames,'acetate exchange'));
protIndex   = find(strcmpi(ecModel_batch.rxnNames,'prot_pool_exchange'));
exchIndexes = [glucIndexEC;oxyIndex;CO2Index;ethIndex;aceIndex;protIndex];
%Rxn Id for biomass pseudoreaction
bioRxn   = 'r_4041';
%Set minimal media constraints and growth as objective function
ecModel  = changeMedia_batch(ecModel_batch,c_source,'Min');
ecModel  = setParam(ecModel,'obj',{bioRxn},1);
solution = solveLP(ecModel,1);
%Get a max. growth rate to delimit range for chemostat simulations
if ~isempty(solution.f)
    maxGrowth = solution.x(strcmpi(ecModel.rxns,bioRxn));
    disp(['The maximum growth rate for the ecModel_batch on ' c_source ' minimal media is: ' num2str(maxGrowth) ' [g/gDw h]'])
end
%block production of metabolites according to reported experimental data
%DOI: 10.1128/AEM.64.11.4226-4233.1998
ecModel = setParam(ecModel,'ub',ecModel.rxns(aceIndex),0.6); %acetate
ecModel = setParam(ecModel,'ub','r_2033',0.05); %pyruvate
ecModel = setParam(ecModel,'ub','r_1808',0.15); %glycerol
ecModel = setParam(ecModel,'ub','r_1631',0); %acetaldehyde
ecModel = setParam(ecModel,'ub','r_1793',0);%formate

%initialize some variables
D_crit  = [];
results = [];
i=0;
figure
for subopt_growth=0:(maxGrowth/50):maxGrowth
    tempModel = ecModel;
    %Fix suboptimal growth rate
    tempModel = setParam(tempModel,'lb',{bioRxn},subopt_growth);
    %Set minimization of glucose uptake as objective function
    tempModel = setParam(tempModel,'obj',Csource_id,-1);
    solution  = solveLP(tempModel,1);
    if ~isempty(solution.x) 
        %Fix optimal glucose uptake and set protein usage as objective to
        %minimize
        tempModel = setParam(tempModel,'lb',Csource_id,0.999*solution.x(glucIndexEC));
        tempModel = setParam(tempModel,'ub',Csource_id,1.001*solution.x(glucIndexEC));
        %Set minimization of total protein usage as objective function
        tempModel = setParam(tempModel,'ub','prot_pool_exchange',1000);
        tempModel = setParam(tempModel,'obj','prot_pool_exchange',-1);
        %Solve!
        sol  = solveLP(tempModel,1);
        if ~isempty(sol.x)
            exchangeVector = sol.x(exchIndexes);
            newRow  = [subopt_growth, exchangeVector'];
            results = [results; newRow];
            if isempty(D_crit) && (sol.x(protIndex)>= 0.999*ecModel.ub(protIndex))
                D_crit = subopt_growth;
            end
        end
    end
    i = i+1;
end
%Plot results
names = {'Glucose' 'Oxygen' 'CO2' 'Ethanol' 'Acetate' 'prot_{pool}' 'D_{crit}'};
yyaxis left
plot(results(:,1),results(:,2),results(:,1),results(:,3),results(:,1),results(:,4),results(:,1),results(:,5),results(:,1),results(:,6),'LineWidth',3)
set(gca,'FontSize',14)
xlabel('Dilution rate [1/h]','FontSize',18)
ylabel('Exchange fluxes [mmol/gDw h]','FontSize',18)
xlim([0 max(results(:,1))])

yyaxis right
plot(results(:,1),results(:,7),'LineWidth',3)
ylabel('Protein pool usage [g_{prot}/gDw]','FontSize',18)

x1 = repelem(D_crit,100);
y1 = linspace(0,0.12,100);
hold on
plot(x1,y1,'LineWidth',1)
legend(names,'Location','northwest','FontSize',16)
hold off
cd ..