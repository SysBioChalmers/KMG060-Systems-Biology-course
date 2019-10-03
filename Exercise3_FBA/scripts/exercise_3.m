%% Exercise 3 solution
%
% Benjamin Sanchez
% Ivan Domenzain    2019-10-02

%% Problem 1: Finding a flux distribution
% 1.1 Setting up the model

% First, define all relevant variables:
load('../models/smallModel.mat'); % The model is stored in the 'model' subfolder
% To simplify later commands, define all relevant variables based on the
% information inside the model.
S   = model.S; % The S-matrix
b   = model.b; % Vector indicating which metabolites are boundary metabolites or not
c   = model.c; % Vector indicating which reaction is the objective function or not
LB  = model.lb; % Vector with lower bounds
UB  = model.ub; % Vector with upper bounds
REV = model.rev; % Vector indicating which reactions are reversible

%%
% A) The number of reactions can be obtained in different ways: the model
% field model.rxns contains the list of reactions, the length of this field
% equals then the number of reactions.
num_rxns = length(model.rxns); % Here we store the output in 'num_rxns'

% However, some of the vectors that we extracted from the model above would
% be the same length: 'c' corresponds to the reactions. So their length can
% also be used.

% Also, the size S-matrix is indicative of the number of reactions and
% metabolites.

% If you want MATLAB to show a neat message stating for instance the number
% reactions, you can do this with the code below. Notice that you first
% have to convert the numerical value stored in num_rxns to a character
% string using the num2str function.
disp(['Number of reactions: ' num2str(num_rxns)])

% B) Modify the commands above for metabolites instead of reactions.

% C) This has something to do with the degrees of freedom...

%%
% 1.2 Sanity and consistency checks
% A) How many reversible reactions are there in the model? 'REV' contains
% information about reversibility, where '1' stands for true (=reversible)
% and '0' stands for false (=irreversible). If you first find which
% reactions are reversible, then you can take the length of this result to
% get the number of reversible reactions:
rev_Rxns = find(REV); % Find out which reactions are noted as reversible
numb_rev_Rxns = length(rev_Rxns); % Find out how many reactions are noted as reversible
disp(['Number of reversible reactions: ' num2str(numb_rev_Rxns)])

% Mathematical confirmation: which structure in the model contains
% information about reversibility (besides REV [or model.rev])? LB and UB!
% If LB = 0 while REV = 1, there is something wrong: the reaction is
% supposed to be reversible, but the LB doesn't allow negative flux. Use
% this to find inconsistencies in the reversibility annotation:
inconsistencies = rev_Rxns(logical(REV(rev_Rxns).*LB(rev_Rxns))==0);
disp(['Reversibility inconsistency found in: ' model.rxnNames{inconsistencies}])

% B) Two of the glycolysis reactions have incorrect reversibility, which?
% Check the reversiblity in literature.
%
% C) Where in the model is reversibility information kept? Correct this for
% the reactions identified above.

%%
% 1.3 Determined problem:
% A) model.b is a vector where indicating in this case which metabolites are
% boundary metabolites. If b=0 for a certain metabolite, than this is not a
% boundary metabolite, but rather an internal metabolite. For internal
% metabolites, the steady-state assumption holds, where Sv=0. Inspecting b,
% it seems that none of the metabolite are boundary metabolites, b is an
% all-zero vector. Basically, Sv=b. 
% You can solve this for v using:
v = S\b

% B) Check the rank of the matrix, and compare this to the number of
% reactions. Can you recalculate the degrees of freedom and compare to what
% you got from 1.1C?


%%
% 1.4 Undetermined problem
% 
% A) We start by finding the positions of the relevant exchange fluxes.
% With strcmp() we can find strings in the vector of reaction names. 
pos(1) = find(strcmp(model.rxnNames,'Glucose exchange'));
pos(2) = find(strcmp(model.rxnNames,'O2 exchange'));
pos(3) = find(strcmp(model.rxnNames,'Biomass exchange'));
pos(4) = find(strcmp(model.rxnNames,'CO2 exchange'));
pos(5) = find(strcmp(model.rxnNames,'Acetate exchange'));
pos(6) = find(strcmp(model.rxnNames,'Ethanol exchange'));
pos(7) = find(strcmp(model.rxnNames,'Glycerol exchange'));
%%
% Both glucose and oxygen uptake rates are blocked.
disp(['glcEX: LB = ' num2str(LB(pos(1))) '; UB = ' num2str(UB(pos(1)))])
disp(['o2EX: LB = '  num2str(LB(pos(2))) '; UB = ' num2str(UB(pos(2)))])
% If you're not convinced by this, think about in which direction the
% reaction is defined (see exercise instructions)

% Change this by decreasing the lower bounds of each reaction:
LB(pos(1)) = -1; % This is for glucose. Modify this line to make the
% relevant change to O2 exchange
%%
% Optimize for growth and plot the exchange fluxes:
v = maximize(c,S,b,LB,UB); % A little function that runs linprog and
% inverts the solution, so that the fluxes are in the correct direction.
% We can specify labels for each bar in the plot:
flux_names = {'glcEX','o2EX','bioEX','co2EX','acEX','ethEX','glyEX'};
% Here we plot, using the data from v, we only plot those fluxes that are
% mentioned in pos, we also specify the labels and title. Two other
% parameters are not used yet [], but we'll get to that later..
barplot(v,pos,flux_names,[],'Exchange Fluxes',[])

% The next line is just to make sure that the next plots are not completely
% overlapping each other.
base = get(gcf,'Position'); base = base(3);

%%
% B) ATP/NADH/NADPH

% We will plot all related fluxes (multiplied by the corresponding
% stoichiometry) in barplots for simplifying the analysis. The following
% two functions have been defined to simplify the progress.

[flux_ATP,plot_ATP] = plottedflux('ATP',model.mets,S,v);
barplot(flux_ATP,plot_ATP,model.rxns(plot_ATP),[],'Normalized Fluxes',base+200)

% Repeat this for NADH and NADPH, by modifying the two commands above.

%% Problem 2: Testing objective functions
% 2.1 Fix glucose / Max growth rate

% The objective function should already be set to growth rate, while the
% glucose consumption should also be fixed. We can confirm this:
objFunc = find(c); % Find which reaction is the objective function
disp(['Reaction ' model.rxnNames{objFunc} ' is defined as objective function'])
LB(pos(1)) = -1; % Ensure that glucose consumption is -1.
va         = maximize(c,S,b,LB,UB); % We store the result in vector va, to be used for plotting later

% 2.2 Now you need to modify the objective function, in a similar way as
% changing UB or LB. Then, store the result from maximized acetate
% production in vector vb.
c          = zeros(num_rxns,1); % First, set c=0 for all reactions
c(pos(5))  = +1; % And only set c=1 for the reaction we want as objective function
vb         = maximize(c,S,b,LB,UB); % We store the result in vector va, to be used for plotting later

% 2.3 Change the objective function as above and store the result in vc.
vc          = maximize(c,S,b,LB,UB);

% The results from the three different objective functions are now plotted
% in a single graph, with the two lines below. Note that we use va, vb and
% vc as input.
var_names = {'Max growth','Max acetate','Max ATP maintenance'};
barplot([va,vb,vc],pos,flux_names,var_names,'Normalized Fluxes',base)

%% Problem 3: Changing growth conditions
% 3.1 Adding additional reactions

% a. Alternative aldehyde dehydrogenase

% The reaction to add is:
% Acetaldehyde + NAD -> Acetate + NADH
% First, we need to identify on which position these metabolites are located
pos_ACA  = strcmp('ACA_c',model.mets);
pos_NAD  = strcmp('NAD_c',model.mets);
pos_AC   = strcmp('AC_c',model.mets);
pos_NADH = strcmp('NADH_c',model.mets);
% Then we generate the new column for this reaction, that will be added to
% the S-matrix later on
ald_deh  = -pos_ACA + -pos_NAD + pos_AC + pos_NADH;

% b. Glyoxylate cycle
 
% The reactions to add are:
% Isocitrate -> Glyoxylate + Succinate
% Acetyl-CoA + Glyoxylate -> Malate + CoA

% Which can be lumped together in just one:
% Isocitrate + Acetyl-CoA -> Succinate + Malate + CoA
% Again, we need to identify the location of the metabolites
pos_ICIT  = strcmp('ICI_m',model.mets);
pos_ACCOA = strcmp('ACCOA_m',model.mets);
pos_SUCC  = strcmp('SUC_m',model.mets);
pos_MAL   = strcmp('MAL_m',model.mets);
pos_COA   = strcmp('COA_m',model.mets);
% And again generate the new column for the S-matrix
glyox_cyc = -pos_ICIT + -pos_ACCOA + pos_SUCC + pos_MAL + pos_COA;

% Updating model components
% 
% The vector b remains the same (no new metabolites) and we will changed
% the objective function later so we don't need to change the vector c yet.
% We just need to update model.rxns (with two extra names), S (with 2 extra
% columns), and LB and UB (with 2 extra rows):
model.rxns = [model.rxns;'ALD2';'GLYOX'];
S          = [S,ald_deh,glyox_cyc];
LB         = [LB;0;0];
UB         = [UB;1000;1000];

%%
% 3.2 Comparing growth conditions
 
% a. Glucose - Aerobic
% Glucose consumption and growth are already correctly constrained, so we
% only change the objective function to optimize growth:
c         = zeros(size(model.rxns));
c(pos(3)) = +1;
va        = maximize(c,S,b,LB,UB);

% b. Glucose - Anaerobic
% We repeat the simulation but with no oxygen consumption:
LB(pos(2)) = 0;
vb         = maximize(c,S,b,LB,UB);

% c. Ethanol - Aerobic
% We block glucose consumption, we unblock oxygen consumption, & we fix the
% ethanol consumption to be the same molar amount of carbon as the previous
% simulations:
LB(pos(1)) = 0;
LB(pos(2)) = -1000;
LB(pos(6)) = ; % Fill in the correct molar amount of ethanol!
vc         = maximize(c,S,b,LB,UB);

% d. Ethanol - Anaerobic
% We repeat the simulation but with no oxygen consumption:
LB(pos(2)) = 0;
vd         = maximize(c,S,b,LB,UB);

% Plot exchange fluxes
var_names  = {'Gluc-Aerobic','Gluc-Anaerobic','EtOH-Aerobic','EtOH-Anaerobic'};
barplot([va,vb,vc,vd],pos,flux_names,var_names,'Fluxes',base+100)

%% 4 Additional exercise In silico gene deletions

rxnGeneMat = model.rxnGeneMat;
grRules    = model.grRules;
%Introduce missing elements for the new reactions
[~,genes]  = size(rxnGeneMat);
rxnGeneMat = [rxnGeneMat; zeros(1,genes); zeros(1,genes)];
grRules    = [grRules; {''}; {''}];
% Allow glucose and oxygen uptakes and block ethanol consumption
LB(pos(1)) = -1;
LB(pos(2)) = -10;
LB(pos(6)) = 0;

% 4.1) Find essential genes
essential = findEssentialGenes(rxnGeneMat,c,S,b,LB,UB,grRules);

disp('Essential genes for aerobic growth on glucose: ')
disp(essential)