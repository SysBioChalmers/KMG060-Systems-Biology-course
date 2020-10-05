function model = minimal_Y6(model,Csource,cSourceUptake)
% change Y6 model media to minimal - ammonium, glucose, oxygen,
% phosphate, sulphate
% Bicarbonate production is blocked to get rid of conflict due to bicarbonate and carbon dioxide
% are regarded as equivalent
% the function is from:https://doi.org/10.1371/journal.pcbi.1004530

% start with a clean slate: set all exchange reactions to upper bound = 1000
% and lower bound = 0 (ie, unconstrained excretion, no uptake)

%Ivan Domenzain 0ctober 2020

exchangeRxns = findExcRxns(model);
model.lb(exchangeRxns) = 0;
model.ub(exchangeRxns) = 1000;

desiredExchanges = {'r_1654'; ... % ammonium exchange
                    'r_1992'; ... % oxygen exchange
                    'r_2005'; ... % phosphate exchange
                    'r_2060'; ... % sulphate exchange
                    'r_1861'; ... % iron exchange, for test of expanded biomass def
                    'r_1832'; ... % hydrogen exchange
                    'r_2100'; ... % water exchange
                    'r_4593'; ... % chloride exchange
                    'r_4595'; ... % Mn(2+) exchange
                    'r_4596'; ... % Zn(2+) exchange
                    'r_4597'; ... % Mg(2+) exchange
                    'r_2049'; ... % sodium exchange
                    'r_4594'; ... % Cu(2+) exchange
                    'r_4600'; ... % Ca(2+) exchange
                    'r_2020' };   % potassium exchange
              
                
model = setBounds(model,desiredExchanges,true);

blockedExchanges = {'r_1663'; ... % bicarbonate exchange
                    'r_4062'; ... % lipid backbone exchange
                    'r_4064'};    % lipid chain exchange
                
model = setBounds(model,blockedExchanges,false);

idx = find(strcmpi(model.rxns,Csource));     % Csource exchange
model = setBounds(model,model.rxns(idx),true,cSourceUptake);

end

function model = setBounds(model,exch_ids,allowUptk,flux)
if nargin<4
    flux = -1000;
end
for i = 1:length(exch_ids)
    rxnID = exch_ids{i};
    index = find(strcmpi(model.rxns,rxnID));
    if ~isempty(index)
        if allowUptk
            model.lb(index) = flux;
        else
            model.lb(index) = 0;
            model.ub(index) = 0;
        end
    else
        warning(['rxn ' rxnID ' was not found in the model'])
    end
end
end
