function [fluxTable,enzTable_abs,enzTable_rel] = compare_ecModels_fluxDist(models,outputName)
if nargin<2
    outputName = [];
end        
enzTable_rel = [];
fluxTable    = [];
enzTable_abs = [];
for i = 1:length(models)
    model = models{i};
    fileName = ['../../models/prot_constrained/ecYeastGEM_' model '.mat'];
    if isfile(fileName)
        load(fileName)
        %Get enzymes 
        if i==1
            overlapProts = ecModelP.enzymes;
        end
        [rxnsTable,absUsages,relUsages] = get_fluxDist_table(ecModelP);
        if isempty(fluxTable)
            fluxTable = table(rxnsTable.rxns,rxnsTable.rxnNames,rxnsTable.formulas,'VariableNames',{'rxns' 'rxnNames' 'formulas'});
        end
        if isempty(enzTable_abs)
            enzTable_abs = table(absUsages.enzymes,'VariableNames',{'enzymes'});
        end
        
        [iA,iB] = ismember(overlapProts,relUsages.enzymes);
        overlapProts = overlapProts(iA);
        if isempty(enzTable_rel)
            enzTable_rel = table(relUsages.enzymes(iB),'VariableNames',{'enzymes'});
        end
        eval(['fluxTable.flux_' model ' = rxnsTable.flux;'])
        eval(['enzTable_abs.usage_' model ' = absUsages.abs_usage;'])
        eval(['enzTable_rel.usage_' model ' = relUsages.rel_usage(iB);'])
    end
end
newDir = '../../results/proteomics_integration';
mkdir(newDir)
fluxTable.grRules    = rxnsTable.grRules;
fluxTable.subSystems = rxnsTable.subSystems;
enzTable_abs.genes      = absUsages.genes;
enzTable_abs.shortNames = absUsages.shortNames;
enzTable_abs.subSystems = absUsages.subSystems;
enzTable_rel.genes      = relUsages.genes(iB);
enzTable_rel.shortNames = relUsages.shortNames(iB);
enzTable_rel.subSystems = relUsages.subSystems(iB);
if isempty(outputName)
    newFile = [newDir '/' outputName '_rxnsTable.txt'];
    writetable(fluxTable,newFile,'Delimiter', '\t','QuoteStrings',false);
    newFile = [newDir '/' outputName '_absEnz_Table.txt'];
    writetable(enzTable_abs,newFile,'Delimiter', '\t','QuoteStrings',false);
    newFile = [newDir '/' outputName '_relEnz_Table.txt'];
    writetable(enzTable_rel,newFile,'Delimiter', '\t','QuoteStrings',false);
end
end