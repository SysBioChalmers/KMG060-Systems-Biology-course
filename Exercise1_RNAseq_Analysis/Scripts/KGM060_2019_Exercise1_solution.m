%% KMGO60 Systems Biology course
%  Computational exercise #1
%  Omics data analysis (RNAseq data)
%
%  Christoph B?rlin
%  Ivan Domenzain
%
%  Last edited. 2019-09-06

%% Step 1: Load RNAseq Data and normalize the data

%Load RNAseq data
rawCounts = readtable('../data/Saccharomyces_RNAseq_RAW_counts.csv');

%Create logical vectors for easy access of the different conditions
%now rawCounts(:,HiT) outputs only the three columns with the data for
%the high temperature stress condition.
samples = rawCounts(:,2:end).Properties.VariableNames;
ref     = startsWith(samples,'Control');
HiT     = startsWith(samples,'HighTemperature');
LpH     = startsWith(samples,'lowPH');
Osm     = startsWith(samples,'OsmoPressure');
anox    = startsWith(samples,'anaerobic');
%Get only the counts from the RNAseq data
rawCountsNum = rawCounts{:,2:end};

%Estimate pseudo-reference with geometric mean row by row
pseudoRefSample = geomean(rawCountsNum,2);
%Data normalization
%First, get the positions for all of the genes with non-zero reads
nz = pseudoRefSample > 0;
%Get ratios of expression, dividing each read by its corresponding gene
%average expression across samples
ratios = rawCountsNum(nz,:)./pseudoRefSample(nz);
%Get sizeFactors for normalization of each column (library of gene reads) 
sizeFactors = median(ratios,1);
% normalize the raw counts using the calculated normalization factors
normCounts = rawCountsNum./sizeFactors;
%Let's take a look to the normalization effects on the gene counts
%distributions. Generate boxplots for each of the samples distributions in
%both rawCountsNum and also normCounts
figure
subplot(1,2,1);
boxplot(rawCountsNum)
title('Raw counts')
ylabel('Read counts')
xlabel('Samples')

subplot(1,2,2);
boxplot(normCounts)
title('Normalized counts')
ylabel('Read counts')
xlabel('Samples')
%##################
%TASK:
%
%Can you extract any information from the boxplots? It might be the case
%that plotting counts in the linear scale is not the best visualization
%method, as you can see this dataset spans several orders of magnitude in
%gene read counts, so transforming data into a logarithmic scale might be
%more informative.
%
%Now you should get a new pair of boxplot graphs with the log2 transformed
%datasets. What can you say about the data normalization step, did it yield
%any evident effects?
figure
subplot(1,2,1);
boxplot(log2(rawCountsNum))
title('Raw counts')
ylabel('log_{2} read counts')
xlabel('Samples')

subplot(1,2,2);
boxplot(log2(normCounts))
title('Normalized counts')
ylabel('log_{2} read counts')
xlabel('Samples')
%##################
%% Step 2: PCA analysis

%Perform PCA
[pc, zscores, pcvars] = pca(normCounts');
%Get the explained variance for each component
explained = pcvars*100/sum(pcvars);
%Pick 5 colors to distinguish the 5 different conditions. Colormap returns
%the a mapped normalized RGB values for the specified amount of colorlevels
%in a given palette (jet in this case)
colorScheme = colormap(jet(5));
hold on
%Get a PCA plot showing PC1 and PC2 for all the samples in the 5 conditions
% **Note: PC's are stored in the zscores table in which for each sample 
%         (columns) the 15 first principal components are shown (rows).
condStr = {};
figure
for i=1:5
    %Get the condition name
    nameTMP = strsplit(rawCounts.Properties.VariableNames{(i-1)*3+2},'_');
    PC1     = zscores((i-1)*3+1:i*3,1);
    PC2     = zscores((i-1)*3+1:i*3,2);
    condStr = [condStr;nameTMP{1}];
    scatter(PC1,PC2,50,colorScheme(i,:),'fill');
end
%##################
%TASK:
%Add the amount of variance each principal component explains to the
%corresponding axis labels.
xlabel(['First Principal Component (' num2str(explained(1)) '% variance)']);
ylabel(['Second Principal Component (' num2str(explained(2)) '% variance)']);
title('Principal Component Scatter Plot');
legend(condStr)
hold off
%##################
%% Step 3: Test for differential gene expression

%pick one of the four stress condition that you would like to compare to
%the control
normCountsStress  = normCounts(:,HiT);
normCountsref     = normCounts(:,ref);

%Perform the test for differential expression using a negative binomial model
tLocal = nbintest(normCountsref,normCountsStress,'VarianceLink','LocalRegression');
% correct for multiple testing using the Benjamini Hochberg correction
padj = mafdr(tLocal.pValue,'BHFDR',true);

%create overview table of the results including the mean for both
%conditions, the resulting log2 fold change and the calculated pValue 
meanRef    = mean(normCountsref,2);
meanStress = mean(normCountsStress,2);
log2FC     = log2(meanRef./meanStress);
geneTable  = table(meanRef,meanStress,log2FC,tLocal.pValue,padj);
%Add row and column names
geneTable.Properties.RowNames = rawCounts.GeneName;
geneTable.Properties.VariableNames={'Mean_Ref','Mean_Stress','Log2_FC','pVal','adjPVal'};

%create a Volcano Plot for visual inspection of the data
colorVector            = zeros(length(meanRef),3);
indexDE                = geneTable.adjPVal<=0.01 & abs(geneTable.Log2_FC)>=2;
downReg                = geneTable.adjPVal<=0.01 & (geneTable.Log2_FC)<=-2;
upReg                  = geneTable.adjPVal<=0.01 & (geneTable.Log2_FC)>=2;
colorVector(indexDE,1) = 1;
scatter(geneTable.Log2_FC,-log10(geneTable.adjPVal),30,colorVector,'fill')
xlabel('Log_{2} fold-change');
ylabel('Abs(Log_{10}(adj.Pvalue))');
disp(['There are :' num2str(sum(indexDE)) ' DE genes in total'])
disp(['There are :' num2str(sum(downReg)) ' down-regulated genes in total'])
disp(['There are :' num2str(sum(upReg)) ' up-regulated genes in total'])
%##################
% Optional TASK: (1.5 extra points) 
% 1) How sensitive is your analysis to the chosen differential expression
% thresholds? Try to repeat the DE analysis for different pValues subject
% to constant log2FC thresholds and plot your results as: number of DE
% genes vs. pValue threshold. If you dare to get several curves for
% different log2FC values it would give you an idea on how both parameters
% are afecting your results.
pValIntervals = [0.001 0.005 0.01 0.05 0.1];
log2FCVect    = (1:(5-1)/100:5);
legendStr     = [];
for pV = pValIntervals
    DEvector = [];
    for log2Val = log2FCVect
        indexDE  = geneTable.adjPVal<=pV & abs(geneTable.Log2_FC)>=log2Val;
        numDE    = sum(indexDE);
        DEvector = [DEvector;numDE];
    end
    plot(log2FCVect,DEvector,'LineWidth',3)
    hold on
    legendStr = [legendStr, {['adj.P-val: ' num2str(pV)]}];
end
xlabel('log_{2}fold-change');
ylabel('# of DE genes');
legend(legendStr)
hold off
%##################
%% Step 4: First analysis of DE genes

%Manual inspection of genes with highest differential expression
geneTable  = sortrows(geneTable,'adjPVal','ascend');
top10Table = geneTable(1:10,:);
disp(top10Table)
%##################
%TASK: 
%Do you find any interesting genes or patterns in this top differentially expressed genes?
%Combine the geneTable with the geneDescription table (from Data/GeneDescriptions.csv) for getting a short
%description for each gene.
%
%Open GeneDescriptions file
geneDescriptions  = readtable('../data/GeneDescriptions.csv');
[~,iA,iB]         = intersect(top10Table.Row,geneDescriptions.GeneName);
top10Descriptions = table2cell(geneDescriptions(iB,3));
top10Table.Descriptions = top10Descriptions(iA);
disp(top10Table)
%##################
%% Step 5: Find Associated GO Terms

%load assocition of GO Terms and genes
GoTerms_table = readtable('../data/GoTermsMapping.txt','delimiter','\t');
GoTermsIDs    = unique(GoTerms_table.GoTerm);

%convert it to a map (the matlab variant of a dictionary) for easier access
%A map consists of unique key - value pairs and has the advantage that with
%the key you can easily access the corresponding value. If you have a
%map called "Universities" and you want to look up the value associated to 
%the Key "Chalmers" use Universities('Chalmers') and you would get the 
%data associated to Chalmers.
GoTermsGeneMap = containers.Map();
for i = 1:height(GoTerms_table)
    key = GoTerms_table{i,'GeneName'}{1};
    if isKey(GoTermsGeneMap,key)
        GoTermsGeneMap(key) = [GoTermsGeneMap(key),GoTerms_table{i,'GoTerm'}{1}];
    else
        GoTermsGeneMap(key) = {GoTerms_table{i,'GoTerm'}{1}};
    end
end
fprintf('Number of annotated genes related to functional process is %d.\n',GoTermsGeneMap.Count)
fprintf('Number of unique GO terms associated to annotated genes is %d.\n',numel(unique(GoTerms_table.GoTerm)))
fprintf('Number of gene-GO term associations is %d.\n',numel(GoTerms_table))

%Look up associated GO Terms for the top differentially expressed gene
selectedGene      = geneTable.Properties.RowNames{1};
associatedGoTerms = GoTermsGeneMap(selectedGene);
disp(['Associated GO Terms for gene ',selectedGene])
disp(associatedGoTerms)

%The GO Term IDs are not informative, therefore we need to load the
%descriptions for them
GO = geneont('File','../data/GoTerms.obo');

%Look up the details of the first GO Term
%The lookup function needs just the ID number (in number format) as input
GoTermID         = str2double(associatedGoTerms{1}(4:end));
GoTermName       = GO(GoTermID).Terms.Name;
GoTermDefinition = GO(GoTermID).Terms.Definition;
%output result
disp(['GO Term ',num2str(GoTermID),' - ',GoTermName,' : ',GoTermDefinition])

%##################
% TASK: 
% Create a loop to print the details for each of the associated GO Terms
% for the top1 differentially expressed gene.
%
% What can you learn from it? Compare the knowledge you got with the
% summary paragraph on the Saccharomyces Genome Database (yeastgenome.org)
%##################
for i=1:length(associatedGoTerms)
    GoTermID         = str2double(associatedGoTerms{i}(4:end));
    GoTermName       = GO(GoTermID).Terms.Name;
    GoTermDefinition = GO(GoTermID).Terms.Definition;
    %output result
    disp(['GO Term ',num2str(GoTermID),' - ',GoTermName,' : ',GoTermDefinition])
end
%% Step 6: Find enriched GO Terms in differentially expressed genes

%Select GO Term, in this case the top1 is chosen
rankedPos = 1;
GoTermID  = str2double(associatedGoTerms{rankedPos}(4:end));

%split gene list into DE and non DE genes
indexDE       = geneTable.adjPVal<=0.01 & abs(geneTable.Log2_FC)>=2;
geneListDE    = geneTable.Properties.RowNames(indexDE);
geneListNonDE = geneTable.Properties.RowNames(~indexDE);

%Check enrichment for selected GO Terms
%Count number of occurences for the GO Term in all genes and in all DE
%genes
GoTermCountAll = 0;
GoTermCountDE  = 0;
for key = GoTermsGeneMap.keys
    associatedGoTerms = (GoTermsGeneMap(key{1}));
    %GO terms full IDs are of the form GO:XXXXXXX
    full_GoTerm_Id = sprintf('GO:%07d',GoTermID);
    %Search the GO term ID in the associated GO terms cell array
    if any(strcmp(associatedGoTerms,full_GoTerm_Id))
        %
        if any(strcmp(geneListDE,key{1}))
            GoTermCountDE = GoTermCountDE+1;
        end
        GoTermCountAll = GoTermCountAll+1;
    end
end
    
disp(['GO Term ',(GoTermID),' is in ',num2str(GoTermCountDE),...
    ' out of ',num2str(numel(geneListDE)),' DE Genes and in total there are ',...
    num2str(GoTermCountAll),' occurences in all ',num2str(numel(GoTermsGeneMap.keys)),' Genes'])

%run hypergeometric test to assess significance
pHyperGeo=hygepdf(GoTermCountDE,numel(GoTermsGeneMap.keys),GoTermCountAll,numel(geneListDE));
disp(['The probablity for this is ',num2str(pHyperGeo)])

%##################
% TASK: 
% Create a loop that goes over all existing GoTerms (GoTermsIDs) and
% calculates the enrichment. Don't forget to use a correction for multiple
% testing, e.g. the Benjamini Hochberg Correction as for the DE detection
% (mafdr).
% Sort GO Terms by adjusted pValue and print out every one that has a
% pValue <=0.01
%##################

%Christoph's
GOterms_pValues = [];
keys            = GoTermsGeneMap.keys;
for i=1:length(GoTermsIDs) 
    %GO terms full IDs are of the form GO:XXXXXXX
    full_GoTerm_Id = GoTermsIDs{i};
    GoTermID       = GoTermsIDs{i}(4:end);
    %Check enrichment for the i-th GO term
    %Count number of occurences for the GO Term in all genes and in all DE
    %genes
    GoTermCountAll = 0;
    GoTermCountDE  = 0;
    %Loop through all genes in GOtermsGeneMap
    for j=1:length(keys)
        key =keys{j};
        %Get associated GO terms to gene
        associatedGoTerms = (GoTermsGeneMap(key));
        %Search the GO term ID in the associated GO terms cell array
        if any(ismember(full_GoTerm_Id,associatedGoTerms))
            if any(strcmp(geneListDE,key))
                GoTermCountDE = GoTermCountDE+1;
            end
            GoTermCountAll = GoTermCountAll+1;
        end
    end
    %Compute significance
    pHyperGeo       = hygepdf(GoTermCountDE,numel(GoTermsGeneMap.keys),GoTermCountAll,numel(geneListDE));
    GOterms_pValues = [GOterms_pValues; pHyperGeo];
end
GOTerms_Padj = mafdr(GOterms_pValues,'BHFDR',true);
t            = table(GoTermsIDs,GOterms_pValues,GOTerms_Padj);
t            = sortrows(t,'GOTerms_Padj','ascend');
close all
clc

% %%%%% Ivan's alternative
GOterms_pValues   = [];
GoTermsIDs_table  = GoTerms_table.GoTerm;
totalGenes        = numel(GoTermsGeneMap.keys);
DEnumber          = numel(geneListDE);
for i=1:length(GoTermsIDs)
    GoTermID = GoTermsIDs(i);
    %Check enrichment for the i-th GO term
    %Count number of occurences for the GO Term in all genes and in all DE
    %genes
    GoTermCountAll = 0;
    GoTermCountDE  = 0;
    presence       = find(strcmp(GoTermsIDs_table,GoTermID));
    if ~isempty(presence)
        relatedGenes  = GoTerms_table.GeneName(presence);
        [~,indexesDE] = intersect(geneListDE,relatedGenes);
        GoTermCountDE  = length(indexesDE);
        GoTermCountAll = length(presence);
    end
    pHyperGeo       = hygepdf(GoTermCountDE,totalGenes,GoTermCountAll,DEnumber);
    GOterms_pValues = [GOterms_pValues; pHyperGeo];
end
GOTerms_Padj = mafdr(GOterms_pValues,'BHFDR',true);
t            = table(GoTermsIDs,GOterms_pValues,GOTerms_Padj);
t            = sortrows(t,'GOTerms_Padj','ascend');
Enrichment_threshold = 0.10;
enrichedGOterms      = GoTermsIDs_table(find(GOTerms_Padj<=Enrichment_threshold));