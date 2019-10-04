% Find essential genes
function essential = findEssentialGenes(rxnGeneMat,c,S,b,LB,UB,genes)
essential = {};
obj = find(c);
n = length(c);
for i=1:n
    mutantLB = LB;
    mutantUB = UB;
    isoenzymes = sum(rxnGeneMat(i,:));
    %Block reactions that are encoded by a single gene
    if isoenzymes == 1
        mutantLB(i) = 0;
        mutantUB(i) = 0;
        v = maximize(c,S,b,mutantLB,mutantUB);
        %If biomass production is lower than a tolerance factor then the
        %gene is considered as essential
        if v(obj) < 1E-4
            essential = [essential;genes(i)];
        end
    end
end

end









