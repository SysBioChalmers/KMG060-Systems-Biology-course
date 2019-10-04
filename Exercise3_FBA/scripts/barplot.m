% Bar plot:
function barplot(x,pos,flux_names,var_names,title,length)
figure
bar(x(pos,:))
set(gca,'XTickLabel',flux_names)
ylabel(title)
if ~isempty(var_names)
    legend(var_names,'Location','northwest')
end
box on
if ~isempty(length)
    position    = get(gcf,'Position');
    position(3) = length;
    set(gcf,'Position',position)
end
set(gca, 'LooseInset', get(gca,'TightInset'))

end
