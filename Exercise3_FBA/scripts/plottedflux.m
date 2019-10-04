% Plotted flux:
function [flux,plot] = plottedflux(name,mets,S,fluxes)
pos  = find(~cellfun(@isempty,strfind(mets,name)));
flux = zeros(size(fluxes));
for i = 1:length(pos)
    for j = 1:length(flux)
        flux(j,:) = flux(j,:) + fluxes(j,:).*S(pos(i),j);
    end
end
flux_abs = abs(flux);
flux_sum = sum(flux_abs,2);
plot     = flux_sum > 1e-3;

end