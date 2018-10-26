function [p] = fitFluxes(runs)

fluxes = NaN*ones(1, length(runs)); widths = NaN*ones(1, length(runs));
for j = 1:length(runs)
    fluxes(j) = runs(j).flux;
    widths(j) = runs(j).width;
end

%Really, we only want to fit to the three largest flux values
sortedFlux = sort(fluxes(~isnan(fluxes)));

if length(sortedFlux) > 3
    largestFluxes = sortedFlux(end-2:end);
    
    for i = 1:length(fluxes)
        if ~ismember(fluxes(i),largestFluxes)
            fluxes(i) = NaN; widths(i) = NaN;
        end
        
    end
end

validFluxes = logical((fluxes > 0) .* (~isnan(fluxes)));

%validFluxes = fluxes > 0
% Fit quadratic to data. Ignore flux = 0
% F = p(1)*L^2 + p(2)*L + p(3) = a*L^2 + b*L + c
p = polyfit(widths(validFluxes), fluxes(validFluxes), 2);



end