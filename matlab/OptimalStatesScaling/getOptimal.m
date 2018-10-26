function [optimalVal] = getOptimal(runs, optimalWidth, key)

% Just want to interpolate from runs.key to get the value at width =
% optimalWidth

vals = []; widths = []; 
for j = 1:length(runs)
    if isfield(runs(j), key) && ~isempty(runs(j).(key)) ...
            && ~isnan(runs(j).width) && ~isnan(runs(j).(key))
    vals(end+1) = runs(j).(key);
    widths(end+1) = runs(j).width;
    end
end

if length(widths) > 2
optimalVal = interp1(widths, vals, optimalWidth, 'PCHIP');
else
    optimalVal = NaN;
end

%Really, we only want to fit to the three largest flux values
% sortedVals = sort(vals(~isnan(vals)));
% 
% if length(sortedVals) > 3
%     largestVals = sortedVals(end-2:end);
%     
%     for i = 1:length(vals)
%         if ~ismember(vals(i),largestVals)
%             vals(i) = NaN; widths(i) = NaN;
%         end
%         
%     end
% end
% 
% validVals = logical(~isnan(vals));
% 
% %validFluxes = fluxes > 0
% % Fit quadratic to data. Ignore flux = 0
% % F = p(1)*L^2 + p(2)*L + p(3) = a*L^2 + b*L + c
% p = polyfit(widths(validVals), vals(validVals), 2);
% 
% optimalVal = p(3) - (p(2)^2) / (4*p(1));


end