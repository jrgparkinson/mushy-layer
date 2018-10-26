function makeOptimalPlot(Ra_arr, CR_arr, y, fluxes, xVar, optimalExtrapolationFlag,...
    xlab, yLabel, makeLegend, xOffset, plotType)
%
% xVar = 0 (Ra) or 1 (CR)
RA = 0;
CR = 1;

PLOT_NORMAL = 0;
PLOT_SEMILOGX = 1;
PLOT_LOGLOG = 2;

if nargin < 9
    makeLegend = false;
end

if nargin < 11
    plotType = PLOT_NORMAL;
end




hold on;
leg = {};

if xVar == RA
    x_arr = Ra_arr;
    leg_arr = CR_arr;
    leg_var = CR;
else
    x_arr = CR_arr;
    leg_arr = Ra_arr;
    leg_var = RA;
end

if nargin < 10 || isempty(xOffset)
    xOffset = 0*leg_arr;
end

legendPlots = [];

for leg_i = 1:length(leg_arr)
    
    if isnan(leg_arr(leg_i) )
        continue;
    end
    
    if xVar == RA
        yForParam      = y(leg_i, :) ;
        fluxesForParam = fluxes(leg_i, :);
    else
        yForParam      = y(:, leg_i) ;
        fluxesForParam = fluxes(:, leg_i);
    end
    
    % Only plot values with valid fluxes
    valid= ~isnan(fluxesForParam);
    
    
    
    if sum(valid) > 0
        
        if xOffset(leg_i) ~= 0
            this_x_arr = (x_arr - xOffset(leg_i))/xOffset(leg_i);
        else
            this_x_arr = x_arr;
        end
        
        if plotType == PLOT_SEMILOGX
            legendPlots(end+1) = semilogx(this_x_arr(valid), yForParam(valid));
        elseif plotType == PLOT_LOGLOG
            legendPlots(end+1) = plot(log10(this_x_arr(valid)), log10(yForParam(valid)));
        else
            legendPlots(end+1) = plot(this_x_arr(valid), yForParam(valid));
        end
        
        
        
        % Also add different label for dodgy runs
        for x_arr_i = 1:length(this_x_arr)
            dodgyRun = false;
            
            if xVar == RA
                if optimalExtrapolationFlag(leg_i, x_arr_i) == 1
                    dodgyRun = true;
                end
            else
                if optimalExtrapolationFlag(x_arr_i, leg_i) == 1
                    dodgyRun = true;
                end
            end
            
            if dodgyRun
                
                if plotType == PLOT_SEMILOGX
                    semilogx(this_x_arr(x_arr_i), yForParam(x_arr_i), 'sk');
                else
                    plot(this_x_arr(x_arr_i), yForParam(x_arr_i), 'sk');
                end
                
            end
            
        end
        
        if leg_var == RA
            leg{end+1} = ['$Rm_S=', num2str(leg_arr(leg_i)),'$'];
        else
            leg{end+1} = ['$\mathcal{C}=', num2str(leg_arr(leg_i)),'$'];
        end
    end
end
hold off;
box on;

xlabel(xlab);


ylabel(yLabel, 'Rotation', 0.0);

if makeLegend
    legend(legendPlots, leg, 'Location', 'eastoutside');
end

if plotType == PLOT_SEMILOGX
    set(gca, 'xscale', 'log');
end


end