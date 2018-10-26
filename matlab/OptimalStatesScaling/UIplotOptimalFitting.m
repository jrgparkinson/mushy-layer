% plot optimal fitting
% an interactive plot showing how each property has been fitted to obtain
% the optimal values, for each CR/Ra pair
function UIplotOptimalFitting

clear all;
close all;

set(groot, 'DefaultTextInterpreter', 'latex')
set(groot, 'DefaultLegendInterpreter', 'latex')


data_dir = ['/home/parkinsonjl/mnt/sharedStorage/', ...
    'optimalStates-highRes-newSoluteFlux/'];
data_dir = ['/home/parkinsonjl/mnt/sharedStorage/', ...
    'optimalStates-highRes-Wells/'];

data_dir = ['/home/parkinsonjl/mnt/sharedStorage/', ...
    'optimalStates-highRes-new/'];

allStatesFile = 'allStates.mat';

optimalValsFile = 'optimalVals.mat';


allStates = [];
leg = {};
legPlots = [];
optimalVals = [];

getAllStates([data_dir, allStatesFile]);
getOptimalVals([data_dir, optimalValsFile]);

% Get possible CR and Ra
k = keys(allStates);
CR_arr = [];
Ra_arr = [];
for CR_i = 1:length(k)
    keyPair = k(CR_i);
    keyPair = keyPair{1};
    temp = keyPair(1);  thisC  = temp{1};
    temp = keyPair(2);  thisRa = temp{1};
    
    
    if  ~ismember(thisC, CR_arr)
        CR_arr(end+1) = thisC;
    end
    
    if ~ismember(thisRa, Ra_arr)
        Ra_arr(end+1) = thisRa;
    end
    
end

Ra_arr = sort(Ra_arr); CR_arr = sort(CR_arr);

%global Ra CR fieldToPlot;
initial_Rai = round(length(Ra_arr)/2);
initial_CRi = round(length(CR_arr)/2);



Ra = Ra_arr(initial_Rai);
CR = CR_arr(initial_CRi);

keyPair = k(1);
keyPair = keyPair{1};
temp = keyPair(1);  thisC  = temp{1};
temp = keyPair(2);  thisRa = temp{1};

a_state = allStates(thisC, thisRa);
fields = sort(fieldnames(a_state));
fieldToPlot = 'flux';
initial_fieldi = find(strcmp(fields, fieldToPlot));

% Make figure
f = figure(1);
set(f, 'Position', [100 200 1600 800]);
ax = axes('Parent', f, 'position', [0.35 0.15 0.6 0.7]);


popupRaPos = [100 700 100 50];
popupRa = uicontrol('Style', 'popup',...
    'String', num2cell(Ra_arr),...
    'Position', popupRaPos,...
    'Callback', @setRa);
popupRa.Value = initial_Rai;

txtRa = uicontrol('Style','text',...
    'Position',popupRaPos - [50 0 70 0],...
    'String','Ra:'); %$\mathcal{R}$

popupCRPos =  popupRaPos - [0 50 0 0];
popupCR = uicontrol('Style', 'popup',...
    'String', num2cell(CR_arr),...
    'Position',popupCRPos,...
    'Callback', @setCR);
popupCR.Value = initial_CRi;

txtCR = uicontrol('Style','text',...
    'Position',popupCRPos - [50 0 70 0],...
    'String','CR:'); %$\mathcal{C}$:

popupFieldPos = popupCRPos - [0 50 -100 0];
popupField = uicontrol('Style', 'popup',...
    'String',  fields,...
    'Position', popupFieldPos,...
    'Callback', @setField);
popupField.Value = initial_fieldi;
txtField = uicontrol('Style','text',...
    'Position',popupFieldPos - [120 0 80 0],...
    'String','Field:'); %$\mathcal{C}$:


clearButton = uicontrol('Style', 'pushbutton', 'String', 'Clear Plot',...
    'Position', popupFieldPos - [0 50 0 0],...
    'Callback', @clearPlotsBtn);


% Make initial plot with initial values
clearPlots();
updatePlot();

    function setRa(source, event)
        Ra = str2num(source.String{source.Value});
        fprintf('Ra: %d \n', Ra);
        updatePlot();
    end

    function setCR(source, event)
        CR = str2num(source.String{source.Value});
        fprintf('CR: %1.3f \n', CR);
        updatePlot();
    end

    function setField(source, event)
        fieldToPlot = source.String{source.Value};
        fprintf('Field to plot: %s \n', fieldToPlot);
        clearPlots();
        updatePlot();
    end

    function updatePlot()
        title(['Ra = ', num2str(Ra), ', CR = ', num2str(CR), ', ', fieldToPlot]);
        
        ylabel(fieldToPlot);
        
        % Get data
        if ~isKey(allStates, CR, Ra)
            k = keys(allStates);
            for i=1:length(k)
                thisK = k{i};
               fprintf('%f, %d \n', thisK{1}(1), thisK{2}(1)); 
            end
            fprintf('No data for (CR, Ra) pair \n');
            return;
        end
        runs = allStates(CR, Ra);
        
        % Fit to fluxes
        % Get the fluxes
        fluxes = NaN*ones(1, length(runs));
        widths = NaN*ones(1, length(runs));
        vals = NaN*ones(1, length(runs));
        
        for j = 1:length(runs)
            fluxes(j) = runs(j).flux;
            widths(j) = runs(j).fullWidth;
            vals(j) = runs(j).(fieldToPlot);
        end
        
        % At this point, we convert flux values
        % This should now already be done in postProcess
        %fluxes = -1 - fluxes;
        
        %Really, we only want to fit to the three largest flux values
        sortedFlux = sort(fluxes(~isnan(fluxes)));
        
        largestFluxesWidths = widths;
        largestFluxes = fluxes;
        largestVals = vals;
        
        % Fit up to 5 largest values
        %numLargestFluxes = 5;
        numLargestFluxes = min(length(sortedFlux), 5);
        
        if numLargestFluxes > 2
            largeFlux = sortedFlux(end-(numLargestFluxes-1):end);
            
            for f_i = 1:length(fluxes)
                if ~ismember(fluxes(f_i),largeFlux)
                    largestFluxes(f_i) = NaN;
                    largestFluxesWidths(f_i) = NaN;
                    largestVals(f_i) = NaN;
                end
                
            end
            
        else
            fprintf('Not enough runs');
            return;
        end
        
        
        %validFluxes = logical((fluxes > 0) .* (~isnan(fluxes)));
        % Let fluxes be < 0 for now
        validFluxes  = logical((~isnan(fluxes)));
        
        validLargeFluxes = logical((~isnan(largestFluxes)));
        
        % If we don't have 3 points, skip
        if sum(validFluxes) < 3
            fprintf('Skipping, less than three flux values \n');
            return;
        end
        
        
        %validFluxes = fluxes > 0
        % Fit quadratic to data. Ignore flux = 0
        % F = p(1)*L^2 + p(2)*L + p(3) = a*L^2 + b*L + c
        p = polyfit(largestFluxesWidths(validLargeFluxes), ...
            largestFluxes(validLargeFluxes), 2);
        
        
        %pVal = polyfit(largestFluxesWidths(twoLargestVals), ...
        %    largestVals(twoLargestVals), 1);
        
        
        % The max value of the polynomial is given by c - b^2/(4*a)
        optimalFlux = p(3) - (p(2)^2) / (4*p(1));
        
        % The position of the maximum is -b/(2a)
        optimalFullWidth = - p(2) / (2*p(1));
        
        % Apply a few checks to this
        if optimalFlux < 0 || optimalFlux > 50  ||...
                optimalFullWidth < 0 || optimalFullWidth > 10.0
            optimalFlux = NaN;
            optimalFullWidth = NaN;
        end
        
        
        
        
        hold on;
        
        x = linspace(min(widths), max(widths));
        currentXLim = xlim;
        fprintf('Current x lim (%1.2f, %1.2f) \n', currentXLim(1), currentXLim(2));
        maxWidth = max(max(x, currentXLim(2)));
        minWidth = min(min(x, currentXLim(1)));
        fprintf('New x lim (%1.2f, %1.2f) \n', minWidth, maxWidth);
        
        plot(widths(validFluxes), vals(validFluxes), 'xk');
        %leg{end+1} = ['$\mathcal{C}=', num2str(CR), ', \mathcal{R}=',num2str(Ra), '$'];
        pltFit = plot(largestFluxesWidths, largestVals, 'db');
        %leg{end+1} = ['$\mathcal{C}=', num2str(CR), ', \mathcal{R}=',num2str(Ra), '$ (max flux)'];
        if length(widths) > 2
            if strcmp(fieldToPlot, 'flux')
                
                fittedFluxes = polyval(p,x);
                legPlots(end+1) =  plot(x, fittedFluxes, '-');
                
            else
                validWidthsVals = (~isnan(widths)).*(~isnan(vals));
                
                if sum(validWidthsVals) >= 2
                    
                    % Linear fit seems safest for now
                    optimalVal = interp1(widths(validWidthsVals==1), ...
                        vals(validWidthsVals==1), ...
                     optimalFullWidth, 'PCHIP');
                    
                     plt =  plot(x, x.*0 + optimalVal, '-');
                    
                    col = get( plt, 'Color');
                    legPlots(end+1) = plt;
                    plt2 = plot([optimalFullWidth, optimalFullWidth], [min(vals) max(vals)], ['--']);
                    set(plt2, 'Color', col);
                     set(pltFit, 'Color', col);
                end
            end
            
        end
        
       
        
        leg{end+1} = ['$\mathcal{C}=', num2str(CR), ', \mathcal{R}=',num2str(Ra), '$ (fit)'];
        
        hold off;
        
        box on;
        
        xlim([minWidth*0.9, maxWidth*1.1]);
        
        if strcmp(fieldToPlot, 'flux')
            maxFluxMeasured = max(fluxes(validFluxes));
            minFluxMeasured = min(fluxes(validFluxes));
            
            fluxUpperPlotLimit = NaN;
            fluxLowerPlotLimit = NaN;
            if ~isnan(maxFluxMeasured) && maxFluxMeasured > 0.0
                fluxUpperPlotLimit = maxFluxMeasured*1.1;
            end
            
            if ~isnan(fluxLowerPlotLimit) && fluxLowerPlotLimit > 0.0
                fluxLowerPlotLimit = fluxLowerPlotLimit*0.9;
            end
            
            if ~isnan(fluxLowerPlotLimit) && ~isnan(fluxUpperPlotLimit)
                ylim([fluxLowerPlotLimitfluxUpperPlotLimit]);
            end
            
        end
        
        
        xlabel('$\lambda$');
        box on;
        
        legend(legPlots, leg, 'Location', 'eastoutside');
        
    end

    function clearPlotsBtn(source, event)
        clearPlots()
    end

    function clearPlots()
        cla reset;
        leg = {}; legPlots = [];
        xlabel('$\lambda$');
        box on;
        xlim([0.2 0.3]);
%        legend(legPlots, leg, 'Location', 'eastoutside');
        
    end

    function getAllStates(file)
        load(file);
    end

    function getOptimalVals(file)
        load(file);
    end


end