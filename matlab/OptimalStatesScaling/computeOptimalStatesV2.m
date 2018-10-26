% Load all optimal states
% Indexing is optimalStates(CR, Ra) = array of postProcessed variables

% computeOptimalStatesV2('/home/parkinsonjl/mnt/sharedStorage/optimalStates-Brinkman/')
function computeOptimalStatesV2(output_folder, fixedLe)

close all;

if nargin < 2
    fixedLe = -1;
end

if nargin < 1
    output_folder = 'optimalStates-Brinkman-lowRes';
end


if output_folder(end) ~= '/'
    output_folder = [output_folder, '/'];
end

data_dir = getDataDir(output_folder);

optimalValsFile = 'optimalVals.mat';
allStatesFile = 'allStates.mat';
optimalFoldersFile = 'optimalFolders.csv';

folderEnding = 'steady';
permeabilityType = 'ChiCubedPermeability';

if strcmp(output_folder, 'optimalStates-Brinkman')
    permeabilityType = 'KozenyPermeability';
end

figure_output_dir = data_dir;


printQuality = '-r50';
% Finished setting options for script

noSteadyState = {};

% For each CR, Ra pair, we store the struct containing all postProcessed
% diagnostics
allStates = MapN();

% Add files in directory below to the optimal states
%folders = dir(data_dir);

folders = getAllOutputFolders(data_dir);

optimalFolders = {};

baseStruct = struct('channelHeight',NaN, ...
    'largeLScale',NaN, ...
    'smallLScale',NaN, ...
    'h',NaN, ...
    'H',NaN, ...
    'hPsi',NaN, ...
    'hChi',NaN,...
    'hporosity',NaN,...
    'Hporosity',NaN,...
    'channelWidth',NaN,...
    'heatAdvectionMush',NaN, ...
    'latentHeatMush',NaN, ...
    'TFrameAdvectionMush',NaN, ...
    'heatDiffusionMush',NaN, ...
    'saltAdvectionMush',NaN, ...
    'liquidSalinityFrameMush',NaN, ...
    'solidSalinityFrameMush',NaN, ...
    'saltDiffusionMush',NaN, ...
    'baroclinicTorqueMush',NaN, ...
    'vorticityPermeabilityMush',NaN, ...
    'vorticityDiffusionMush',NaN, ...
    'maxPorosityMush',NaN, ...
    'maxStreamfunctionMush',NaN, ...
    'maxPermeabilityMush',NaN, ...
    'avPorosityMush',NaN, ...
    'avStreamfunctionMush',NaN, ...
    'avPermeabilityMush',NaN, ...
    'heatAdvectionInterior',NaN, ...
    'latentHeatInterior',NaN, ...
    'TFrameAdvectionInterior',NaN, ...
    'heatDiffusionInterior',NaN, ...
    'saltAdvectionInterior',NaN, ...
    'liquidSalinityFrameInterior',NaN, ...
    'solidSalinityFrameInterior',NaN, ...
    'saltDiffusionInterior',NaN, ...
    'baroclinicTorqueInterior',NaN, ...
    'vorticityPermeabilityInterior',NaN, ...
    'vorticityDiffusionInterior',NaN, ...
    'maxPorosityInterior',NaN, ...
    'maxStreamfunctionInterior',NaN, ...
    'maxPermeabilityInterior',NaN, ...
    'avPorosityInterior',NaN, ...
    'avStreamfunctionInterior',NaN, ...
    'avPermeabilityInterior',NaN, ...
    'flux',NaN, ...
    'L2FsVertDiffusion',NaN, ...
    'L2FsVertFluid',NaN, ...
    'L2FsVertFrame',NaN, ...
    'L1FsVertDiffusion',NaN, ...
    'L1FsVertFluid',NaN, ...
    'L1FsVertFrame',NaN, ...
    'L0FsVertDiffusion',NaN, ...
    'L0FsVertFluid',NaN, ...
    'L0FsVertFrame',NaN, ...
    'fullWidth',NaN, ...
    'cellsWide',NaN, ...
    'folder',NaN);

%folders_regex = 'CR(\d+\.?\d?)RaC(\d+)Le(\d+)(\s+)pts.*';
%folders_regex = 'CR(\d+\.?\d+?)RaC(\d+\.?\d+?)Le(\d+)(\w+)pts(\d+)-steady'; %Make sure it ends in 0

for i = 1:length(folders)
    folder = folders{i};
    
    try
        
    try
        inputs = readInputs([folder, 'inputs']);
        St = str2double(inputs.stefan);
        frameAdv = str2double(inputs.nonDimVel);

        CR = str2double(inputs.compositionRatio);
        Ra = str2double(inputs.rayleighComp);
        Le = str2double(inputs.lewis);
        %PermType = str2num(inputs.permeabilityFunction);
        num_cells = inputs.num_cells;
        p = strsplit(num_cells, ' ');
        width = str2num(p{1});
   
    catch e %e is an MException struct
        fprintf('Could not get parameters from inputs file, skipping \n');
        fprintf(1,'The identifier was:\n%s \n',e.identifier);
        fprintf(1,'There was an error! The message was:\n%s \n',e.message);

        continue
    end
    
        % Only choose le=200
        if fixedLe > 0 && Le ~= fixedLe
            continue
        end
        
        CR = round(CR, 3); % Round to three decimal places
        
        fprintf('Importing %s \n', folder);
        
        %postProcess = dlmread([data_dir, '/', folder.name, '/postProcess.csv'], ',', 1, 0);
        postProcessfilename = [folder, '/postProcess.csv'];
        
        %tODO KEEP REWRITING THIS
        
        if ~exist(postProcessfilename, 'file')
            fprintf('  No postProcess.csv \n');
            continue
        end
        
        A = importdata(postProcessfilename, ',', 1);
        
        postProcessStruct = baseStruct;
        
        newStructs = cell2struct(num2cell(A.data), A.colheaders, 2);
        
        fnames = fieldnames(newStructs);
        for fieldname_i = 1:length(fnames)
            a =  fnames(fieldname_i);
            thisFieldName = a{1};
            if isfield(postProcessStruct, thisFieldName)
                postProcessStruct.(thisFieldName) = newStructs.(thisFieldName);
            else
               fprintf('!!! %s is not a field in the base struct', thisFieldName);
            end
        end
        
        postProcessStruct.cellsWide = width;
        postProcessStruct.folder = folder;
        
        fprintf('  Fs = %1.10f \n', postProcessStruct.flux);
        
        % Most recent post processed files have this field, skip the others
        test = find(strcmp(fieldnames(postProcessStruct), 'h'));
        if length(test) == 0
            continue
        end
        
        fprintf('CR=%1.3f, Ra=%d, ', CR, Ra);
        
        if isKey(allStates, CR, Ra)
            fprintf('key exists \n');
            thisState = allStates(CR, Ra);
            
            try 
                thisState(end+1) = postProcessStruct;
            catch E
                theseFields = fieldnames(postProcessStruct);
                existingFields = fieldnames(thisState(1));
                
                help = 0;
            end
            allStates(CR, Ra) = thisState;
        else
            fprintf('key does not exist \n');
            allStates(CR, Ra) = [postProcessStruct];
        end % end if  key exists
        
        
        
   catch e %e is an MException struct
        fprintf('Issue with processing folder \n');
        fprintf(1,'The identifier was:\n%s \n',e.identifier);
        fprintf(1,'There was an error! The message was:\n%s \n',e.message);

        continue
        
    end
        
end % end looop over folders
    



% Get all concentration ratios and rayleigh numbers
k = keys(allStates);
CR_arr = [];
Ra_arr = [];
for i = 1:length(k)
    keyPair = k(i);
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

% Now, for each Ra, C pair, use fluxes to find optimal state then
% interpolate to find optimal values of all other variables.

% Store everything in a new MapN(), where for each CR, Ra pair we have a
% struct containing all the optimal vals
optimalVals = MapN();

% Set up to plot each flux fit
m = 3; % Num rows
n = 3; % Num cols

subPlotsPerWindow = m*n;
window_i = 1;
subplot_i = 1;


for i = 1:length(k)
    key = k(i); key = key{1};
    CR = key{1}; Ra = key{2};
    runs = allStates(CR, Ra);
    fprintf(['CR = ', num2str(CR), ', Ra = ', num2str(Ra), ', files = ',num2str(length(runs)),'\n']);
    
    % Copy the first structure and replace values with NaN
    % will then put all the optimal values in here
    optimal = runs(1);
    
    % Set every element of optimalVals to NaN for now
    fields = fieldnames(optimal);
    for j = 1:numel(fields)
        optimal.(fields{j}) = NaN;
    end
    
    % Now we've got all the fluxes for this (CR, Ra), try and calculate the
    % optimal flux and width
    
    % Need at least three points to fit a quadratic
    if length(runs) > 2
        
        % Get the fluxes
        fluxes = NaN*ones(1, length(runs));
        widths = NaN*ones(1, length(runs));
        cells = NaN*ones(1, length(runs));
        runFolders = {};
        
        for j = 1:length(runs)
            fluxes(j) = runs(j).flux;
            widths(j) = runs(j).fullWidth;
            cells(j) = runs(j).cellsWide;
            runFolders{j} = runs(j).folder;
        end
        
        if sum(isnan(widths)) == length(widths)
            widths = cells;
        end
        
        % At this point, we convert flux values
        % This should now already be done in postProcess
        %fluxes = -1 - fluxes;
        
        %Really, we only want to fit to the three largest flux values
        sortedFlux = sort(fluxes(~isnan(fluxes)));
        
        largestFluxesWidths = widths;
        largestFluxes = fluxes;
        
        smallestFluxes = fluxes;
        smallestFluxesWidths = widths;
        
        smallFlux = [];
        
        smallFluxCounted = 0;
        nonNanFluxes = fluxes(~isnan(fluxes));
        
        sortedWidths = sort(widths);
        for w_i = 1:length(sortedWidths)
            % Get the index of this width in the flux, width arrays
            thisWidth = sortedWidths(w_i);
            index = find(widths == thisWidth, 1);
            thisF = fluxes(index);
            if length(thisF) > 0 
                if ~isnan(thisF)  && thisF > 0
                    if smallFluxCounted < 2
                smallFlux(end+1) = thisF;
                smallFluxCounted = smallFluxCounted + 1;
                    end
                end
            end
            
            
        end
        
        % Fit up to 5 largest values
        %numLargestFluxes = 5;
        numLargestFluxes = min(length(sortedFlux), 5);
        
        if numLargestFluxes > 2
            largeFlux = sortedFlux(end-(numLargestFluxes-1):end);
            
            for f_i = 1:length(fluxes)
                if ~ismember(fluxes(f_i),largeFlux)
                    largestFluxes(f_i) = NaN;
                    largestFluxesWidths(f_i) = NaN;
                end
                
                if ~ismember(fluxes(f_i),smallFlux)
                    smallestFluxes(f_i) = NaN;
                    smallestFluxesWidths(f_i) = NaN;
                end
                
            end
            
        else
            continue;
        end
        
        
        %validFluxes = logical((fluxes > 0) .* (~isnan(fluxes)));
        % Let fluxes be < 0 for now
        validFluxes  = logical((~isnan(fluxes)));
        
        validLargeFluxes = logical((~isnan(largestFluxes)));
        validSmallFluxes = logical((~isnan(smallestFluxes)));
        
        
        % Print out a list of the folders which we use to calculate optimal
        % values. If we need to re process simulations, these are the ones
        % we should do.
        validFolders = runFolders(validLargeFluxes);
        
        for folder_i=1:length(validFolders)
            optimalFolders{end+1} = validFolders{folder_i};
        end
        
        % If we don't have 3 points, skip
        if sum(validFluxes) < 3
            fprintf('Skipping, less than three flux values \n');
            continue
        end
        
        %validFluxes = fluxes > 0
        % Fit quadratic to data. Ignore flux = 0
        % F = p(1)*L^2 + p(2)*L + p(3) = a*L^2 + b*L + c
        p = polyfit(largestFluxesWidths(validLargeFluxes), ...
            largestFluxes(validLargeFluxes), 2);
        
        
        % The max value of the polynomial is given by c - b^2/(4*a)
        optimal.flux = p(3) - (p(2)^2) / (4*p(1));
        
        % The position of the maximum is -b/(2a)
        optimal.fullWidth = - p(2) / (2*p(1));
        
        % Apply a few checks to this
        if optimal.flux < 0 || optimal.flux > 50  ||...
                optimal.fullWidth < 0
            optimal.flux = NaN;
            optimal.fullWidth = NaN;
        end
        
        fprintf('Optimal flux: %1.5f \n', optimal.flux);
        
        
        % check if optimal width falls within bounds of runs we've done
        widthsUsed = largestFluxesWidths(validLargeFluxes);
        if optimal.fullWidth < min(widthsUsed) || optimal.fullWidth > max(widthsUsed)
            optimal.extrapolatedMax = true;
        else
            optimal.extrapolatedMax = false;
        end
        
        % Also find critical width
        % Do this by linear extrapolation from smallest two (valid) fluxes
        if length(smallestFluxesWidths(validSmallFluxes)) > 1
            
            pCrit = polyfit(smallestFluxesWidths(validSmallFluxes), ...
                smallestFluxes(validSmallFluxes), 1);
            
            % Find width where flux = 0
            optimal.criticalWidth = -pCrit(2)/pCrit(1);
            
            if optimal.criticalWidth < 0 || ...
                    optimal.criticalWidth > min(smallestFluxesWidths(validSmallFluxes))
                optimal.criticalWidth = min(smallestFluxesWidths(validSmallFluxes) );
                optimal.flagCritWidth = true;
            else
                optimal.flagCritWidth = false;
            end
            
        elseif length(smallestFluxesWidths(validSmallFluxes)) == 1
            optimal.criticalWidth =  min(smallestFluxesWidths(validSmallFluxes) );
            optimal.flagCritWidth = true;
            
        else
            optimal.criticalWidth = 0.0;
            optimal.flagCritWidth = true;
        end
        
        fprintf('L_{crit} = %1.5f \n', optimal.criticalWidth);
        
        % Let's make some plots to see what's going on
        % First check if we need to make a new window
        if subplot_i > subPlotsPerWindow || ...
                subplot_i == 1 && window_i == 1
            if  window_i ~= 1
                % Save figure
                set(h(window_i),'Units','Inches');
                pos = get(h(window_i),'Position');
                set(h(window_i),'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
                print(gcf, '-dpdf', [figure_output_dir, 'optStates',num2str(window_i),'.pdf'], printQuality);
            end
            
            
            subplot_i = 1;
            window_i = window_i + 1;
            
            
            
            h(window_i) = figure(); %set(gcf, 'Visible', 'off');
            set(h(window_i), 'Position', [100, 100, 1400, 800]);
        end
        
        subplot(m, n, subplot_i);
        hold on;
        
        plot(widths(validFluxes), fluxes(validFluxes), 'xk');
        plot(largestFluxesWidths, largestFluxes, 'db');
        if length(widths) > 2
            x = linspace(min(widths), max(widths));
            fittedFluxes = polyval(p,x);
            plot(x, fittedFluxes, '-');
            maxWidth = max(x);
            minWidth = min(x);
        end
        
        maxFlux = max(fluxes(validFluxes));
        plot([optimal.criticalWidth,optimal.criticalWidth],[0, maxFlux*1.2]);
        hold off;
        % Add extra blank line to title to avoid overlap with scientific
        % notation
        titleString = ['$\mathcal{C}=', num2str(CR), ', Ra=', num2str(Ra), '$'];
        if optimal.extrapolatedMax
            titleString = ['!', titleString, '!'];
        end
        
        if maxFlux < 0.012
            title({titleString, ''});
        else
            title(titleString);
        end
        
        
        
        xlabel('$L$');
        ylabh = ylabel('$F$');
        set(ylabh,'rotation',0)
        set(ylabh,'Units','normalized');
        set(ylabh,'position',get(ylabh,'position') - [0.1 0 0]);
        
        limits = [minWidth*0.9, maxWidth*1.1];
        if limits(2) <= limits(1) || isnan(limits(1)) || isnan(limits(2))
            if isnan(limits(1))
                limits(1) = 0;
            end
                
            limits(2) = limits(1) + 1;
        end
        
        xlim(limits);
        maxFluxMeasured = max(fluxes(validFluxes));
        minFluxMeasured = min(fluxes(validFluxes));
        
        fluxUpperPlotLimit = NaN;
        fluxLowerPlotLimit = NaN;
        if ~isnan(maxFluxMeasured) && maxFluxMeasured > 0.0
            fluxUpperPlotLimit = maxFluxMeasured*1.1;
        end
        
        if ~isnan(fluxLowerPlotLimit) && fluxLowerPlotLimit > 0.0
            fluxLowerPlotLtimit = fluxLowerPlotLimit*0.9;
        end
        
        if ~isnan(fluxLowerPlotLimit) && ~isnan(fluxUpperPlotLimit)
            ylim([fluxLowerPlotLimit fluxUpperPlotLimit]);
        end
        
        box on;
        
        % Increment this
        subplot_i = subplot_i + 1;
        
        % Calculate other optimal values by interpolating between values
        % near optimal width
        
        % Checks
        fieldsBetweenZeroOne = {'maxPorosityInterior', ...
            'avPorosityInterior'};
        
        % Lengths should be > 0
        fieldsGreaterThanZero = {'h','H'};
        
        for j = 1:numel(fields)
            % Don't fit flux and width, we already have these from the
            % quadratic
            thisField = fields{j};
            
            if strcmp(thisField, 'flux') || strcmp(thisField, 'fullWidth') || ...
                    strcmp(thisField, 'folder')
                continue
            end
            
            vals = [];
            for runs_i = 1:length(runs)
                try
                    vals(runs_i) = runs(runs_i).(thisField);
                catch E
                    vals(runs_i) = NaN;
                   err = 0; 
                end
            end
            
            
            validWidthsVals = (~isnan(widths)).*(~isnan(vals));
            
            if sum(validWidthsVals) >= 2
                
                % Linear fit seems safest for now
                try
                    optimal.(fields{j}) = interp1(widths(validWidthsVals==1), ...
                        vals(validWidthsVals==1), ...
                        optimal.fullWidth, 'PCHIP');
                catch E
                    optimal.(fields{j}) = NaN;
                end
                
            end
            
            if any(strcmp(fieldsBetweenZeroOne, fields{j}))
                if  optimal.(fields{j}) < 0 || optimal.(fields{j})  > 1
                    optimal.(fields{j}) = NaN;
                end
            end
            
            if any(strcmp(fieldsGreaterThanZero, fields{j}))
                if  optimal.(fields{j}) < 0
                    optimal.(fields{j}) = NaN;
                end
            end
            
            
        end
        
        
        
        % Find the folder nearest to the optimal cell width
        % i.e. optimal.folder = '' runs(run_i).cellsWide
        
        optimalCellWidth = optimal.cellsWide;
        
        thisCR = CR;
        thisRa = Ra;
        Le = 200;
        
        if isnan(optimalCellWidth)
            optimal.folder = '';
            
        else
            
            % Find folder nearest to this cell width
            roundedWidth = round(optimalCellWidth);
            testWidth = roundedWidth;
            
            searchResult = [];
            search_i = 1;
            
            CR_str = sprintf('%0.2f', thisCR+1);
            if strcmp(CR_str(end), '0')
                CR_str = CR_str(1:end-1);
            end
            
            %testName = ['CR',CR_str,'RaC',num2str(thisRa),'Le',num2str(Le),permeabilityType,'pts',num2str(testWidth),'-',folderEnding];
            %fprintf('%s \n', testName);
            
            a = runFolders(1);
            parts =strsplit(a{1}, 'pts');
            testSuffix = parts{1};
            
            testName = [testSuffix, 'pts', num2str(testWidth),'-',folderEnding];
            
            while 7~=exist([testName],'dir')
                
                %fprintf('Testing %s \n', testName);
                
                
                %testName = ['CR',CR_str,'RaC',num2str(thisRa),'Le',num2str(Le),permeabilityType,'pts',num2str(testWidth),'-',folderEnding];
                testName = [testSuffix, 'pts', num2str(testWidth),'-',folderEnding];
                
                %fprintf('%d\n', exist([data_dir, '/', testName],'dir'));
                
                % Generate new test width
                if mod(search_i, 2) == 0
                    testWidth = roundedWidth - ceil(search_i/2);
                else
                    testWidth = roundedWidth + ceil(search_i/2);
                end
                
                search_i = search_i + 1;
                
                
                if search_i < 0 || search_i > 500
                    testName = a{1};
                    break;
                end
                
            end
            
            %fprintf('%s \n', testName);
            optimal.folder = testName;
            
        end  
        
        
    else
        fprintf('  Length(runs) = %d \n', length(runs));
    end
    
    % Save final window
    % Save figure
    
    % Store optimal vals in MapN, indexed by (CR, Ra) value
    fprintf('Done, flux = %1.5f, full width = %1.3f \n', optimal.flux, optimal.fullWidth);
    optimalVals(CR, Ra)= optimal;
    
    
end


% Save data before making plots in case that fails
fprintf('Saving data');
save([data_dir, optimalValsFile], 'optimalVals')
save([data_dir, allStatesFile], 'allStates')

fid = fopen([data_dir, optimalFoldersFile], 'w');
for row = 1:length(optimalFolders)
    fprintf(fid, '%s\n', optimalFolders{row});
end
fclose(fid);



try
	%Output figures. Wrap this in a try/catch as it may fail on some systems

	set(h(window_i),'Units','Inches');
	pos = get(h(window_i),'Position');
	set(h(window_i),'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
	print(gcf, '-dpdf', [figure_output_dir, 'optStates',num2str(window_i),'.pdf'], printQuality);

	outputFiles = {};
	for i=2:window_i
	   outputFiles{end+1} = [figure_output_dir, 'optStates', num2str(i), '.pdf'];
	end
	delete([figure_output_dir, 'optStatesAll.pdf']);
	append_pdfs([figure_output_dir, 'optStatesAll.pdf'], outputFiles{:});

catch e %e is an MException struct
        fprintf('Could not save figures \n');
        fprintf(1,'The identifier was:\n%s \n',e.identifier);
        fprintf(1,'There was an error! The message was:\n%s \n',e.message);

        
    end




end
