% Load all optimal states
% Indexing is optimalStates(CR, Ra) = struct([runs], optimal
% state properties if appropriate)
% where runs = struct(dir, flux, width etc.)
clear all;
close all;

set(groot, 'DefaultTextInterpreter', 'latex')
set(groot, 'DefaultLegendInterpreter', 'latex')

data_dir = '/media/parkinsonjl/DATA/optimalStates-highRes/';
%data_dir = '/media/parkinsonjl/FREECOM HDD/optimalStates-lowRes/';
%data_dir = '/home/parkinsonjl/mnt/sharedStorage/AW002_PARKINSON_MUSH/optimalStates-highRes-new/';
%data_dir = '/home/parkinsonjl/mnt/sharedStorage/AW002_PARKINSON_MUSH/optimalStates-highRes-newDiags2/';

data_dir = '/home/parkinsonjl/mnt/sharedStorage/AW002_PARKINSON_MUSH/optimalStates-highRes-newSoluteFlux/';

forceRecalculateDiags = false;
forceRecalculateFlux = false;
allowCalcDiags = false; % don't do this for now - bugs in matlab code

printQuality = '-r300';

% Keep track of running folders which haven't reached steady state
noSteadyState = {};

if exist([data_dir, 'optimalStates.mat'], 'file') == 2
    load([data_dir, 'optimalStates.mat'], 'optimalStates')
else
    optimalStates = MapN();
end

if forceRecalculateFlux ||forceRecalculateDiags
    optimalStates = MapN();
end
% Uncomment to reload all data - this takes a long time! Be sure you want
% to do this!
%optimalStates = MapN();

blankState = struct('runs', [], 'optimalFlux', [], 'optimalWidth', [], ...
    'optimalHeight', [], 'optimalAvPerm', [], 'optimalAvPoros', [], ...
    'optimalMaxPsi', [], 'optimalConfH', []);

% Add files in directory below to the optimal states
folders = dir(data_dir);

% Check that all our runs still exist
% Also remove duplicates, if there are any
k = keys(optimalStates);
for i = 1:length(k)
    key = k(i); key = key{1};
    CR = key{1}; Ra = key{2};
    state = optimalStates(CR, Ra);
    
    % Check the state contains all the keys it should do
    blankKeys = fieldnames(blankState);
    for blank_i = 1:length(blankKeys)
        blankKey = blankKeys{blank_i};
        if ~isfield(state, blankKey)
            state.(blankKey) = [];
        end
    end
    
    run_i = 1;
    numRuns = length(state.runs);
    while run_i <= numRuns
        run = state.runs(run_i);
        full_dir = [data_dir, run.dir.name];
        %fprintf(['Checking that ', full_dir, ' still exists \n']);
        if exist(full_dir) ~= 7
            state.runs(run_i) = [];
            numRuns = numRuns - 1;
            fprintf(['Deleting runs with directory ', run.dir.name, '\n']);
        end
        
        % Check for duplicates
        for run_j = 1:length(state.runs)
            if run_i == run_j
                continue;
            end
            
            if run_i > numRuns
                break;
            end
            
            if state.runs(run_i).width == state.runs(run_j).width
                state.runs(run_j) = [];
                numRuns = numRuns - 1;
                break;
            end
            
        end
        
        run_i = run_i + 1;
    end
    
    if length(state.runs) == 0
        remove(optimalStates, CR, Ra);
    else
        optimalStates(CR, Ra) = state;
    end
    
end

% Temp hack:
%remove(optimalStates, 1.1, 350);
blankRuns = struct('dir', '', ...
    'flux', [], ...
    'width', [], ...
    'H', [], ...
    'averagePerm', [], ...
    'averagePoros', [], ...
    'maxPsi', [], ...
    'confinementH', []);

%folders_regex = 'CR(\d+\.?\d?)RaC(\d+)Le(\d+)(\s+)pts.*';
folders_regex = 'CR(\d+\.?\d+?)RaC(\d+\.?\d+?)Le(\d+)(\w+)pts(\d+)-0'; %Make sure it ends in 0
for i = 1:length(folders)
    folder = folders(i);
    
    
    
    newRun = blankRuns;
    newRun.dir = folder;
    
    [~,matches,~]  = regexp(folder.name, folders_regex, 'match', 'tokens', 'tokenExtents');
    if ~isempty(matches)
        match = matches{1};
        CR = str2double(match{1});
        Ra = str2double(match{2});
        Le = str2double(match{3});
        PermType = match{4};
        width = str2double(match{5});
        
        
        
        if isKey(optimalStates, CR, Ra)
            thisState = optimalStates(CR, Ra);
            
            %Check we don't already have this folder before overwriting it
            alreadyExists = false;
            for run_i = 1:length(thisState.runs)
                thisWidth = thisState.runs(run_i).width;
                if isequal(thisState.runs(run_i).dir.name,folder.name)
                    
                    alreadyExists = true;
                end
            end
            
            % check how recently pout.0 was edited
            % pout.0 is written to frequently i.e. every minute
            % so, if it hasn't been modified in a short period of 
            % time we must have finished the run
            steadyState = false;        
            poutInfo = dir([data_dir, folder.name, '/pout.0']);          
            now = datetime('now');
                 
            if abs(datenum(now) - poutInfo.datenum) > 0.001
                steadyState = true;
            else
                noSteadyState{end+1} = folder.name;
                fprintf([folder.name, ' has not reached steady state \n']);
            end
                

            if ~alreadyExists && ...
                    steadyState %exist(timetable, 'file') == 2
                thisState.runs(end+1) = newRun;
            end
            
            optimalStates(CR, Ra) = thisState;
        else
            
            newState =  blankState;
            newState.runs = [newRun];
            optimalStates(CR, Ra) = newState;
        end
    end
end


% Get all concentration ratios and rayleigh numbers
k = keys(optimalStates);
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


% Now, for each Ra, C pair, try and calculate optimal fluxes, wavelengths
% etc.

optimalFluxMat = NaN*ones(length(CR_arr), length(Ra_arr));
optimalWidthMat = optimalFluxMat;
optimalPorosMat = optimalFluxMat;
optimalPermMat = optimalFluxMat;
optimalPsiMat = optimalFluxMat;
optimalHeightMat = optimalFluxMat;
optimalConfHMat = optimalFluxMat;

for i = 1:length(k)
    key = k(i); key = key{1};
    CR = key{1}; Ra = key{2};
    state = optimalStates(CR, Ra);
    fprintf(['CR = ', num2str(CR), ', Ra = ', num2str(Ra), ', files = \n']);
    for run_i = 1:length(state.runs)
        %run = state.runs(run_i);
        run_dir = [data_dir, state.runs(run_i).dir.name];
        fprintf(['   ', state.runs(run_i).dir.name, '\n'])
        
        
        
        if ~isfield(state.runs(run_i), 'flux') || isempty(state.runs(run_i).flux) ...
                || forceRecalculateFlux
            % Let's calculate the fluxes and widths and store these
            
            
            inputs = readInputs([run_dir, '/inputs']);
            
            
            
            % Get (half) width for this run
        if isfield(inputs, 'domain_width')
             state.runs(run_i).width = str2double(inputs.domain_width);
        elseif isfield(inputs, 'domain_length')
            state.runs(run_i).width = str2double(inputs.domain_length);
        end
            
            % Get flux for this run
            
            
            diags = getDiagnostics([run_dir, '/diagnostics.out']);
             if isfield(diags, 'Fs_vertical_av')
                 state.runs(run_i).flux = diags.Fs_vertical_av(end);
             elseif isfield(diags, 'Fs_bottom')
                state.runs(run_i).flux = diags.Fs_bottom(end);
            else
                % Try pout file
                pout = Pout([run_dir, '/pout.0']);
                if isempty(pout.fluxBottom) || ~pout.steadyState
                    state.runs(run_i).flux = NaN;
                else
                    state.runs(run_i).flux = pout.fluxBottom(end);
                end
             end
             
             % Convert from calculated flux to flux where S=0 => F =0
             if  ~isnan(state.runs(run_i).flux)
                %state.runs(run_i).flux = state.runs(run_i).flux + CR;
                % Add one to convert to wells units, then multiply by -1 to
                % make positive
                
                % approximation to the flux in at the bottom of the domain
                % i.e. -1*V*(fraction of domain not covered by plume)
                %       = - (fraction of domain not covered by plume)
                inputFlux = -1.0; 
                state.runs(run_i).flux = inputFlux-(state.runs(run_i).flux);
             end 
            
             % A few checks
             if  state.runs(run_i).flux < 0
                  state.runs(run_i).flux  = NaN;
             elseif  state.runs(run_i).flux > 10
                  state.runs(run_i).flux  = NaN;
             end
            
            
        end
        
        % Might want to get other diagnostics here
        % average: permeability, streamfunction, porosity
        % mushy layer depth
        if (~isfield(state.runs(run_i), 'averagePerm') ...
                || isempty(state.runs(run_i).averagePerm) ...
                || forceRecalculateDiags) ...
                && allowCalcDiags
            
            pout = Pout([run_dir, '/pout.0']);
            inputs = readInputs([run_dir, '/inputs']);
            
            
            dim = 2; subcycled = true;
            
            % Skip if we don't appear to have any timesteps
            if length(pout.timesteps) == 0
                continue
            end
            frame = pout.timesteps(end);
            plot_prefix = inputs.plot_prefix;
            
            fprintf(['Final file = ', plot_prefix, num2str(frame), '\n']);
            
            output = MushyLayerOutput(dim, frame, [run_dir, '/'],...
                plot_prefix, subcycled);
            
            % Skip if this file doesn't exist
            if length(output.levelArray) == 0
                output = MushyLayerOutput(dim, frame+1, [run_dir, '/'],...
                    plot_prefix, subcycled);
                if length(output.levelArray) == 0
                    output = MushyLayerOutput(dim, frame+1, [run_dir, '/'],...
                        plot_prefix, subcycled);
                    if length(output.levelArray) == 0
                        continue
                    end
                end
            end
            
            %[chanWidth, chanDepth, numChannels, H] = output.channelGeometry(); % H is all I want here
            chanWidth = NaN; chanDepth = NaN; numChannels = NaN; H = NaN; % Some issues with output.channelGeometry() at the moment
            
            
            %[maxMushyU(file_i), maxMushyV(file_i), maxSpeed(file_i)] = output.maxMushyVel();
            % [maxLiqU(file_i), maxLiqV(file_i), ~] = output.maxLiqVel();
            %[channelDepthPorosity, averagePorosity,...
            %     halfPorosityDepth] = output.porosityMetrics(-0.5);
            %    [channelDepthPerm, averagePerm] = output.permMetrics();
            %   maxMushySaltAdv = output.maxMushyAdvectionSrcTerm();
            
            
            porosity = output.dataForComp(output.components.Porosity);
            psi = output.getStreamfunction(2000, 1).';
            permeability =  output.dataForComp(output.components.Permeability);
            w =  output.dataForComp(output.components.yAdvectionvelocity);
            
            dx = output.finest_dx();
            
            minPoros = 0.01;
            maxPoros = 0.9;
            
            permMush = permeability; psiMush = psi; porosMush = porosity;
            
            permMush(porosity < minPoros) = NaN; permMush(porosity > maxPoros) = NaN;
            psiMush(porosity < minPoros) = NaN; psiMush(porosity > maxPoros) = NaN;
            porosMush(porosity < minPoros) = NaN; porosMush(porosity > maxPoros) = NaN;
            
            psiMax = max(max(abs(psiMush)));
            
            % Get a new measure for the height, H
            % Need to find all points where a) w=0 and b) psi < psi_max / 5
            % Then find the lowest such point in the domain.
            psiCopy = psiMush;
            psiCopy = psi;
            wCopy = w;
            %tolerance = 0.05;
            %wCopy(w < -tolerance) = NaN;
            %wCopy(w > tolerance) = NaN;
            psiCopy(abs(psiCopy) > abs(psiMax)/5) = NaN;
            
            %psiCopy(w < -tolerance) = NaN;
            %psiCopy(w > tolerance) = NaN;
            [width, height] = size(psi);
            for j=1:height
                w_j = abs(w(:, j));
                smallestw = min(w_j);
                turningPoint = (w_j <= smallestw);
                for width_i = 1:width
                    if turningPoint(width_i) == 0
                        psiCopy(width_i, j) = NaN;
                    end
                end
            end
            
            % To test this, let's make a quick plot
            
            %             [X,Y] = output.domainGrid();
            %             figure();
            %
            %             pAxes = axes; %axis equal;
            %             colormap(pAxes,cool);
            %
            %             %h = pcolor(X,Y,psiCopy.');
            %             h = pcolor(X,Y,psiCopy.');
            %             set(h, 'EdgeColor', 'none');
            %             cAxes = axes; %axis equal;
            %             colormap(cAxes,[0 0 0]);
            %             contour(X,Y,(psi).', 50);
            %             axis(cAxes,'off')
            %             cAxes.Visible = 'off';
            %
            %             linkaxes([pAxes, cAxes]);
            %
            %             pause;
            
            
            % Get position of eutectic
            [eutecticx, eutecticy]= find(porosity < 0.001);
            eutectic_j = min(min(eutecticy));
            
            [confinementx, confinementy] = find(~isnan(psiCopy));
            confinement_j = min(min(confinementy));
            
            confinementHeight = (eutectic_j - confinement_j)*dx;
            
            
            state.runs(run_i).averagePerm =nanmean(nanmean(permMush)) ;
            state.runs(run_i).averagePoros = nanmean(nanmean(porosMush));
            state.runs(run_i).maxPsi = psiMax;
            state.runs(run_i).H = H;
            state.runs(run_i).confinementH = confinementHeight;
            
            % Replace data in optimalStates
            %state.runs(run_i) = run;
        end
        
        fprintf(['       flux = ', num2str(state.runs(run_i).flux), '\n'])
        fprintf(['       width = ', num2str(state.runs(run_i).width), '\n'])
    end
    
    % Now we've got all the fluxes for this (CR, Ra), try and calculate the
    % optimal flux and width
    
    runs = state.runs;
    
    % Need at least three points to fit a quadratic
    if length(runs) > 2
        
        p = fitFluxes(runs);
        
        % The max value of the polynomial is given by c - b^2/(4*a)
        state.optimalFlux = p(3) - (p(2)^2) / (4*p(1));
        
        % The position of the maximum is -b/(2a)
        state.optimalWidth = - p(2) / (2*p(1));
        
        % Calculate other optimal values
        state.optimalHeight = getOptimal(runs, state.optimalWidth, 'H');
        state.optimalAvPerm = getOptimal(runs, state.optimalWidth,'averagePerm');
        state.optimalAvPoros = getOptimal(runs, state.optimalWidth, 'averagePoros');
        state.optimalMaxPsi = getOptimal(runs, state.optimalWidth, 'maxPsi');
        state.optimalConfH = getOptimal(runs, state.optimalWidth, 'confinementH');
        
        
    end
    
    
    % Update any data we've changed
    optimalStates(CR, Ra)= state;
    
    % Fill out mega matrix with optimalFlux and optimalWidth so we can plot
    
    Ra_i = find(Ra_arr == Ra);
    CR_i = find(CR_arr == CR);
    if isfield(state, 'optimalFlux') && ~isempty(state.optimalFlux)
        % For old data:
%         if state.optimalFlux > 2
%             optimalFluxMat(CR_i, Ra_i) = state.optimalFlux - 2;
%         else
%             optimalFluxMat(CR_i, Ra_i) = state.optimalFlux - 1;
%         end
optimalFluxMat(CR_i, Ra_i) = state.optimalFlux;
        
    else
        optimalFluxMat(CR_i, Ra_i) = NaN;
    end
    
    % Checks
    if optimalFluxMat(CR_i, Ra_i) < 0.0 || optimalFluxMat(CR_i, Ra_i) > 10
        optimalFluxMat(CR_i, Ra_i) = NaN;
    end
    
    if isfield(state, 'optimalWidth') && ~isempty(state.optimalWidth)
        % For old data
        %optimalWidthMat(CR_i, Ra_i) = state.optimalWidth*4;
        optimalWidthMat(CR_i, Ra_i) = state.optimalWidth;
    else
        optimalWidthMat(CR_i, Ra_i) = NaN;
    end
    
    % Checks
    if optimalWidthMat(CR_i, Ra_i) < 0.0 || optimalWidthMat(CR_i, Ra_i) > 10.0
        optimalWidthMat(CR_i, Ra_i) = NaN;
    end
    
    if isfield(state, 'optimalAvPerm') && ~isempty(state.optimalAvPerm)
        optimalPermMat(CR_i, Ra_i) = state.optimalAvPerm;
    else
        optimalPermMat(CR_i, Ra_i) = NaN;
    end
    
    if isfield(state, 'optimalAvPoros') && ~isempty(state.optimalAvPoros)
        optimalPorosMat(CR_i, Ra_i) = state.optimalAvPoros;
    else
        optimalPorosMat(CR_i, Ra_i) = NaN;
    end
    
    if isfield(state, 'optimalMaxPsi') && ~isempty(state.optimalMaxPsi)
        if state.optimalMaxPsi > 0
            optimalPsiMat(CR_i, Ra_i) = state.optimalMaxPsi;
        else
            optimalPsiMat(CR_i, Ra_i) = NaN;
        end
    else
        optimalPsiMat(CR_i, Ra_i) = NaN;
    end
    
    if isfield(state, 'optimalHeight') && ~isempty(state.optimalHeight)
        if state.optimalHeight < 100 && state.optimalHeight > 0
            optimalHeightMat(CR_i, Ra_i) = state.optimalHeight;
        else
            optimalHeightMat(CR_i, Ra_i) = NaN;
        end
    else
        optimalHeightMat(CR_i, Ra_i) = NaN;
    end
    
    if isfield(state, 'optimalConfH') && ~isempty(state.optimalConfH)
        optimalConfHMat(CR_i, Ra_i) = state.optimalConfH;
    else
        optimalConfHMat(CR_i, Ra_i) = NaN;
    end
    
    fprintf('\n');
    
end

save([data_dir, 'optimalStates.mat'], 'optimalStates')


% Rescale things
Ra_crit = 0;
Ra_arr = (Ra_arr - Ra_crit);

CR_arr = CR_arr -1;

% Now we can plot flux vs width for a particular (CR, Ra)
m = 3; % Num rows
n = 3; % Num cols

subPlotsPerWindow = m*n;

window_i = 1;

while subPlotsPerWindow*(window_i-1) < length(k)
    
    h(window_i) = figure();
    set(h(window_i), 'Position', [100, 100, 1400, 800]);
    
    start_k = subPlotsPerWindow*(window_i-1) + 1;
    end_k = min(start_k + subPlotsPerWindow - 1, length(k));
    
    for i = start_k:end_k
        subplot_i = i-start_k + 1;
        key = k(i); key = key{1};
        CR = key{1}; Ra = key{2};
        state = optimalStates(CR, Ra);
        
        runs = state.runs;
        
        fluxes = NaN*ones(1, length(runs)); widths = NaN*ones(1, length(runs));
        
        for j = 1:length(runs)
            fluxes(j) = runs(j).flux;
            widths(j) = runs(j).width;
        end
        
        maxWidth = max(widths); minWidth = min(widths);
        
        subplot(m, n, subplot_i);
        hold on;
        plot(widths, fluxes, 'x');
        if length(widths) > 2
            p = fitFluxes(runs);
            fluxMax = -p(2)/(2*p(1));
            x = linspace(min(min(widths), fluxMax*0.9), ...
                max(max(widths), fluxMax*1.1));
            fittedFluxes = polyval(p,x);
            plot(x, fittedFluxes, '-');
            maxWidth = max(x);
            minWidth = min(x);
        end
        hold off;
        title(['$\mathcal{C}=', num2str(CR-1), ', Ra=', num2str(Ra), '$']);
        xlabel('$L$'); ylabel('$F_0$');
        xlim([minWidth*0.9, maxWidth*1.1]);
        maxFluxMeasured = max(fluxes);
        minFluxMeasured = min(fluxes);
        
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
        
        box on;
        
    end
    
    window_i = window_i + 1;
    
end


% Also want to plot how optimal state metrics vary

%CR_leg = {};
%for i = 1:length(CR_arr)
%    CR_leg{end+1} = ['CR', num2str(CR_arr(i))];
%end

hOptimal = figure();
set(hOptimal, 'Position', [100 100 1400 900]);
m = 2; % Num rows
n = 2; % Num cols

% In the first panel, plot F vs Ra for each CR
subplot(m, n, 1);
hold on;
CR_leg = {};
for i = 1:length(CR_arr)
    fluxForRa = optimalFluxMat(i, :);
    actualFluxes= ~isnan(fluxForRa);
    if sum(actualFluxes) > 0
        
        plot(Ra_arr(actualFluxes), fluxForRa(actualFluxes), 'x-');
        CR_leg{end+1} = ['$\mathcal{C}=', num2str(CR_arr(i)),'$'];
    end
end
hold off;
box on;
xlabel('$Ra$'); ylabel('$F_0$');
%legend(CR_leg);

% In the second panel, plot width vs Ra for each CR
subplot(m, n, 2);
hold on;

for i = 1:length(CR_arr)
    widthForRa = optimalWidthMat(i, :);
    actualWidths= ~isnan(widthForRa);
    if sum(actualWidths) > 0
        plot(Ra_arr(actualWidths), widthForRa(actualWidths), 'x-');
    end
    
end
hold off;
xlabel('$Ra$'); ylabel('$L$');
box on;
legend(CR_leg, 'Location', 'eastoutside');

% In the third panel, plot F vs CR for each Ra
subplot(m, n, 3);
hold on;
Ra_leg = {};
for i = 1:length(Ra_arr)
    fluxForRa = optimalFluxMat(:, i);
    actualFluxes= ~isnan(fluxForRa);
    if sum(actualFluxes) > 0
        
        plot(CR_arr(actualFluxes), fluxForRa(actualFluxes), 'x-');
        Ra_leg{end+1} = ['$Ra=', num2str(Ra_arr(i)),'$'];
    end
end
hold off;
box on;
xlabel('$C$'); ylabel('$F_0$');
%legend(Ra_leg);

% In the fourth panel, plot width vs CR for each Ra
subplot(m, n, 4);
hold on;
Ra_leg = {};
for i = 1:length(Ra_arr)
    widthForRa = optimalWidthMat(:, i) ;
    actualWidths= ~isnan(widthForRa);
    if sum(actualWidths) > 0
        
        plot(CR_arr(actualWidths), widthForRa(actualWidths), 'x-');
        Ra_leg{end+1} = ['$Ra=', num2str(Ra_arr(i)),'$'];
    end
end
hold off;
box on;
xlabel('$C$'); ylabel('$L$');
legend(Ra_leg, 'Location', 'eastoutside');


set(hOptimal,'Units','Inches');
pos = get(hOptimal,'Position');
set(hOptimal,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf, '-dpdf', ['optimalFluxes.pdf'], printQuality);


colors = {'m', 'c', 'r', 'g', 'b', 'k', 'y'};



%Figure just to analyse L scaling
hL = figure();
set(hL, 'Position', [200 200 1500 1000]);
m = 2; n = 2;

subplot(m, n, 1);
hold on;
for i = 1:length(Ra_arr)
    for j = 1:length(CR_arr)
        
        CR = CR_arr(j);
        Ra = Ra_arr(i);
        St = 5.0;
        
        Pi = optimalPermMat(j, i);
        chi = optimalPorosMat(j, i);
        psi = optimalPsiMat(j, i);
        flux = optimalFluxMat(j, i) - 1;
        L = optimalWidthMat(j, i);
        H = optimalHeightMat(j, i);
        
        maxChi = max(max(optimalPorosMat));
        minH = min(min(optimalHeightMat));
        
        % want to ensure this quantity is > 1. Scale it
        StHphi =St*H*(1-chi);
        StHphi = StHphi / (minH*(1-maxChi));
        
        plot(Ra, L, ['x', colors{1+mod(j, length(colors))}]);
        
    end
end
hold off;
box on;
xlabel('$Ra$'); ylabel('$L$');


subplot(m, n, 2);
hold on;
for i = 1:length(Ra_arr)
    for j = 1:length(CR_arr)
        
        CR = CR_arr(j);
        Ra = Ra_arr(i);
        St = 5.0;
        
        Pi = optimalPermMat(j, i);
        chi = optimalPorosMat(j, i);
        psi = optimalPsiMat(j, i);
        flux = optimalFluxMat(j, i);
        L = optimalWidthMat(j, i);
        H = optimalHeightMat(j, i);
        
        maxChi = max(max(optimalPorosMat));
        minH = min(min(optimalHeightMat));
        
        % want to ensure this quantity is > 1. Scale it
        StHphi =St*H*(1-chi);
        StHphi = StHphi / (minH*(1-maxChi));
        
        plot(CR, L, ['x', colors{1+mod(j, length(colors))}]);
        
    end
end
hold off;
box on;
xlabel('$C$'); ylabel('$L$');


subplot(m, n, 3);
hold on;
for i = 1:length(Ra_arr)
    for j = 1:length(CR_arr)
        
        CR = CR_arr(j);
        Ra = Ra_arr(i);
        St = 5.0;
        
        Pi = optimalPermMat(j, i);
        chi = optimalPorosMat(j, i);
        psi = optimalPsiMat(j, i);
        flux = optimalFluxMat(j, i);
        L = optimalWidthMat(j, i);
        H = optimalHeightMat(j, i);
        
        maxChi = max(max(optimalPorosMat));
        minH = min(min(optimalHeightMat));
        
        % want to ensure this quantity is > 1. Scale it
        StHphi =St*H*(1-chi);
        StHphi = StHphi / (minH*(1-maxChi));
        
        plot(H, L, ['x', colors{1+mod(j, length(colors))}]);
        
    end
end
hold off;
box on;
xlabel('$H$'); ylabel('$L$');



subplot(m, n, 4);
hold on;
CR_leg = {};
for j = 1:length(CR_arr)
    
    CR = CR_arr(j);
    Ra = Ra_arr;
    St = 5.0;
    
    Pi = optimalPermMat(j, :);
    chi = optimalPorosMat(j, :);
    psi = optimalPsiMat(j, i);
    flux = optimalFluxMat(j, :);
    L = optimalWidthMat(j, :);
    H = optimalHeightMat(j, :);
    
    maxChi = max(max(optimalPorosMat));
    minH = min(min(optimalHeightMat));
    
    % want to ensure this quantity is > 1. Scale it
    %StHphi =St*H*(1-chi);
    %StHphi = StHphi / (minH*(1-maxChi));
    
    porosScale = (CR./Ra).^(1/6);
    
    hasL = ~isnan(L);
    if sum(hasL) > 0
        CR_leg{end+1} = ['$\mathcal{C}=', num2str(CR),'$'];
        
        %plot(1./(St*(1-porosScale(hasL))), L(hasL), ['x-', colors{1+mod(j, length(colors))}]);
        plot(1./(St*(1-chi(hasL))), L(hasL), ['x-', colors{1+mod(j, length(colors))}]);
        
    end
    
end
legend(CR_leg, 'Location', 'eastoutside');
hold off;
box on;
xlabel('$1/(St \phi)$');
%xlabel('1/[St(1-(CR/Ra)^{1/6})]');
ylabel('$L$');


% Some scaling argument plots
% Want to plot:
% 1) psi vs L (delta chi) conc ratio
% 2) F vs St*CR*(1-chi)^2
hScaling = figure();
set(hScaling, 'Position', [50 50 1600 1000]);
m = 3; n = 3;


subplot(m, n, 1);
hold on;
for i = 1:length(Ra_arr)
    %for j = 1:length(CR_arr)
    
    CR = CR_arr.';
    St = 5.0;
    
    Ra = Ra_arr(i);
    
    chi = optimalPorosMat(:, i);
    psi = optimalPsiMat(:, i);
    flux = optimalFluxMat(:, i) ;
    L = optimalWidthMat(:, i);
    
    hasPsi = ~isnan(psi);
    
    chiScale = (CR/Ra).^(1/6);
    Lscale = 1./(St*(1-chiScale));
    
    if sum(hasPsi) > 0
        plot(CR(hasPsi).*(1-chi(hasPsi)).*L(hasPsi), psi(hasPsi), ['x-', colors{1+mod(i, length(colors))}]);
        %plot(CR(hasPsi).*(1-chiScale(hasPsi)).*Lscale(hasPsi), psi(hasPsi),...
        %    ['x-', colors{1+mod(i, length(colors))}]);
    end
    % end
end
hold off;
box on;
xlabel('$L  \phi  \mathcal{C}$'); ylabel('max$(\psi)$');

subplot(m, n, [2:3]);
hold on;
Ra_leg = {};
for i = 1:length(Ra_arr)
    %for j = 1:length(CR_arr)
    
    CR = CR_arr;
    St = 5.0;
    
    
    psi = optimalPsiMat(:, i);
    actualPsi = ~isnan(psi);
    
    if sum(actualPsi) > 1
        
        %chi = optimalPorosMat(j, i);
        
        %flux = optimalFluxMat(j, i);
        %L = optimalWidthMat(j, i);
        
        Ra_leg{end+1} = ['$Ra=', num2str(Ra_arr(i)),'$'];
        
        plot(CR(actualPsi)/St, psi(actualPsi), ['x-', colors{1+mod(i, length(colors))}]);
    end
    
    %end
end
legend(Ra_leg, 'Location', 'eastoutside');
hold off;
box on;
xlabel('$\mathcal{C}/\mathcal{S}$'); ylabel('max$(\psi)$');


subplot(m, n, 4);
hold on;

%for i = 1:length(Ra_arr)
for j = 1:length(CR_arr)
    
    CR = CR_arr(j);
    Ra = Ra_arr;
    St = 5.0;
    
    
    Pi = optimalPermMat(j, :);
    chi = optimalPorosMat(j, :);
    psi = optimalPsiMat(j, :);
    flux = optimalFluxMat(j, :) ;
    L = optimalWidthMat(j, :);
    
    
    % Ignore CR=0.75, Ra=100 value for now as it's dodgy
    if CR == 0.75
        Pi(Ra==100) = NaN;
    end
    
    hasPi = ~isnan(Pi);
    
    
    
    
    if sum(hasPi) > 0
        
        plot(sqrt(CR./Ra(hasPi)), Pi(hasPi), 'x-');
        
    end
    
end
%end
axis([0 0.1 0.05 0.25]);
hold off;
box on;
xlabel('$\sqrt{\mathcal{C}/Ra}$'); ylabel('$\Pi$');



subplot(m, n, [5:6]);
hold on;
CR_leg = {};
for j = 1:length(CR_arr)
    
    CR = CR_arr(j);
    Ra = Ra_arr;
    St = 5.0;
    
    Pi = optimalPermMat(j, :);
    chi = optimalPorosMat(j, :);
    psi = optimalPsiMat(j, i);
    flux = optimalFluxMat(j, :) ;
    L = optimalWidthMat(j, :);
    H = optimalHeightMat(j, :);
    
    maxChi = max(max(optimalPorosMat));
    minH = min(min(optimalHeightMat));
    
    % want to ensure this quantity is > 1. Scale it
    %StHphi =St*H*(1-chi);
    %StHphi = StHphi / (minH*(1-maxChi));
    
    Ra = Ra - 60;
    
    porosScale = (CR./Ra).^(1/6);
    
    hasL = ~isnan(L);
    if sum(hasL) > 0
        CR_leg{end+1} = ['$\mathcal{C}=', num2str(CR),'$'];
        plot(1./(St*(1-porosScale(hasL))), L(hasL), 'x-');
        
    end
    
end
legend(CR_leg, 'Location', 'eastoutside');
hold off;
box on;
%xlabel('1/\phi');
xlabel('$1/[\mathcal{S}(1-(\mathcal{C}/Ra)^{1/6})]$');
ylabel('$L$');


subplot(m, n, 7);
hold on;
%for i = 1:length(Ra_arr)
for j = 1:length(CR_arr)
    
    CR = CR_arr(j);
    Ra = Ra_arr;
    St = 5.0;
    chi = optimalPorosMat(j, :);
    flux = optimalFluxMat(j, :) ;
    L = optimalWidthMat(j, :);
    
    porosScale = (CR./Ra).^(1/6);
    
    Lscale = (1/St)./(1-porosScale);
    
    H = optimalHeightMat(j, :);
    
    
    %         hasFlux = ~isnan(flux);
    %         if sum(hasFlux) > 0
    %
    %         plot((chi(hasFlux) + CR*(1-chi(hasFlux)))./L(hasFlux), ...
    %             flux(hasFlux),  'x-');
    %
    %         end
    
    hasH = ~isnan(H);
    if sum(hasH) > 0
        
        plot(Ra(hasH), ...
            St*H(hasH),  'x-');
        
    end
    
end
%end
hold off;
box on;
%xlabel('St*(1-chi)*[chi + CR*(1-chi)]'); ylabel('F_0');
%xlabel('St*(CR + +(CR-1)*chi^2 + (1-2*CR)*chi]'); ylabel('F_0');
xlabel('$Ra$'); ylabel('$H \mathcal{S}$');

subplot(m, n, [8:9]);
hold on;
CR_leg = {};
%for i = 1:length(Ra_arr)
for j = 1:length(CR_arr)
    
    CR = CR_arr(j);
    Ra = Ra_arr;
    St = 5.0;
    averagePoros = optimalPorosMat(j, :);
    flux = optimalFluxMat(j, :) ;
    
    hasflux = ~isnan(flux);
    
    Ra_crit = 60-CR*2; % Smaller for larger CR
    Ra = Ra - Ra_crit;
    
    H = optimalHeightMat(j, :);
    
    
    porosScale = (CR./Ra).^(1/6);
    
    H = (H./min(H.*(1-porosScale)))/St;
    
    Lscale = H./sqrt(St*H.*(1-porosScale) -1);
    Lscale =  1./(St*(1-porosScale.^2));
    
    
    if sum(hasflux) > 0
        
        CR_leg{end+1} = ['$\mathcal{C}=', num2str(CR),'$'];
        
        %plot((1-porosScale(hasflux)).^2*CR*St, flux(hasflux), ['x-', colors{1+mod(j, length(colors))}]);
        
        plot(CR*St*((1-porosScale(hasflux)).^2), flux(hasflux), ['x-']);
        %plot(CR./Lscale(hasflux), flux(hasflux), ['x-']);
        
    end
end
%end
legend(CR_leg, 'Location', 'eastoutside')
hold off;
box on;

%xlabel('$\mathcal{S} \mathcal{C} (1-(\mathcal{C}/[Ra-Ra_{crit}])^{1/6} )^2$'); ylabel('$F_0$');
%xlabel('$\mathcal{C} (1-\chi_{scale})/L_{scale}$'); ylabel('$F_0$');
xlabel('$\mathcal{S}*\mathcal{C}*(1-(\mathcal{C}/Ra^*)^{1/6} )^2$'); ylabel('$F_0$');


set(hScaling,'Units','Inches');
pos = get(hScaling,'Position');
set(hScaling,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf, '-dpdf', ['scalingPlots.pdf'], printQuality);




% Some more scaling plots

hScaling2 = figure();
set(hScaling2, 'Position', [50 50 1600 1000]);
m = 2; n = 2;


subplot(m, n, 1);
hold on;

CR_leg = {};
%for i = 1:length(Ra_arr)
for j = 1:length(CR_arr)
    
    CR = CR_arr(j);
    Ra = Ra_arr;
    St = 5.0;
    averagePoros = optimalPorosMat(j, :);
    flux = optimalFluxMat(j, :) ;
    
    hasflux = ~isnan(flux);
    
    Ra_crit = 20 + 3/CR; % Smaller for larger CR
    %Ra_crit = 0;
    Ra = Ra - Ra_crit;
    
    H = optimalHeightMat(j, :);
    
    
    porosScale = (CR./Ra).^(1/6);
    
    H = (H./min(H.*(1-porosScale)))/St;
    
    Lscale = H./sqrt(St*H.*(1-porosScale) -1);
    Lscale =  1./(St*(1-porosScale.^2));
    
    
    if sum(hasflux) > 0
        
        CR_leg{end+1} = ['$\mathcal{C}=', num2str(CR),'$'];
        
        %plot((1-porosScale(hasflux)).^2*CR*St, flux(hasflux), ['x-', colors{1+mod(j, length(colors))}]);
        
        plot(Ra(hasflux).^(1/6), flux(hasflux)/CR, ['x-']);
        
        %plot(1 - 2*(CR./Ra(hasflux).^(1/6)), flux(hasflux)/CR, ['x-']);
        
    end
end
%end
legend(CR_leg, 'Location', 'eastoutside')


hold off;
box on;
ylabel('$F/\mathcal{C}$');
xlabel('$(Ra-Ra_{crit})^{1/6}$'); 
%xlabel('$1 - 2 (\mathcal{C} / (Ra-Ra_{crit}))^{1/6}$'); 



% Provide some info on currently runnning simulations
for i = 1:length(noSteadyState)
   fprintf('%s has not reached steady state \n', noSteadyState{i}); 
end

