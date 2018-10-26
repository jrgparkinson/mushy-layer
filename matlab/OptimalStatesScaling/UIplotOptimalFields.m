% GUI to display the spatial structure of fields for
% the chosen (CR, Ra) pairs

function UIplotOptimalFields

%clear all;
%close all;

katzUnits = true;
data_dir = ['/home/parkinsonjl/mnt/sharedStorage/', ...
    'optimalStates-highRes-new/'];

% katzUnits = false;
% data_dir = ['/home/parkinsonjl/mnt/sharedStorage/', ...
%    'optimalStates-highRes-Wells/'];

allStatesFile = 'allStates.mat';
optimalValsFile = 'optimalVals.mat';

% This is fixed
Le = 200;

% Possibly make these changeable options
minPorosity = 0.05;
maxPorosity = 0.97;

% First get the available optimal states
allStates = [];
optimalVals = [];

fprintf('Loading optimal states... ');

loadData([data_dir, allStatesFile]);
loadData([data_dir, optimalValsFile]);

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

% Finished loading data
fprintf(' done \n');

% Choose some initial Ra/CR
initial_Rai = round(length(Ra_arr)/2);
initial_CRi = round(length(CR_arr)/2);

Ra = Ra_arr(initial_Rai);
CR = CR_arr(initial_CRi);


% Make GUI
f = figure();
set(f, 'Position', [100 200 1600 1000]);
%ax = axes('Parent', f, 'position', [0.35 0.15 0.6 0.7]);

popupRaPos = [70 700 100 50];
popupRa = uicontrol('Style', 'popup',...
    'String', num2cell(Ra_arr),...
    'Position', popupRaPos,...
    'Callback', @setRa);
popupRa.Value = initial_Rai;

txtRa = uicontrol('Style','text',...
    'Position',popupRaPos - [50 0 50 0],...
    'String','Ra:'); %$\mathcal{R}$

popupCRPos =  popupRaPos - [0 30 0 0];
popupCR = uicontrol('Style', 'popup',...
    'String', num2cell(CR_arr),...
    'Position',popupCRPos,...
    'Callback', @setCR);
popupCR.Value = initial_CRi;

txtCR = uicontrol('Style','text',...
    'Position',popupCRPos - [50 0 50 0],...
    'String','CR:'); %$\mathcal{C}$:

txtMaxPorosity = uicontrol('Style','text',...
    'Position',txtCR.Position + [-10 -50 80 0],...
    'String',['Max porosity (',num2str(maxPorosity),'):']); %$\mathcal{C}$:

sldMaxPorosity = uicontrol('Style', 'slider',...
    'Min',0,'Max',0.97,'Value',maxPorosity,...
    'Position', txtMaxPorosity.Position + [0 -30 -20 0],...
    'Callback', @changeMaxPorosity);

txtInfo = uicontrol('Style','text',...
    'Position',sldMaxPorosity.Position + [0 -70 0 0],...
    'String',''); %$\mathcal{C}$:



m=3; n=4; % For subplots

% Make these function global
plotFile = [];
inputs = [];
fields = struct();
T = []; chi = [];

updatePlot(true);


    function changeMaxPorosity(source, event)
        maxPorosity = source.Value;
        fprintf('Max Porosity: %d \n', maxPorosity);
        txtMaxPorosity.String = ['Max porosity (',num2str(maxPorosity),'):'];
        drawnow;
        
        updatePlot();
    end

    function setCR(source, event)
        CR = str2num(source.String{source.Value});
        fprintf('CR: %d \n', CR);
        
        updatePlot(true);
    end

    function setRa(source, event)
        Ra = str2num(source.String{source.Value});
        fprintf('Ra: %d \n', Ra);
        
        updatePlot(true);
    end

    function updatePlot(newCRRa)
        
        if nargin < 1
            newCRRa = false;
        end
        
        if ~isKey(optimalVals, CR, Ra)
            txtInfo.String = 'No data for parameters selected';
            return;
        else
            optimal = optimalVals(CR, Ra);
            if strcmp(optimal.folder, '') || strcmp(optimal.folder, 'NaN')
                txtInfo.String = 'No data for parameters selected';
                return;
            end
        end
        
        plotFolder = [data_dir, optimal.folder];
        
        
        txtInfo.String = 'Loading data... ';
        drawnow;
        
        fprintf('Loading data for (CR, Ra) pair... \n ');
        fprintf('Folder: "%s" \n', optimal.folder);
        
        tic; 
        
        if newCRRa
            plotFile = getFinalPlotFile(plotFolder);
            
            if length(plotFile.levelArray()) == 0
                fprintf(' could not load plot file\n');
                txtInfo.String = 'No data available for these parameters';
                return
            end
            
            %txtInfo.String = 'Getting data for new plot...';
            
            inputs = readInputs([plotFolder, '/inputs']);
            %diagFile = getDiagnostics([data_dir, folderName, '/diagnostics.out']);
            
            
            
            T = plotFile.dataForComp(plotFile.components.Temperature).';
            chi = plotFile.dataForComp(plotFile.components.Porosity).';
            
            
            
            
            % Convert to Wells Units
            if katzUnits
            T = T - 1;
            end
        end
        
        fprintf('Finished getting data, toc = %1.3f \n', toc);
        txtInfo.String = 'Processing data... ';
        drawnow;
        
        St = str2double(inputs.stefan);
        frameAdv = str2double(inputs.nonDimVel);
        dx = plotFile.finest_dx();
        
        [fields.heatAdvection, fields.heatDiffusion, fields.latentHeat, fields.TFrameAdvection, ...
            fields.saltAdvection, fields.saltDiffusion, fields.liquidSalinityFrame, fields.solidSalinityFrame, ...
            fields.vorticityDiffusion, fields.baroclinicTorque, fields.vorticityPermeability] = ...
            computeFields(plotFile, frameAdv, St, CR, Ra, Le);
        fields.chi = chi;
        
         fprintf('Finished computing fields data, toc = %1.3f \n', toc);
        % Apply porosity filter
        
        mushyLayer = (chi > minPorosity).*(chi < maxPorosity).*(T > -0.95);
        
       
        
        
        fieldsList = fieldnames(fields);
        for f_i = 1:length(fieldsList)
            name= fieldsList{f_i};
            
            fields.(name)(mushyLayer~=1) = NaN;
            
            % porosityMush = chi; porosityMush(mushyLayer ~= 1) = NaN;
            
        end
        
        fprintf('Finished applying porosity filter, toc = %1.3f \n', toc);
       
        
        titles = fieldsList;
        titles{1} = '$\mathbf{U}  \cdot \nabla \theta$';
        titles{2} = '$\nabla^2 \theta$';
        titles{3} = '$\mathcal{S} d\chi/dz$';
        titles{4} = '$d\theta/dz$';
        titles{5} = '$-\mathbf{U}  \cdot \nabla \theta$';
        titles{6} = '$\nabla \cdot \chi \nabla \theta$';
        titles{7} = '$d(\theta \chi)/dz$';
        titles{8} = '$\mathcal{C} d \chi/dz$';
        titles{9} = '$\nabla^2 \psi$';
        titles{10} = '$Ra \Pi d\theta/dx$';
        titles{11} = '$\nabla \Pi \cdot \nabla \psi / \Pi$';
        titles{12} = '$\chi$';
        
         % Get axis limits from mushyLayer region
        [YI,XI] = find(mushyLayer==1);
        minY = min(YI)*dx; maxY = max(YI)*dx;
        minX = min(XI)*dx; maxX = max(XI)*dx;
        
        x = 0:dx:(maxX-minX);
        y = 0:dx:(maxY-minY);
        
        maxX = maxX - minX;
        maxY = maxY - minY;
        
        minX = 0; minY = 0;
        
        
        
        [X,Y] = meshgrid(x,y);
        
        for f_i = 1:length(fieldsList)
            fieldname = fieldsList{f_i};
            thisField = fields.(fieldname);
            
            
            largestVal = max(max(abs(thisField)));
            smallestVal = min(min(thisField));
            
            subplot(m, n, f_i);
            
            
            
            plotField = thisField(min(YI):max(YI), min(XI):max(XI));
            if strcmp(fieldname, 'chi')
                contour(X, Y, plotField);
                ax = gca; ax.CLim = [smallestVal largestVal];
                colormap(gca, 'jet');
                
            else
                h = pcolor(X, Y, plotField);
                %set(h, 'EdgeColor', 'none');
                ax = gca; ax.CLim = [-largestVal largestVal]; % Symmetric colorbars
                colormap(gca, bluewhitered);
            end
            
            
            
            colorbar(); set(gca,'YDir','normal')
            title(titles{f_i}); % blank first line to move plots down slightly
            box on;
            set(ax,'Layer','top')
            axis image;
           % xlim([minX, maxX]);
            ylim([minY, maxY]);
            
           
            xMaxStr = sprintf('%1.2f', maxX);
            ax.XTick =  [minX maxX];
            ax.XTickLabel = {'0', xMaxStr};
            
            yMaxStr = sprintf('%1.2f', maxY);
            ax.YTick =  [maxY];
            ax.YTickLabel = {yMaxStr}; % No need to label origin as already done in x direction
            
            
            %set(gca, 'OuterPosition', [0, 0.76, 0.49, 0.23])
            
        end
        
        
        
        txtInfo.String = ['Displaying data for Ra=',num2str(Ra),', CR=',num2str(CR),'.'];
        
        
    end


    function clearPlots()
        cla reset;
    end

    function loadData(file)
        load(file);
    end

end