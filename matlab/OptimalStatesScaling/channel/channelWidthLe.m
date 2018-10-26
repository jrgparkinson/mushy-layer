% Analyse observed channel width as a function of Lewis number
close all;

optimalRa = true;

data_dir = getDataDir('saltDiffusionChannelWidth/');

if optimalRa
    Ra_data_dir = getDataDir('optimalStates-highRes-new/');
    load([Ra_data_dir, 'optimalVals.mat']);
    
else
    Ra_data_dir = data_dir; % g
end

figure_output_dir = ['/home/parkinsonjl/convection-in-sea-ice/MushyLayer/', ...
    'matlab/Scaling/channel/'];

getData = true;

St = 5;

CRForLeVariation = 1.0;
RaForLeVariation = 60;

CRForRaVariation = 0.5;
LeForRaVariation = 200;

if getData
    
    % Load all folders
    folders  = dir(data_dir);
    
    Le_arr = NaN*ones(length(folders), 1);
    
    a_arr = Le_arr;
    h_arr = Le_arr;
    F_arr = Le_arr;
    a_uncertainty = Le_arr;
    
    streamfunction = Le_arr;
    
    unsteadyRuns = 0*ones(length(folders), 1);
    
    folders_regex = 'CR(\d+\.?\d+?)RaC(\d+\.?\d+?)Le(\d+\.?\d+?)(\w+)pts(\d+)-(\w+)';
    for i = 1:length(folders)
        folder = folders(i);
        
        [~,matches,~]  = regexp(folder.name, folders_regex, 'match', 'tokens', 'tokenExtents');
        if ~isempty(matches)
            match = matches{1};
            CR = str2double(match{1}) - 1; % This is where we transform conc ratio vals
            Ra = str2double(match{2});
            Le = str2double(match{3});
            PermType = match{4};
            width = str2double(match{5});
            ending = match{6};
            
            % These states appear to have no flow or we don't want to
            % consider them
            if Le > 1500 || Ra ~= RaForLeVariation || CR ~= CRForLeVariation
                continue
            end
            
            
            
            fprintf('Loaded Le = %d \n', Le);
            
            
            if strcmp(ending, 'steady')
                unsteadyRuns(i) = false;
            else
                unsteadyRuns(i) = true;
            end
            
            
            Le_arr(i) = Le;
            
            % Get parameters
            filename = [data_dir, '/', folder.name, '/postProcess.csv'];
            
            if ~exist(filename, 'file')
                fprintf('  No postProcess.csv \n');
                continue
            end
            
            A = importdata(filename, ',', 1);
            
            
            postProcessStruct = cell2struct(num2cell(A.data), A.colheaders, 2);
            
            inputs = readInputs( [data_dir, '/', folder.name, '/inputs']);
            
            if isfield(inputs, 'domain_width')
                dom_width = str2num(inputs.domain_width);
            else
                dom_width = str2num(inputs.domain_length);
            end
            
            
            num_cells = inputs.num_cells;
            parts= strsplit(num_cells, ' ');
            Nx = parts{1};
            Nx = str2num(Nx);
            
            dx = dom_width/Nx;
            
            a_uncertainty(i) = dx/2;
            
            
            streamfunction(i) = postProcessStruct.avStreamfunctionMush;
            
            height = NaN;
            if isfield(postProcessStruct, 'channelHeight')
                height = postProcessStruct.channelHeight;
                height = postProcessStruct.h;
            end
            
            h_arr(i) = height;
            
            width = NaN;
            if isfield(postProcessStruct, 'channelWidth')
                if postProcessStruct.channelWidth < 0.1
                    width = postProcessStruct.channelWidth;
                else
                    fprintf('Channel width of %1.5f is unrealistic \n', postProcessStruct.channelWidth);
                end
            end
            
            % empirically observe that channel is one grid cell for these
            % lewis numbers
            if Le >= 5000
                width = dx;
            end
            
            a_arr(i) = width;
            
            flux = NaN;
            if isfield(postProcessStruct, 'flux')
                flux = postProcessStruct.flux;
            end
            F_arr(i) = -flux;
            
            
            
            % for j=1:length(dSldz)
            %  [T, Sl, psi, a_pred(j, i), Sa, x] = analyticSolChannel(Le_arr(i),...
            %      Ra, dSldz(j), T0, S0);
            % end
            
            
            
            
        end
        
    end
    
    
    
    % Also get data for a ~ Ra
    a_Ra_arr  = []; H_Ra_arr = []; channelHeight_Ra_arr = []; confinement_Ra_arr = [];
    Ra_arr = [];
    folders = dir(Ra_data_dir);
    
    for f_i = 1:length(folders)
        folder = folders(f_i);
        [~,matches,~]  = regexp(folder.name, folders_regex, 'match', 'tokens', 'tokenExtents');
        if ~isempty(matches)
            match = matches{1};
            CR = round(str2double(match{1}) - 1, 3); % This is where we transform conc ratio vals
            Ra = str2double(match{2});
            Le = str2double(match{3});
            PermType = match{4};
            width = str2double(match{5});
            ending = match{6};
            
            if round(CR,3) ~= CRForRaVariation
               % fprintf('%1.2f - %1.2f = %1.10f \n', CR, CRForRaVariation, CR-CRForRaVariation);
                continue
            end
            
            if round(Le) ~= LeForRaVariation
                continue
            end
            
            % Only do optimal folders:
            if optimalRa
                if isKey(optimalVals,CR,Ra)
                    
                    optVals = optimalVals(CR,Ra);
                    optFolder = optVals.folder;
                    if ~strcmp(optFolder, folder.name)
                        continue
                    end
                else
                    continue
                end
                
            else
                
                % Only do for fixed width:
                if width ~= 128
                    %if width ~= 96
                    %fprintf('Ignore Ra=%d, width=%d \n', Ra, width);
                    continue
                end
                
                if find(Ra_arr == Ra) > 0
                    continue
                end
                
            end
            
            fprintf('Loaded Ra = %d \n', Ra);
            
            % Get parameters
            filename = [Ra_data_dir, '/', folder.name, '/postProcess.csv'];
            
            if ~exist(filename, 'file')
                fprintf('  No postProcess.csv \n');
                continue
            end
            
            A = importdata(filename, ',', 1);
            postProcessStruct = cell2struct(num2cell(A.data), A.colheaders, 2);
            
            inputs = readInputs( [Ra_data_dir, '/', folder.name, '/inputs']);
            
            
            
            Ra_arr(end+1) = Ra;
            
            width = NaN;
            if ~isfield(postProcessStruct, 'channelWidth')
                makePlots = false; redoRuns = true;
                postProcessFolder(folder.name, makePlots, redoRuns);
                
                % Get post processed data again
                A = importdata(filename, ',', 1);
                postProcessStruct = cell2struct(num2cell(A.data), A.colheaders, 2);
            end
            
            if isfield(postProcessStruct, 'channelWidth')
                if postProcessStruct.channelWidth < 0.15
                    width = postProcessStruct.channelWidth;
                else
                    fprintf('Channel width of %1.5f is unrealistic \n', postProcessStruct.channelWidth);
                end
            else
                fprintf('No channelWidth field \n');
                
            end
            
            %height = postProcessStruct.channelHeight;
            
            H_Ra_arr(end+1) = postProcessStruct.H;
            channelHeight_Ra_arr(end+1) = postProcessStruct.channelHeight;
            confinement_Ra_arr(end+1) = postProcessStruct.h;
            
            a_Ra_arr(end+1) = width;
        end
        
        
    end
    
end % end if get data


h = figure();
set(h, 'Position', [100 400 1600 500]);
n = 4;
m = 1;
subplot(m, n, 1);

%yyaxis left;
%plot(log1010(Le_arr), a_arr, 'x');
% ylabel('a');

aPowerLe = 2;

x = log10(Le_arr*St);

%y = log10(h_arr./(a_arr.^aPowerLe));%
y = log10(1./(a_arr.^aPowerLe)); %h_arr

notnan = ~isnan(x).*~isnan(y);
% Also ignore very large Le runs

%toFit = notnan.*(Le_arr < 50000).*(~unsteadyRuns);
toFit = notnan.*(Le_arr < 50000);
xFit = x(toFit==1); yFit = y(toFit==1);

xPlot = x(notnan==1); %y=y(notnan==1);

%p = polyfit(xFit,yFit, 1);
%p2 = p; p2(1) = 0.5; p2(2) = p2(2) + 0.08;
%log10_a_fit = polyval(p, xPlot);
%log10_a_fit2 = polyval(p2, xPlot);

alog10_uncertainty = a_uncertainty./a_arr;

xPred = sort(xFit);
xPred = [xPred(1) xPred(end)];
predMin = min(yFit)*0.8;
predMax =predMin + abs(xPred(2)-xPred(1));
%predMin = predMax - 0.5*abs(xPred(2)-xPred(1));%predMax*((xPred(1)/xPred(2))^(0.5));
log10_predictedScaling = [predMin predMax];
pred2 = [predMin*1.15 predMin*1.15 + 0.75*abs(xPred(2)-xPred(1))];

hold on
%data = errorbar(x, y, alog10_uncertainty, 'x', 'LineWidth', 2.0);
data = plot(x, y, 'x', 'LineWidth', 2.0);
%fit = plot(xPlot, log10_a_fit, '-');

predictedScaling = plot(xPred, log10_predictedScaling, ':');
predictedScaling = plot(xPred, pred2, ':');

legend({'Data', '$\sim Le$', '$\sim Le^{3/4}$'}, 'Location', 'northwest');

daspect([1 1 1]);

thisAx = gca;
axPos = thisAx.Position;
ylabh = ylabel(['log$_{10}(1/a^',num2str(aPowerLe),')$']); %'Rotation', 0, 'units', 'normalized', ...
    %'Position', [axPos(1)-0.3 0.45 0]);
xlabel('log$_{10}(Le \, St)$');

box on;
title(['$\mathcal{C} = ',num2str(CRForLeVariation),', Rm_S = ',num2str(RaForLeVariation),'$']);



subplot(m, n, 2);

Ra_star = (Ra_arr - 30)/30;

x = log10(Ra_arr);
%x = log10(1./Ra_star);
y = log10(a_Ra_arr); %h_Ra_arr
a_SqrtH = log10(a_Ra_arr ./sqrt(H_Ra_arr));

a_SqrtChannelH = log10(sqrt(channelHeight_Ra_arr) ./ a_Ra_arr);

%y = a_SqrtChannelH;
y = log10(1./a_Ra_arr.^2);
%y = log10(channelHeight_Ra_arr./a_Ra_arr.^2);


xValid = sort(x(~isnan(x)));
xPred = [xValid(1) xValid(end)];
aPred = [min(y)*1 min(y)*1+2*(xPred(2)-xPred(1))];
aPred3 = [min(y)*1 min(y)*1+1.5*(xPred(2)-xPred(1))];
aPred2 = [min(y)*1 min(y)*1+(xPred(2)-xPred(1))];

error = max(a_uncertainty)./a_Ra_arr;

hold on;
%errorbar(x, y, error, 'x', 'LineWidth', 2.0);
plot(x, y, 'x', 'LineWidth', 2.0);
%plot(log10(1./Ra_arr), log10(a_Ra_arr), 'x');
%plot(x, a_SqrtChannelH-0.5, 'x');
%plot(x, a_SqrtH-0.5, 'x');


plot(xPred, aPred, ':');
plot(xPred, aPred3, ':');
plot(xPred, aPred2, ':');

%x = linspace(1.5,2.5,20); y=2+0.5*(x-2);
%plot(x,y, '-')
%plot(x, x, '-');


hold off;

xlabel('log$_{10}(Rm_S)$');

thisAx = gca;
axPos = thisAx.Position;
ylabel('log$_{10}(1/a^2)$'); % 'Rotation', 0, 'units', 'normalized', ...
    %'Position', [axPos(1)-0.7 0.45 0]);

legend('Data', ...
    '$\sim Rm_S^{2}$', ...
     '$\sim Rm_S^{3/2}$', ...
    '$\sim Rm_S$', ...
    'Location', 'northwest');

box on;

title(['$\mathcal{C} = ',num2str(CRForRaVariation),', Le = ',num2str(LeForRaVariation),'$']);

%set(h,'Units','Inches');
%pos = get(h,'Position');
%set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%print(gcf, '-dpdf', [figure_output_dir, 'ChannelWidthLeDarcy.pdf'], '-r50');



subplot(m, n, 3);
% plot(log10(Ra_arr), log10(channelHeight_Ra_arr));
%  ylabel('$log10(h_{channel})$');
%  
 plot(log10(Ra_arr), log10(1./confinement_Ra_arr));
 ylabel('$log10(1/h)$');

 
xlabel('$log10(Ra)$');

subplot(m, n, 4);
plot(Le_arr, h_arr, 'x');
xlabel('$Le$'); ylabel('$h_{channel}$');

