% Analyse observed channel width as a function of Lewis number
close all;

set(groot, 'defaultAxesFontName', 'times');

getData = true;
if getData
     clear all;
     getData = true;
end

variableHscaling = false;
St = 5;

data_dir = getDataDir('saltDiffusionChannelWidth/');

%figure_output_dir = ['/home/parkinsonjl/convection-in-sea-ice/MushyLayer/', ...
%    'matlab/Scaling/channel/'];
figure_output_dir = data_dir;


set(groot, 'DefaultAxesLineStyleOrder', 'x-|x--|x:');

%Ra = 50;
%dSldz = []; %[3 4 5];
%T0 =0.5;
%S0 = T0+0.01;

RA_MAX = 20000;
LE_MAX = 5000;
RA_MIN = 700;

if getData
    
    % Load all folders
    folders  = dir(data_dir);
    
  
    a = MapN();
    h = MapN();
    unsteadyRuns = MapN();
    
    folders_regex = 'CR6.0RaC(\d+\.?\d+?)Le(\d+\.?\d+?)(\w+)pts(\d+)-(\w+)';
    for i = 1:length(folders)
        folder = folders(i);
        
        [~,matches,~]  = regexp(folder.name, folders_regex, 'match', 'tokens', 'tokenExtents');
        if ~isempty(matches)
            match = matches{1};
            %CR = str2double(match{1}) - 1; % This is where we transform conc ratio vals
            Ra = str2double(match{1});
            Le = str2double(match{2});
            PermType = match{3};
            Nx = str2double(match{4});
            ending = match{5};
            
            if Ra > RA_MAX || Le > LE_MAX || Ra < RA_MIN
                continue
            end
%             if CR ~= 6.0
%                 continue
%             end
            
            
            
            fprintf('Loaded Le = %d, Ra = %d, \n', Le, Ra);
            
            
            if strcmp(ending, 'steady')
                unsteadyRuns(Ra,Le) = false;
            else
                unsteadyRuns(Ra,Le) = true;
            end
            
            
           
            
            % Get parameters
            filename = [data_dir, '/', folder.name, '/postProcess.csv'];
            
            if ~exist(filename, 'file')
                fprintf('  No postProcess.csv \n');
                continue
            end
            
            A = importdata(filename, ',', 1);
            
            
            postProcessStruct = cell2struct(num2cell(A.data), A.colheaders, 2);
            
            inputs = readInputs( [data_dir, '/', folder.name, '/inputs']);
            

            chanWidth = NaN;
            if isfield(postProcessStruct, 'channelWidth')
                if postProcessStruct.channelWidth < 0.2
                    chanWidth = postProcessStruct.channelWidth;
                else
                    fprintf('Channel width of %1.5f is unrealistic \n', postProcessStruct.channelWidth);
                end
            end
            
            % If height is much less than width then this isn't a channel
            if postProcessStruct.channelHeight < 0.1*chanWidth
                chanWidth = NaN;
            end
            
            a(Ra,Le) = chanWidth;
           % h(Ra,Le) = postProcessStruct.h;
            h(Ra,Le) = postProcessStruct.channelHeight;
            %h(Ra,Le) = postProcessStruct.H;
            
            
        end
        
    end
    
    
end % end if get data

k = keys(a);
Le_arr = [];
Ra_arr = [];
for i = 1:length(k)
    keyPair = k(i);
    keyPair = keyPair{1};
    temp = keyPair(2);  thisLe  = temp{1};
    temp = keyPair(1);  thisRa = temp{1};
    
    
    if  ~ismember(thisLe, Le_arr)
        Le_arr(end+1) = thisLe;
    end
    
    if ~ismember(thisRa, Ra_arr)
        Ra_arr(end+1) = thisRa;
    end
    
end

Ra_arr = sort(Ra_arr); Le_arr = sort(Le_arr);

% Create matrix of a data
a_mat = NaN*ones(length(Ra_arr), length(Le_arr));
h_mat = NaN*ones(length(Ra_arr), length(Le_arr));
for Le_i=1:length(Le_arr)
    Le = Le_arr(Le_i);    
    for Ra_i=1:length(Ra_arr)
        Ra = Ra_arr(Ra_i);
        
        if isKey(a,Ra,Le)
            a_mat(Ra_i, Le_i) = a(Ra,Le);
            h_mat(Ra_i, Le_i) = h(Ra,Le);
        end
        
    end
end


% Construct scaling comparison line for left hand plot
aRaQuarter = 1.2*max(max(a_mat))+0*Ra_arr;
for i=2:length(aRaQuarter)
   aRaQuarter(i) =  aRaQuarter(i-1)*(Ra_arr(i)/Ra_arr(i-1))^(-0.25);
end


h = figure();
set(h, 'Position', [200 200 1200 500]);

m = 1; n=2;
subplot(m, n, 1);

hold on;
leg = {};

aPower = 4;

for Le_i=1:length(Le_arr)
    Le = Le_arr(Le_i);
    
    aForLe = a_mat(:, Le_i);
    hForLe = h_mat(:, Le_i);
    
     
    valid = ~isnan(aForLe);
    
    if sum(valid) < 2
        continue
    end
        
    if variableHscaling
    %plot(log10(Ra_arr(valid)),log10(hForLe(valid)./(aForLe(valid).^aPower)));
    plot(log10(Ra_arr(valid)),log10(hForLe(valid)./(aForLe(valid).^aPower)));
    
    else
       plot(log10(Ra_arr(valid)),log10((1/St)./(aForLe(valid).^aPower)));
     
    end
    
    leg{end+1} = ['Le=',num2str(Le)];
    
end

if variableHscaling
RaQuarter = plot([3 4], [3 4], ':');
leg{end+1} = '$ \sim Rm_S$';
RaQuarter = plot([3 4], [3 3.75], ':');
leg{end+1} = '$ \sim Rm_S^{3/4}$';

ylab = ['log$_{10}(h/a^',num2str(aPower),')$'];
else
    RaQuarter = plot([3 4], [2.8 3.8], ':');
leg{end+1} = '$ \sim Rm_S$';
%RaQuarter = plot([3 4], [3 3.75], ':');
%leg{end+1} = '$ \sim Rm_S^{3/4}$';

ylab = ['log$_{10}[1/(St \, a^',num2str(aPower),')]$'];
end

hold off;
legend(leg, 'Location', 'eastoutside');

axLeft = gca;

ylabLeft = ylabel(ylab); %'Rotation', 0, 'Units', 'normalized' 
xlabel('log$_{10}(Rm_S)$');
box on;

%daspect([1 1 1]);


subplot(m, n, 2);

hold on;
leg = {};

for Ra_i=1:length(Ra_arr)
    Ra = Ra_arr(Ra_i);
    aForRa = a_mat(Ra_i, :);
    hForRa = h_mat(Ra_i, :);
    
    valid = ~isnan(aForRa);
    
    if Ra == 6000 || Ra > 10000 || Ra < 1000
        continue
    end
    
    if sum(valid) < 2
        continue
    end

    if variableHscaling
    %plot(log10(Le_arr(valid)),log10(1./(aForRa(valid).^2))); %'-x'
    plot(log10(Le_arr(valid)),log10(hForRa(valid)./(aForRa(valid).^4))); %'-x'
    else       
    plot(log10(Le_arr(valid)),log10((1/St)./(aForRa(valid).^4))); %'-x'
    end
    
    leg{end+1} = ['$Rm_S$=',num2str(Ra)];
    
end

if variableHscaling
    LeZero = plot([1.9 3.9], [2.75 2.75], ':');
leg{end+1} = '$\sim Le^{0}$';
LeQuarter = plot([0.9 1.55], [2.1 2.75], ':');
leg{end+1} = '$ \sim Le$';

ylab = 'log$_{10}(h/a^4)$';
else
        LeZero = plot([1.8 3.8], [2.5 2.5], ':');
leg{end+1} = '$\sim Le^{0}$';
LeQuarter = plot([0.8 1.3], [1.9 2.4], ':');
leg{end+1} = '$ \sim Le$';

ylab = 'log$_{10}[1/(St a^4)]$';
    
end

hold off;
legend(leg, 'Location', 'eastoutside');

axRight = gca;
ylabright = ylabel(ylab); %, 'Rotation', 0, 'Units', 'normalized' 
xlabel('log$_{10}(Le)$');
box on; 
%daspect([1 0.5 1]);


% Now sort out label positions
%axLeftPos = axLeft.Position;
%ylabLeft.Position = [axLeftPos(1)-0.5 0.5 0.0];

%axRightPos = axRight.Position;
%ylabright.Position = [axRightPos(1)-0.9 0.5 0.0];

fname = 'channelWidthLeRaDarcyBrinkman';
if ~variableHscaling
    fname = [fname, 'NoH'];
end

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf, '-dpdf', [figure_output_dir, fname, '.pdf'], '-r100');

