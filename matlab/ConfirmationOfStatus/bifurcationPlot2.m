function bifurcationPlot2

close all;
clear all;


% Darcy:
sets(1).upper = 'highLeBifurcationDarcy';
sets(1).lower = 'optimalStates-hysteresis';
sets(1).Le = 2e10;
sets(1).CR = 1.0; sets(1).Ra = 200; sets(1).perm = 'ChiCubedPermeability';


sets(2).upper = 'optimalStates-highRes-new';
sets(2).lower = 'optimalStates-hysteresis';
sets(2).Le = 200;
sets(2).CR = 1.0; sets(2).Ra = 200; sets(2).perm = 'ChiCubedPermeability';

plotName = 'bifurcationComparisonDarcy.pdf';

% These runs keep crashing

% Darcy-brinkman:
% sets(1).upper = 'highLeBifurcationDB';
% sets(1).lower = 'optimalStates-hysteresis';
% sets(1).Le = 2e10;
% sets(1).CR = 1.0; sets(1).Ra = 10000; sets(1).perm = 'KozenyPermeability';
% 
% 
% sets(2).upper = 'optimalStates-Brinkman';
% sets(2).lower = 'optimalStates-hysteresis';
% sets(2).Le = 200;
% sets(2).CR = 1.0; sets(2).Ra = 10000; sets(2).perm = 'KozenyPermeability';
% 
% 
% sets(1).upper = 'highLeBifurcationDB';
% sets(1).lower = 'optimalStates-hysteresis';
% sets(1).Le = 2e10;
% sets(1).CR = 5.0; sets(1).Ra = 4000; sets(1).perm = 'KozenyPermeability';
% 
% 
% sets(2).upper = 'optimalStates-Brinkman';
% sets(2).lower = 'optimalStates-hysteresis';
% sets(2).Le = 200;
% sets(2).CR = 5.0; sets(2).Ra = 4000; sets(2).perm = 'KozenyPermeability';

diagnosticToPlot = 'maxStreamfunctionMush';
diagnosticToPlot = 'flux';

diagnosticToPlotLabel = '$F$';


h = figure();
h.Position = [300 300 400 300];

hold on;

for set_i = 1:length(sets)
    upper_branch = getDataDir( sets(set_i).upper);
lower_branch = getDataDir( sets(set_i).lower);
CR =  sets(set_i).CR; Ra= sets(set_i).Ra; perm =  sets(set_i).perm; 
Le= sets(set_i).Le;


   upper_data = getDat(upper_branch, CR, Ra, perm, diagnosticToPlot);
    lower_data = getDat(lower_branch, CR, Ra, perm, diagnosticToPlot);
    
 
% Tidy up data
upper_data = tidy(upper_data, CR);
lower_data = tidy(lower_data, CR);

% Find position of maximum
%maxFlux = max(upper_data(:, 2));
%maxFlux_i = find(upper_data(:, 2) == maxFlux);
%maxFluxLambda = upper_data(maxFlux_i, 1);


try
plot(upper_data(:, 1), upper_data(:, 2), 'x');
%plot(lower_data(:, 1), lower_data(:, 2), 'x');
catch e
    
end




end

hold off;



% maxF = round(max(upper_data(:, 2)), 1) + 0.1; 
ax = gca;
ax.YTick = [0 0.3 0.6];
ax.YTickLabel = {'0', '', '0.6'};
% ylim([0 maxF]);
% xlim([0 max(upper_data(:, 1)*1.05)]);

box on; 
xlabel('$\lambda$');
ylab = ylabel(diagnosticToPlotLabel);
prevPos = ylab.Position;
ylab.Position = [prevPos(1)+0.03 prevPos(2)];
 legend({'Le=2 $\times 10^{10}$', 'Le=200'}, 'Location', 'southeast');

 
 set(h,'Units','Inches');
  pos = get(h,'Position');
  set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
  %print(h,[lower_branch, '/bifurcationCR',num2str(CR),'Ra',num2str(Ra),'Le',num2str(Le),'.pdf'],'-dpdf','-r0')
  print(h,[lower_branch, '/', plotName],'-dpdf','-r0')
  
  fprintf('Saved to %s \n', [lower_branch, '/', plotName]);


end


function dat = getDat(data_dir, CR, Ra, perm, field)
if nargin < 5
    field = 'flux';
end

dat = [];

folders = dir(data_dir);
for f_i = 1:length(folders)
    folder = folders(f_i);
    
    folders_regex = 'CR(\d+\.?\d+?)RaC(\d+\.?\d+?)Le(\d+e?\d+)(\w+)pts(\d+)-(\w+)'; %Make sure it ends in 0
    [~,matches,~]  = regexp(folder.name, folders_regex, 'match', 'tokens', 'tokenExtents');
    
    if ~isempty(matches)
        match = matches{1};
        thisCR = str2double(match{1}) - 1; % This is where we transform conc ratio vals
        thisRa = str2double(match{2});
        thisLe = str2double(match{3});
        thisPermType = match{4};
        thiswidth = str2double(match{5});
        
        if CR == thisCR && Ra == thisRa && strcmp(thisPermType, perm)
            % Get the flux and width
            filename = [data_dir, '/', folder.name, '/postProcess.csv'];
            
            if ~exist(filename, 'file')
                fprintf('  No postProcess.csv \n');
                continue
            end
            
            A = importdata(filename, ',', 1);
            
            postProcessStruct = cell2struct(num2cell(A.data), A.colheaders, 2);
            
            dat(end+1, :) = [postProcessStruct.fullWidth, postProcessStruct.(field)];
        end
    end
    
    
end

end


function data = tidy(data, CR)
% 1) Sort
data = sortrows(data);

% 2) If flux is less than this we assume it's 0.
minFlux = 1e-4*CR;

try
for i = length(data)-1:-1:1
    if data(i, 2) < minFlux ||  data(i+1, 2) < minFlux
        data(i, 2) = 0.0;
    end
end

% 3) Extend to lambda = 0

if data(1, 2) < minFlux && data(1, 1) > 0.0
    lambda =  data(1, 1);
    dlambda = data(2, 1) - data(1, 1);
    if dlambda > 0
        while lambda >= 0
            data(end+1, :) = [lambda, 0.0];
            lambda = lambda - dlambda;
        end
    end
end
catch e
    
end
% 4) Sort again
data =  sortrows(data);


end