function bifurcationPlot

close all;
clear all;

%upper_branch = getDataDir('optimalStates-Brinkman');
%lower_branch = getDataDir('optimalStates-hysteresis');
%CR = 10; Ra=3000; perm = 'KozenyPermeability';

% upper_branch = getDataDir('highLeBifurcationDarcy');
% lower_branch = getDataDir('optimalStates-hysteresis');
% CR = 1.0; Ra=200; perm = 'ChiCubedPermeability'; Le=2e10;

upper_branch = getDataDir('optimalStates-highRes-new');
lower_branch = getDataDir('optimalStates-hysteresis');
CR = 0.25; Ra=250; perm = 'ChiCubedPermeability'; Le=200;

savePlots = true;
getData = true;
image_output_folder = '/home/parkinsonjl/Documents/Talks/CommonMedia/'; %lower_branch;


if getData
    upper_data = getDat(upper_branch, CR, Ra, perm);
    lower_data = getDat(lower_branch, CR, Ra, perm);
    
    
    
end % end if get data

% Tidy up data
upper_data = tidy(upper_data, CR);
lower_data = tidy(lower_data, CR);

% Find position of maximum
maxFlux = max(upper_data(:, 2));
maxFlux_i = find(upper_data(:, 2) == maxFlux);
maxFluxLambda = upper_data(maxFlux_i, 1);


h = figure();
h.Position = [300 300 300 300];

hold on;
plot(upper_data(:, 1), upper_data(:, 2), 'x');
%plot(lower_data(:, 1), lower_data(:, 2), 'x');

hold off;


%xa = [maxFluxLambda maxFluxLambda];
%ya = [maxFlux*0.85 maxFlux*0.99];
%[xaf,yaf] = ds2nfu(xa,ya);
%annotation('textarrow', xaf,yaf,'String','(\lambda_O, F_O)');

maxF = round(max(upper_data(:, 2)), 1); 
if maxF < max(upper_data(:, 2))
    maxF = maxF +  + 0.1;
end
ax = gca;
ax.YTick = [0 maxF/2 maxF];
ax.YTickLabel = {'0', '', num2str(maxF)};
ylim([0 maxF]);
xlim([0 max(upper_data(:, 1)*1.05)]);

box on; 
xlabel('$L$');
ylabel('$F$', 'Rotation', 0.0);

if savePlots
    
set(h,'Units','Inches');
 pos = get(h,'Position');
 set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
 
 
 saveFile = [image_output_folder, '/bifurcationCR',num2str(CR),'Ra',num2str(Ra),'Le',num2str(Le)];
 
 print(h,[saveFile,'.svg'],'-dsvg','-r0')
fprintf('Saved to %s \n', saveFile);
end

end


function dat = getDat(data_dir, CR, Ra, perm)
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
            
            dat(end+1, :) = [postProcessStruct.fullWidth, postProcessStruct.flux];
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