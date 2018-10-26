% Compare optimal channel spacing with wide domain results
clear all; close all;

CR = 0.25;
Da = 1.0;
RA_MAX = 400;
RA_MIN = 70;

optFolder = getDataDir('optimalStates-highRes-new/');
wideFolder = getDataDir('channelSpacing/');

load([optFolder, 'optimalVals.mat']);

% Get widths for this CR
k = keys(optimalVals);
Rm_width_opt = [];
Rm_width_wide = [];

for CR_i = 1:length(k)
    keyPair = k(CR_i);
    keyPair = keyPair{1};
    temp = keyPair(1);  thisCR  = temp{1};
    temp = keyPair(2);  thisRa = temp{1};
    
    
    if  thisCR == CR && thisRa >= RA_MIN
        
        optVals = optimalVals(thisCR, thisRa);
        %width_arr(end+1) =
        %width = 1.0/thisRa;
        width = optVals.fullWidth;
        if ~isnan(width)
        Rm_width_opt(end+1, :) = [thisRa width];
        end
    end

end


% Now get wide domain data
folders = dir(wideFolder);
for f_i = 1:length(folders)
   
    folder = folders(f_i);
    
    fprintf('%s \n', folder.name);
    
    folders_regex = 'CR(\d+\.?\d+?)RaC(\d+\.?\d+?)Le(\d+)(\w+)pts(\d+)-(\w+)'; %Make sure it ends in 0
    [~,matches,~]  = regexp(folder.name, folders_regex, 'match', 'tokens', 'tokenExtents');
    
    if isempty(matches)
        fprintf('Could not match regex\n');
    else
        match = matches{1};
        thisCR = str2double(match{1}) - 1; % This is where we transform conc ratio vals
        thisRa = str2double(match{2});
        Le = str2double(match{3});
        PermType = match{4};
        width = str2double(match{5});
        
        if CR ~= thisCR || thisRa > RA_MAX
            fprintf('Not in the parameter range we are considering, skip\n');
            continue;
        end
        
        spacingFile = [wideFolder, folder.name, '/channelSpacing.csv'];
        
        if exist(spacingFile, 'file') ~= 2
            fprintf('No channelSpacing.csv file \n');
            continue
        end
        
        try
            %M = csvread(spacingFile, 1);
            %spacing = M(end,2);
           
            T = readtable(spacingFile);
             
            if ismember('spacing', T.Properties.VariableNames)

                spacing = T.spacing(end);
                
                spacing = nanmean(T.spacing);
                %spacing  = nanmedian(T.spacing);

                domWidth = (width/1024)*1.6;
                Nchan = round(domWidth/spacing);
                err = abs(domWidth/(Nchan-1) - spacing);
                
                % increase error to account for varying channel spacing
                %err = err*(1+nanstd(T.spacing));
                
                widthVariation = max(T.spacing) - min(T.spacing);
                
                %err = err + widthVariation/2;
                %err = widthVariation/2;
                
                %spacing = (max(T.spacing) + min(T.spacing))/2;
                
                %spacing = min(T.spacing);
                
                err = (nanstd(T.spacing))/2;

                Rm_width_wide(end+1, :) = [thisRa, spacing, err];
                
                fprintf('Got channel spacing = %1.5f \n', spacing);
            end
            
        catch e
            fprintf(1, 'There was an error: \n %s', e.message);
            %fprintf('Error caught \n');
            % ignore errors
        end
        
    end
    
end

Rm_width_wide = sortrows(Rm_width_wide);


h = figure();
h.Position = [200 200 700 400];
hold on;
plot(Rm_width_opt(:, 1), Rm_width_opt(:, 2));
errorbar(Rm_width_wide(:, 1), Rm_width_wide(:, 2), Rm_width_wide(:, 3), ...
 'x', 'LineWidth', 2.0);
hold off;

xlabel('$Rm$'); ylabel('$\lambda$', 'Rotation', 0.0);
legend({'Optimal widths', 'Wide domain'});

box on;


set(h,'Units','Inches');
 pos = get(h,'Position');
 set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
 print(h,[wideFolder, '/channelSpacingComparison.pdf'],'-dpdf','-r0')
