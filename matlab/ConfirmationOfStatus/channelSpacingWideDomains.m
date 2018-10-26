close all; clear all;

data_dir = getDataDir('channelSpacing/');

folders = dir(data_dir);

%spacing = NaN*ones(length(folders), 1);

h = figure();

leg = {};

hold on;
for f_i = 1:length(folders)
   
    folder = folders(f_i);
    
    fprintf('%s \n', folder.name);
    
    folders_regex = 'CR(\d+\.?\d+?)RaC(\d+\.?\d+?)Le(\d+)(\w+)pts(\d+)-(\w+)'; %Make sure it ends in 0
    [~,matches,~]  = regexp(folder.name, folders_regex, 'match', 'tokens', 'tokenExtents');
    
    if ~isempty(matches)
        match = matches{1};
        CR = str2double(match{1}) - 1; % This is where we transform conc ratio vals
        Ra = str2double(match{2});
        Le = str2double(match{3});
        PermType = match{4};
        width = str2double(match{5});
        
        
        
        spacingFile = [data_dir, folder.name, '/channelSpacing.out'];
        
        if exist(spacingFile, 'file') ~= 2
            fprintf('No channelSpacing.out file \n');
            continue
        end
        
        try
            M = csvread(spacingFile, 1);
            
            
            t = M(:,1) - M(1,1);
            spacing = M(:,2);
            
            plot(t, spacing);
            leg{end+1} = ['Ra = ', num2str(Ra)];
            
        catch e
            fprintf('Error caught \n');
            % ignore errors
        end
        
    end
    
end
hold off;
xlabel('$t$');
ylabel('$\lambda$', 'Rotation', 0.0);
box on; 
legend(leg);