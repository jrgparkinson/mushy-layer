baseDir = '/home/parkinsonjl/mnt/sharedStorage/optimalStates-Brinkman/';

commitChanges = false;

%da_folders = dir(baseDir);


%Da_name = da_folders(i).name;

% Assume name looks like 'Da5e-5'
Da = 5e-3;
Da_name = 'Da5e-3';


CR_Ra_folders = dir([baseDir, Da_name]);

% Just do this one for now

fprintf('Processing Da = %e \n', Da);

for j=1:length(CR_Ra_folders)
    n =  CR_Ra_folders(j).name;
    fprintf('  Processing %s \n', n);
    
    % Now need to find CR and Ra
    CR_regex = 'CR(\d+\.?\d+?)RaC(\d+\.?\d+?e?\+?(\d+)?)'; %Make sure it ends in steady
    [~,matches,~]  = regexp(n, CR_regex, 'match', 'tokens', 'tokenExtents');
    if ~isempty(matches)
        match = matches{1};
        CR = str2double(match{1});
        Ra = str2double(match{2});
        
        fprintf('    CR = %1.3f, RaC=%.2e \n', CR, Ra);
        
        CRname = sprintf('CR%1.3f', CR);
        Ra_name = sprintf('RaC%.2e', Ra);
        
        % Now need to make folder to move this to, if it doesn't exist
        CR_folder = [baseDir, Da_name, '/', CRname];
        
        if exist(CR_folder, 'dir') ~= 7
            mkdir(CR_folder);
        end
        
        % Now Ra_folder bit
        Ra_folder = [CR_folder, '/', Ra_name];
        if exist(Ra_folder, 'dir') ~= 7
            mkdir(Ra_folder);
        end
        
        
        %Move this folder into the new directory
        oldFolder = [baseDir, Da_name, '/', n];
        newFolder = [Ra_folder, '/', n];
        fprintf('    Move %s to %s \n', oldFolder, newFolder);
        
        if commitChanges
            movefile(oldFolder, newFolder);
        end
        
    end
    
    
end


% Having moved all folders to the right place, now rename them
CR_folders = dir([baseDir, Da_name]);

for j=1:length(CR_folders)
    thisCR_folder = CR_folders(j).name;
    if length(thisCR_folder) < 5 || ~strcmp(thisCR_folder(1:2), 'CR')
        continue
    end
    
    Ra_folders = dir([baseDir, Da_name, '/', thisCR_folder]);
    
    for k=1:length(Ra_folders)
        thisRa_folder = Ra_folders(k).name;
        if length(thisRa_folder) < 5 || ~strcmp(thisRa_folder(1:3), 'RaC')
            
            % bug fix
            if length(thisRa_folder) >3 && strcmp(thisRa_folder(1:2), 'Ra') && ~strcmp(thisRa_folder(1:3), 'RaC')
                newName = strrep(thisRa_folder, 'Ra','RaC');
                oldDir = [baseDir, Da_name, '/', CR_folders(j).name, '/',thisRa_folder];
                newDir = [baseDir, Da_name, '/', CR_folders(j).name, '/',newName];
                
                if commitChanges
                    movefile(oldDir, newDir);
                end
                
            end
            
            continue
            
        end
        
        width_folders = dir([baseDir, Da_name, '/', CR_folders(j).name, '/',thisRa_folder]);
        
        for l=1:length(width_folders)
            thisName = width_folders(l).name;
            regex = '.+(pts.*)';
            [~,matches,~]  = regexp(thisName, regex, 'match', 'tokens', 'tokenExtents');
            if ~isempty(matches)
                match = matches{1};
                
                Nx = match{1};
                newName = Nx;
                fprintf('    rename %s to %s \n', thisName, newName);
                
                rootDir =  [baseDir, Da_name, '/', CR_folders(j).name, '/',thisRa_folder,'/'];
                oldDir = [rootDir, thisName];
                newDir = [rootDir, newName];
                
                if commitChanges
                    movefile(oldDir, newDir);
                end
                
            end
        end
        
    end
    
end










%folders_regex = ['CR(\d+\.?\d+?)RaC(\d+\.?\d+?)Le(\d+e?\.?\d+?)(\w+)pts(\d+)-',ending]; %Make sure it ends in steady
%csvFileName = 'postProcess.csv';

%[~,matches,~]  = regexp(folderName, folders_regex, 'match', 'tokens', 'tokenExtents');
%if ~isempty(matches)

