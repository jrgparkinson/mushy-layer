% Load all runs with the directory specified and compute various
% diagnostics, writing the output into a csv file postProcess.out

% To run from command line:
%{ 
matlab -nodisplay -nosplash -nodesktop -r "postProcessAllFolders"

matlab -nodisplay -nosplash -nodesktop -r "cd ~/convection-in-sea-ice/MushyLayer/matlab/PostProcess/; 
postProcessAllFolders('/network/group/aopp/oceans/AW002_PARKINSON_MUSH/optimalStates-highRes-new/'); exit;"

%}

function postProcessAllFolders(data_dir, katzUnits, makePlots, redoRuns, requireSteady)

if nargin < 5
     requireSteady = true;
end

if nargin < 4  
    redoRuns = false; 
end 

if nargin < 3 
    makePlots = false;
end

if nargin < 2 
    katzUnits = true;
end

if nargin < 1
    data_dir = ['/home/parkinsonjl/mnt/sharedStorage/', ...
    'optimalStates-highRes-new/'];
end

if makePlots
    close all;
end

%data_dir = ['/home/parkinsonjl/mnt/sharedStorage/', ...
%    'optimalStates-highRes-Wells/'];
%katzUnits = false;
ending = '*';
if requireSteady
    ending = 'steady';
end

folders = dir([data_dir, '*-',ending]);

% Get optimal folders
if exist([data_dir, 'optimalFolders.csv'], 'file')
      optimalFolders=textread([data_dir, 'optimalFolders.csv'],'%s');
else
    optimalFolders = [];
end

%foldersToDelete = folders.name;
foldersToDelete = [];
%foldersToDelete = optimalFolders;

 for i = 1:length(foldersToDelete)
    folder = foldersToDelete(i);
    fname = [data_dir, '/', folder{1}, '/postProcess.csv'];
    if exist(fname, 'file')
        delete(fname);
        fprintf('Deleted %s \n', fname);
    else
        fprintf('File does not exist: %s \n', fname);
    end
     %delete fname
 end



numFolders = length(folders);
foldersProcessed = 0;
foldersLeft = numFolders;

tic;

for i = 1:length(folders)
    folder = folders(i);
    
    
    % Skip folders not in steady state
    % By definition should be steady state
    steadyState = false;
    poutInfo = dir([data_dir, folder.name, '/pout.0']);
    now = datetime('now');
    
    try 
        %if abs(datenum(now) - poutInfo.datenum) > 0.001
        folder_parent_dir = [];
        
            postProcessFolder(folder.name, makePlots, redoRuns, folder_parent_dir, katzUnits, requireSteady);
        %else
        %    fprintf([folder.name, ' has not reached steady state \n']);
       % end
    
    catch e
       fprintf('Could not process folder, \n%s \n',e.identifier);
       fprintf(1,'Error message:\n%s \n',e.message); 
    end
    
    timeElapsed = toc;
    foldersProcessed = foldersProcessed + 1;
    foldersLeft = foldersLeft  - 1;
    
    avTimePerFolder = timeElapsed/foldersProcessed;
    timeLeft = avTimePerFolder*foldersLeft;
    
    remainingMins = floor(timeLeft/60);
    remainingSecs = round(timeLeft - remainingMins*60);
    
    fprintf('%d folders left (%d total), elapsed time %1.5fs, estimated time left %1.5fs (%d min %d s) \n', ...
        foldersLeft, numFolders, timeElapsed, ...
        timeLeft, remainingMins, remainingSecs);
    
     if makePlots
        pause;
    end
    
end

end