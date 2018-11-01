% Function to compute error metric for an AMR run.
% Requires a uniform mesh 'fine' solution
% Computes:
% - L1,L2, Max err
% - Run time
% 
% And saves to a file in the AMR run directory for easy future access
%

function computeAMRerror(AMRDir, HighResFile, redoRuns, saveFilename, plotDir)

% Default values:
if nargin < 2
    base_dir = '/home/parkinsonjl/mnt/sharedStorage/AMRConvergenceTest/DBVariablePorosityGaussian1proc/';
    AMRDir = [base_dir,  'AMR-Subcycle-Reflux-Freestream0.48-MaxLevel1-DBVariablePorosity-64-ref2--0/'];
    HighResDir = [base_dir, 'Uniform-DBVariablePorosity-256--0/'];
    HighResFile =   getFinalPlotFile(HighResDir);
end

if nargin < 3
    redoRuns = true;
end

if nargin < 4
    saveFilename = 'errorAnalysis.mat';
end

if nargin < 5 
    plotDir = 'errPlots';
end

fprintf('Analysing error for %s \n', AMRDir);

saveFile = fullfile(AMRDir, saveFilename);

if exist(saveFile, 'file') == 2 && ~redoRuns
    fprintf('Folder already analysed \n');
    return
end
    
% Initialise everything in case analysis fails
err = struct();
err.L1 = NaN;
err.L2 = NaN;
err.Max= NaN;
err.Sum = NaN;
runTime = NaN;

% Work out what directory the run's are stored in
dirParts = strsplit(AMRDir, '/');
finalFolder = dirParts{end};
i = 1;
while length(finalFolder) == 0
    finalFolder = dirParts{end-i};
    i = i+1;
end
runDir = strrep(AMRDir, finalFolder, '');



% Try and do the error calculation
try
   
% Compute run time
runTime = getRuntime(AMRDir);

% Load data
AMRPlotFile =   getFinalPlotFile(AMRDir);

% Check we loaded succesfully
if length(AMRPlotFile.levelArray) == 0
    fprintf('No plot file exists \n');
    return
end

% Get comparison object
exactSol = ChomboCompare(HighResFile);

% Get list of components
compNames = fieldnames(AMRPlotFile.components);

% Compute error

compsToCompute = [AMRPlotFile.components.xAdvectionvelocity, ...
   AMRPlotFile.components.Temperature, ...
   AMRPlotFile.components.xDarcyvelocity];

% Compute error for all components
compsToCompute = 1:1:length(compNames);

for comp_i = 1:length(compsToCompute)
    
    
comp = compsToCompute(comp_i); %AMRPlotFile.components.Temperature; %AMRPlotFile.components.xAdvectionvelocity; %AMRPlotFile.components.Temperature;
compName = compNames{comp}; 
fprintf('Differencing component %s \n', compName);

try

MLdiff =  exactSol.diff(AMRPlotFile, [comp], ...
                            [comp]);
                                                  
                   
[err.(compName).L1, err.(compName).L2, err.(compName).Max, err.(compName).Sum] = AMRSum(MLdiff, comp);
 
 catch e %e is an MException struct
        fprintf(1,'The identifier was:\n%s \n',e.identifier);
        fprintf(1,'There was an error! The message was:\n%s \n',e.message);
        % more error handling...
        [err.(compName).L1, err.(compName).L2, err.(compName).Max, err.(compName).Sum] = [NaN, NaN, NaN, NaN];
 
end
    
    


% Plot error                        
figSize = [100 100 1200 700];
  
                        
h = figure();
h.Visible = 'off';
set(h, 'Position', figSize);
[X,Y] = MLdiff.grid();
dat = MLdiff.dataForComp(comp);
pcolor(X,Y,dat.');

axis equal;

title(finalFolder);
box on;
c = colorbar();
c.Label.String = compNames{comp};
xlabel('x');
ylabel('y');

errPlotFolder = [runDir, '/',plotDir,'/'];
if ~exist(errPlotFolder, 'dir')
    mkdir(errPlotFolder)
end

compFolder = [errPlotFolder, compNames{comp}, '/'];
if ~exist(compFolder, 'dir')
    mkdir(compFolder)
end

errPltFile = fullfile(compFolder, [finalFolder, '.png']);
fprintf('Saving error plot to %s \n', errPltFile);
 print(h,errPltFile,'-dpng','-r150')
 


 
%  hTemp = figure();
%  set(hTemp, 'Position', figSize);
%  thisComp = AMRPlotFile.components.Temperature;
% [X,Y] = AMRPlotFile.grid();
% dat = AMRPlotFile.dataForComp(thisComp);
% pcolor(X,Y,dat.');
% 
% title(finalFolder);
% box on;
% c = colorbar();
% c.Label.String = compNames{thisComp};
% xlabel('x');
% ylabel('y');
% 
% errPltFile = [runDir, compNames{thisComp}, '-', finalFolder, '.png'];
% fprintf('Saving error plot to %s \n', errPltFile);
%  print(hTemp,errPltFile,'-dpng','-r100')
 
% Close plot window once we're done
 close(h);

end
    
  % Catch any errors that occured  
catch e %e is an MException struct
        fprintf(1,'The identifier was:\n%s \n',e.identifier);
        fprintf(1,'There was an error! The message was:\n%s \n',e.message);
        % more error handling...
end
    

% Print out the result
fprintf('Runtime: %1.3e \n', runTime);
fprintf('============= \n');

highResFilename = HighResFile.filename;
save(saveFile, 'runTime', 'err', 'highResFilename');
end
