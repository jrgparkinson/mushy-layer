% Richardson error means comparing 32x32 with 64x64, 64x64 with 128x128 and
% so on. This gives an estimate of the rate of convergence.

function computeRichardsonError(output_folder, prefix, redoAnalysis)

if nargin < 3 
    redoAnalysis = true;
end

if nargin < 2
   prefix = 'Uniform-'; 
end

%output_folder = getDataDir(data_dir);

folders = dir([fullfile(output_folder, prefix), '*']);


parfor i=1:length(folders)
   folderName = folders(i).name;
   coarseDir = fullfile(output_folder, folderName);
   
   numCells = getNumCells(folderName);
   
   finerRes = numCells*2;
   finerFolder = strrep(folderName, num2str(numCells), num2str(finerRes));
   
   
   
   if exist(fullfile(output_folder, finerFolder), 'dir')
      % Compare the two
      %coarseML = getFinalPlotFile([output_folder, folderName, '/']);
      err_file_name = 'richardsonError.mat';
      err_file = fullfile(coarseDir, err_file_name);
      
      % If we already have the error file, and we don't want to redo
      % analysis, skip this folder
      if exist(err_file, 'file') && ~redoAnalysis
          continue
      end
      
      fineML = getFinalPlotFile(fullfile(output_folder, finerFolder));
      
      computeAMRerror(coarseDir, fineML, true, err_file_name, 'richardsonErr')
   end
    
end



end