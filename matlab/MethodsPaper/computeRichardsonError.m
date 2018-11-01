% Richardson error means comparing 32x32 with 64x64, 64x64 with 128x128 and
% so on. This gives an estimate of the rate of convergence.

function computeRichardsonError(data_dir, prefix)

if nargin < 1
    data_dir = '/AMRConvergenceTest/DBVariablePorosityGaussian1proc-t1.6-v2/';
end

if nargin < 2
   prefix = 'Uniform-DBVariablePorosity-'; 
end

output_folder = getDataDir(data_dir);

folders = dir([output_folder, prefix, '*']);

parfor i=1:length(folders)
   folderName = folders(i).name;
   coarseDir = fullfile(output_folder, folderName);
   
   numCells = getNumCells(folderName);
   
   finerRes = numCells*2;
   finerFolder = strrep(folderName, num2str(numCells), num2str(finerRes));
   
   
   
   if exist(fullfile(output_folder, finerFolder), 'dir')
      % Compare the two
      %coarseML = getFinalPlotFile([output_folder, folderName, '/']);
      fineML = getFinalPlotFile(fullfile(output_folder, finerFolder));
      
      computeAMRerror(coarseDir, fineML, true, 'richardsonError.mat', 'richardsonErr')
   end
    
end



end