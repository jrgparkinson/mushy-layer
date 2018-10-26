function compileNu

chi=0.4;
Da = '0.01';
cfl = -1;
Ra = {'1e3','1e4','1e5','5e5'};
NuLeBars = [1.01, 1.41, 3.17, 5.24];

UniformPrefix = 'Uniform-convectionDB-';

%UniformPrefix = 'Uniform-convectionDarcyBrinkman-';
%AMRPrefix = 'VariableMesh2SubcycleRefluxFreestream0.45-convectionDB-';

%prefix = 'AMRConvergenceTestConvectionDB-chi0.4-Da1e-06-Ra';
%Ra = {'1e7','1e8','1e9'};

chi='0.4';
Da = '1e-06';
cfl = '0.001';
Ra = {'1e7','1e8','1e9'};
NuLeBars = [1.08, 3.07, 12.9];
AMRPrefix = 'VariableMesh2SubcycleRefluxFreestream0.9-ref2-convectionDB-';


for r = 5:7
    
res = 2^r;

fprintf('Nx = %d \n', res);
for i=1:length(Ra)
    outputFolder = getOutputFolder(Ra{i}, chi, Da, cfl);
    
%folder= getDataDir(['AMRConvergenceTest/ConvectionDB/', params,    Ra{i}]);
folder= getDataDir(['AMRConvergenceTest/ConvectionDB/', outputFolder]);

% Get latest sub folder
j = 0;
subFolder = [folder,'/', UniformPrefix,num2str(res),'--', num2str(j)];
while exist(subFolder, 'dir') == 7
    j = j+1;
    subFolder = [folder,'/', UniformPrefix,num2str(res),'--', num2str(j)];
end
subFolder = [folder,'/', UniformPrefix,num2str(res),'--', num2str(j-1)];
[NuUniform, UniformRunning]= getNuFolder(subFolder);


%subFolder = [folder,'/',AMRPrefix,num2str(res),'-ref2--0'];

j = 0;
subFolder = [folder,'/',AMRPrefix,num2str(res),'--', num2str(j)];
while exist(subFolder, 'dir') == 7
    j = j+1;
    subFolder = [folder,'/',AMRPrefix,num2str(res),'--', num2str(j)];
end
subFolder = [folder,'/',AMRPrefix,num2str(res),'--', num2str(j-1)];
[NuAMR, AMRRunning] = getNuFolder(subFolder);


UniformStr = ' ';
AMRStr = ' ';
if UniformRunning
    UniformStr = '*';
end

if AMRRunning
    AMRStr = '*';
end


fprintf('Ra=%s, Nu=%1.2f %s (Uniform) / %1.2f%s (AMR)  / %1.2f (Le Bars & Worster) \n', ...
    Ra{i}, NuUniform, UniformStr, NuAMR, AMRStr, NuLeBars(i));
 
end

end


end



function [Nu, stillRunning] = getNuFolder(folder)

%fprintf('Processing subfolder %s \n', folder);

Nu = NaN;
stillRunning = false;

diagFileExists = false;
diagFile = [folder,'/diagnosticsLatest.out'];
if exist(diagFile, 'file') == 2
    diagFileExists = true;
else
    diagFile = [folder,'/diagnostics.out'];
    if exist(diagFile, 'file') == 2
       diagFileExists = true;
    end
end

if diagFileExists
    
    diags = getDiagnostics(diagFile);
    
%     recentMean = mean(diags.Nusselt(round(end-end/10):end));
%     recentStd = std(diags.Nusselt(round(end-end/10):end));
% 
%     recentDiff = abs(diags.Nusselt(end) - diags.Nusselt(end-1));
% 
%     err = max(recentDiff, recentStd);
    
     %fprintf('Ra=%s, Nu=%1.5f +/- %1.5f \n', Ra, diags.Nusselt(end), err);
    Nu = diags.Nusselt(end); 
    
    FileInfo = dir(diagFile);
    TimeStamp = FileInfo.datenum;
    
    ageDays = now - TimeStamp; % in days
    ageMinutes = ageDays*24*60;
    ageSeconds = ageDays*24*3600; % in seconds
    
    
    if ageMinutes > 60
        stillRunning = false;
    else
        stillRunning = true;
    end
        

end
end

function folder = getOutputFolder(Ra, chi, Da, cfl)
if cfl < 0
    folder = ['chi',chi,'-Da',Da,'-Ra',Ra];
else
    folder = ['chi',chi,'-Da',Da,'-Ra',Ra, '-cfl',cfl];
end

end