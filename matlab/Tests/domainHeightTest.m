% Load a series of files with different heights and compare flux
% diagnostics
close all;

base_dir = ['/home/parkinsonjl/mnt/sharedStorage/', ...
    'domainHeightTestSameThetaInf/'];

katzUnits = true;
    
folders = dir([base_dir, '*-0']);

numFolders = length(folders);

domainHeight = NaN*ones(numFolders, 1);
steadyFlux = domainHeight;

for i = 1:numFolders
    folder = folders(i);
    folderName = folder.name;
    
    %poutInfo = dir([data_dir, folder.name, '/pout.0']);
    
    
    inputs = readInputs([base_dir, folderName, '/inputs']);
    diags = getDiagnostics([base_dir, folderName, '/diagnostics.out']);
    
    
    calcFlux = true;
    
    if calcFlux
        finalPlotFile = getFinalPlotFile([base_dir, folderName]);

        frameAdv = str2num(inputs.nonDimVel);
        Le = str2num(inputs.lewis);
        thisFlux = finalPlotFile.computeVerticalSoluteFlux(frameAdv, Le, katzUnits);

    else
        
        thisFlux = NaN;
        if isfield(diags, 'Fs_vertical_av')
            thisFlux = diags.Fs_vertical_av(end);
        elseif isfield(diags, 'Fs_bottom')
            thisFlux = diags.Fs_bottom(end);
        end
        
        if katzUnits
            thisFlux = -1.0 - thisFlux;
        else
            thisFlux = -thisFlux;
        end

    end
    
    numCells = inputs.num_cells;
    cells = strsplit(numCells, ' ');
    Nx = str2num(cells{1});
    Nz = str2num(cells{2});
    
    if isfield(inputs, 'domain_width')
        domainWidth = str2num(inputs.domain_width);
    else
        domainWidth = str2num(inputs.domain_length);
    end
    
    steadyFlux(i) = thisFlux;
    domainHeight(i) = domainWidth*Nz/Nx;
    
end

figure();
plot(domainHeight, steadyFlux, '-x');
xlabel('$H$');
ylabel('$F$');
title('Salt flux dependence on domain height');

savefig('domainHeightTest.fig');