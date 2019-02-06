function analyseVariablePorosityTest(base_dir, Nzs, redoAnalysis, runAnalysis, uniform_prefix,...
    compName, errType)
close all;


if nargin < 1
    fprintf('**Warning - the directory containing the data has not been specified** \n');
    base_dir = 'FixedPorousHole';
    base_dir = '/home/parkinsonjl/mnt/sharedStorage/TestDiffusiveTimescale/FixedPorousHole-1proc/';
end

if nargin < 2
    Nzs = [16,32,64,128];
end

if nargin < 3
    redoAnalysis = false;
end

if nargin < 4
    runAnalysis = false;
end

if nargin < 5
    uniform_prefix = 'Uniform-DBVariablePorosity-';
end

if nargin < 6
   compName = 'Porosity';
end

if nargin < 7
    errType = 'L2';
    % other options: L1, Max
end

[~, name] = system('hostname'); 
name = strtrim(name);

if strcmp(name, 'atmlxlap005')

% For running the fixed porous test problem on laptop:
base_dir = '/home/parkinsonjl/mnt/sharedStorage/TestDiffusiveTimescale/FixedPorousHole-1proc/';
Nzs = [16,32,64,128,256,512];
redoAnalysis = false;
runAnalysis = false;
uniform_prefix = 'Uniform-DBVariablePorosity-';
compName  = 'xDarcyvelocity';
errType = 'L2';

% base_dir = '/home/parkinsonjl/mnt/sharedStorage/TestDiffusiveTimescale/PorousMushyHole-t0.00015/';
% Nzs = [16,32,64,128,256,512];
% redoAnalysis = false;
% runAnalysis = false;
% uniform_prefix = 'Uniform-PorousMushyHole-';
% compName = 'Porosity';
% errType = 'L2';

end

fineNumCells = Nzs(end);
fine_res_dir = [uniform_prefix,num2str(fineNumCells),'--0'];

% For uniform grids
if runAnalysis
    computeRichardsonError(base_dir, uniform_prefix)
end


folders = dir(base_dir);

for j=1:length(Nzs)-1
   Nz = Nzs(j);
    
   thisErrStruct = struct();
   
   thisErrStruct.numCells = Nz;
   thisErrStruct.singleLevel = NaN;
   thisErrStruct.richardson = NaN;
   thisErrStruct.rate = NaN;
   thisErrStruct.nref2 = NaN;
   thisErrStruct.nref4 = NaN;
   thisErrStruct.nref22 = NaN;
      
   errStruct(j) = thisErrStruct;
end

if runAnalysis
    tic;
    fprintf('Loading high res file \n');
    HighResFile = getFinalPlotFile(fullfile(base_dir, fine_res_dir));
    fprintf('Loaded high res file (%1.3f seconds) \n', toc);

    tic;
    
   
    for i=1:length(folders)
        folder = folders(i);
        folderName = folder.name;
        
        fullFolder = fullfile(base_dir, folderName);

        if exist(fullFolder, 'dir') ~= 7
            continue
        end

        if length(folderName) > 3 && ~strcmp(folderName, fine_res_dir) ...
                && ~strcmp(folderName, 'errPlots') ...
                && ~strcmp(folderName, 'richardsonErr')
            fprintf('Considering %s \n', folderName);
            skip = false;
            thisNumCells = getNumCells(folderName);
            if thisNumCells >= fineNumCells
                continue
            end

            fprintf('%s \n', folder.name);
            computeAMRerror(fullFolder, HighResFile, redoAnalysis);
        end
    end

    fprintf('Analysed all folders (%1.3f seconds) \n', toc);

end

% Now that all files have been analysed, plot the results
runs = getRuns(base_dir);

thisRunTime = NaN*ones(length(runs), 10);
    thisError = thisRunTime;
    thisNumCells = thisRunTime;

for i=1:length(runs)
    folders = dir([fullfile(base_dir, runs{i}), '*']);
    
    for j=1:length(folders)
        folder_dir = fullfile(base_dir, folders(j).name);
        fprintf('%s \n', folders(j).name);
        
        errorFile = fullfile(folder_dir, 'errorAnalysis.mat');
        if exist(errorFile, 'file') ~= 2
            continue
        end
        
        try
        
        load(errorFile);
        
        numCells = getNumCells(folders(j).name);
        
        
        
        ref = 1;
        maxLev = 0;
       
        
        [refTokens,refMatches] = regexp(folders(j).name,'ref(\d)','tokens','match');
        [levTokens,levMatches] = regexp(folders(j).name,'MaxLevel(\d)','tokens','match');
        
        if length(refMatches) > 0
            ref = str2num(refTokens{1}{1});
        end
        
        if length(levMatches) > 0
            maxLev = str2num(levTokens{1}{1});
        end
        
        totalRefinement = ref^maxLev;
        
        if numCells*totalRefinement >= fineNumCells
            continue
        end
        
        %exponent = 1+log(numCells*totalRefinement)/log(2) - log(Nzs(1))/log(2);
        % Indexing is based on base level
        exponent = 1+log(numCells)/log(2) - log(Nzs(1))/log(2);
        
        if exponent < 1
            fprintf('Skipping Nz < min(Nz specified) - skipping');
            continue;
        end
        
        
        if numCells*totalRefinement == fineNumCells/2 && maxLev<2 ...
                && (length(strfind(folders(j).name, 'AMR') ) > 0 || length(strfind(folders(j).name, 'Uniform') ) > 0)
            
            pout = Pout(fullfile(folder_dir, 'pout.0'));
            
            ncells = pout.pointsUpdated;
            
            index = log(ref)/log(2) + 1;
            if ref == 1
           performance(index).n_ref = 0;
            else
                performance(index).n_ref = ref;
            end
           performance(index).ncells = ncells;
           performance(index).runtime = runTime;
           performance(index).folder = folders(j).name;
        end
        
        
        
        
        thisRunTime(i, j) = runTime;
        if isfield(err, compName)
            thisError(i, j) = err.(compName).(errType);
        end
        thisNumCells(i, j) = numCells;
        
        
        thisErrStruct = errStruct(exponent);
        
        if strfind(folders(j).name, 'Uniform')
           thisErrStruct.singleLevel = thisError(i, j);
           
           
           % Also get richardson error
            richardsonErrorFile = fullfile(folder_dir, 'richardsonError.mat');
            if exist(richardsonErrorFile, 'file') ~= 2
                continue
            end

            load(richardsonErrorFile);
            thisErrStruct.richardson = err.(compName).(errType);
           
        elseif length(strfind(folders(j).name, 'AMR')) > 0 ...
                &&  length(strfind(folders(j).name, 'ref2')) > 0  ...
                &&  length(strfind(folders(j).name, 'MaxLevel1')) > 0
            thisErrStruct.nref2 = thisError(i, j);
            
        elseif  length(strfind(folders(j).name, 'AMR') ) > 0 ...
                &&  length(strfind(folders(j).name, 'ref2') ) > 0 ...
                &&  length(strfind(folders(j).name, 'MaxLevel2')) > 0
            thisErrStruct.nref22 = thisError(i, j);
            
        elseif  length(strfind(folders(j).name, 'AMR')  ) > 0 ...
                &&  length(strfind(folders(j).name, 'ref4')) > 0
                
            thisErrStruct.nref4 = thisError(i, j);
        end
        
        
        
        % Save struct
        errStruct(exponent) = thisErrStruct;
        
                      % Catch any errors that occured  
            catch e %e is an MException struct
                    fprintf(1,'The identifier was:\n%s \n',e.identifier);
                    fprintf(1,'There was an error! The message was:\n%s \n',e.message);
                    % more error handling...
            end
        
    end
    
    
    
end


makeSummaryFig = false;
if makeSummaryFig
    h = figure();
    set(h, 'Position', [100 100 1400 800]);


    marker = {};
    for i=1:length(runs)
        marker{i} = 'x';

        if runs{i}(1:3) == 'AMR'
            marker{i} = 'o';
        elseif runs{i}(1:8) == 'Variable'
            marker{i} = '+';
        end
    end



    m = 2; n=2;

    subplot(m, n, 1)

    hold on;
    for i=1:length(runs)
    plot(log10(thisNumCells(i, :)), log10(thisRunTime(i, :)), marker{i});
    end
    hold off;
    box on;
    xlabel('num cells'); ylabel('runtime');

    subplot(m, n, 2)

    hold on;
    for i=1:length(runs)
    plot(log10(thisNumCells(i, :)), log10(thisError(i, :)), marker{i});
    end
    hold off;
    box on;
    xlabel('num cells'); ylabel('error');



    subplot(m, n, [3 4])

    hold on;
    for i=1:length(runs)
    plot(log10(thisRunTime(i, :)), log10(thisError(i, :)), marker{i});
    end
    hold off;

    xlabel('runtime'); ylabel('error');


    runsShort = {};
    for i=1:length(runs)
       runsShort{i} = runs{i}(1:10); 
    end

    legend(runs, 'Location', 'eastoutside');
    box on;
end


extrapolateRunTime = false;

if extrapolateRunTime
    figure();
    hold on;
    for i=1:length(runs)
        ncell = log10(thisNumCells(i, :));
        runt = log10(thisRunTime(i, :));
        validVals = 1-isnan(ncell).*isnan(runt);


        extrapncell = sort([ncell(validVals==1) max(ncell) max(ncell)*1.5 max(ncell)*2 max(ncell)*3]);

        P = polyfit(sort(ncell(validVals==1)), sort(runt(validVals==1)), 2);

        extraprun = P(3) + extrapncell.^2*P(1) +  extrapncell*P(2);

        %extrapncell = [ncell(validVals==1) max(ncell) max(ncell)*3];

        %extraprun = interp1(ncell(validVals==1), runt(validVals==1), extrapncell, 'linear', 'extrap');

        plot(log10(thisNumCells(i, :)), log10(thisRunTime(i, :)), marker{i});
        plot(extrapncell, extraprun, '-');

    end
    hold off;
    box on;
    xlabel('num cells'); ylabel('runtime');

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make figure version of tables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(h,'Units','Inches');
h.Position = [2.0 2.0 6.0 2.5];
textFontSize = 10;
legendFontSize = 8;
domainFontSize = 8;

set(0, 'defaultlinelinewidth',1);
set(0, 'defaultaxeslinewidth',1);
set(0, 'defaultpatchlinewidth',1);
set(0, 'defaultAxesFontSize', textFontSize);

pltBottom = 0.18;
pltHeight = 0.77;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


errors = [];
num_cells = [];

for j=1:length(errStruct)
    if j > 1
        errStruct(j).rate = ( errStruct(j-1).richardson /  errStruct(j).richardson ) * ...
            (errStruct(j-1).numCells / errStruct(j).numCells);
   
    end
    
  % s = sprintf(format, errStruct(j).numCells, errStruct(j).richardson, errStruct(j).rate, errStruct(j).singleLevel, ...
   %    errStruct(j).nref2, errStruct(j).nref4, errStruct(j).nref22);
   
   num_cells(end+1) = errStruct(j).numCells;
   errors(end+1, :) = [ errStruct(j).richardson, errStruct(j).singleLevel,  ...
       errStruct(j).nref2, errStruct(j).nref4, errStruct(j).nref22];

end

styles = {'x-', 'x-', 's-', 'd-', 'o-'};

hold on;
for i=1:size(errors, 2)
plot(num_cells, errors(:, i), styles{i});
end

ncells_second_order = [num_cells, num_cells(end)*1.4];
second_order = 1.2 * max(errors(1, :))*(ncells_second_order/num_cells(1)).^(-2) ;
plot(ncells_second_order, second_order, '--');

hold off;
box on;

axLeft = gca;
axLeft.YScale = 'log';
axLeft.XScale = 'log';

axLeft.Position = [0.1 pltBottom 0.48 pltHeight];
%ax.Position(3) = ax.Position(3) + 0.05;

axLeft.XLim(2) = axLeft.XLim(2)*2;
axLeft.YLim(2) = axLeft.YLim(2)*2;

ylab = [format_err_type(errType), ' error (', nice_comp_name(compName), ')'];

xlabel('$1/\Delta x$');
ylabel(ylab);
legend({'Single-level Richarson', 'Single-level 512 difference', '$n_{ref}=2$', '$n_{ref}=4$', '$n_{ref}=(2,2)$', '2nd order'}, ...
    'FontSize', legendFontSize,...
    'Position', [axLeft.Position(1)+axLeft.Position(3)-0.15, axLeft.Position(2)+pltHeight-0.2, 0.01, 0.01]);

xp = axLeft.XLim(1)*0.5;
yp = axLeft.YLim(2);
text(xp, yp, '(a)');



axRight = axes;

normalizedCells = [];
normalizedCPU = [];
nref = [0, 2, 4];

for j=1:length(performance)
    normalizedCells(end+1) = performance(j).ncells/performance(1).ncells;
    normalizedCPU(end+1) = performance(j).runtime/performance(1).runtime;
    
   %fprintf(format, performance(j).n_ref, performance(j).ncells, performance(j).runtime, ...
   %    normalizedCells, normalizedCPU);
end

hold on;
plot(nref, normalizedCells, '-x');
plot(nref, normalizedCPU, '-x');
hold off;

box on;

axRight.Position = [0.66 pltBottom 0.32 pltHeight];

xlabel('Refinement ratio');
legend('Normalized cells', 'Normalized time');
%title('(b)');
text(-0.8, 1.02, '(b)');

if contains(uniform_prefix, 'VariablePorosity')
    fig_num = 6;
else
    fig_num = 8;
end

figureName = fullfile(base_dir, ['Fig', num2str(fig_num), uniform_prefix, compName, errType]);
fprintf('Saved to %s \n', figureName);
print(h,[figureName, '.eps'],'-depsc','-r50')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print tables, first latex then human readable format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('%s error computed for the field: %s \n', errType, compName);

% Error
fprintf('%5s | %15s  | %5s  | %13s | %12s  | %12s  | %12s \n', '1/dz', 'Richardson', 'Rate', '512 diff', 'n_{ref} = 2', 'n_{ref} = 4', 'n_{ref} = (2,2)');

format = '%5d & %15.2e  & %1.2f  &  %13.2e & %12.2e  & %12.2e  & %12.2e \\\\ \n';
fullStr = '';
for j=1:length(errStruct)
    if j > 1
        errStruct(j).rate = ( errStruct(j-1).richardson /  errStruct(j).richardson ) * ...
            (errStruct(j-1).numCells / errStruct(j).numCells);
   
    end
    
   s = sprintf(format, errStruct(j).numCells, errStruct(j).richardson, errStruct(j).rate, errStruct(j).singleLevel, ...
       errStruct(j).nref2, errStruct(j).nref4, errStruct(j).nref22);
   
   fullStr = [fullStr, s];
end


fullStr = strrep(fullStr, 'NaN', '-');
fprintf('%s', fullStr);


% Performance
fprintf('%7s | %15s  | %15s | %20s | %15s  \n', 'n_{ref}', 'Cells advanced', 'CPU time (s)',...
    'Normalized cells advanced', 'Normalized CPU time');

format = '%7d & %15d  & %15.0f & %1.4f & %1.4f \\\\ \n';
for j=1:length(performance)
    normalizedCells = performance(j).ncells/performance(1).ncells;
    normalizedCPU = performance(j).runtime/performance(1).runtime;
    
   fprintf(format, performance(j).n_ref, performance(j).ncells, performance(j).runtime, ...
       normalizedCells, normalizedCPU);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print out table

% Error
fprintf('%5s | %15s  | %5s  | %13s | %12s  | %12s  | %12s \n', '1/dz', 'Richardson', 'Rate', '512 diff', 'n_{ref} = 2', 'n_{ref} = 4', 'n_{ref} = (2,2)');

format = '%5d | %15.2e  | %1.2f  |  %13.2e | %12.2e  | %12.2e  | %12.2e \n';
for j=1:length(errStruct)
    if j > 1
        errStruct(j).rate = ( errStruct(j-1).richardson /  errStruct(j).richardson ) * ...
            (errStruct(j-1).numCells / errStruct(j).numCells);
    
    end
   fprintf(format, errStruct(j).numCells, errStruct(j).richardson, errStruct(j).rate, errStruct(j).singleLevel,  errStruct(j).nref2, errStruct(j).nref4, errStruct(j).nref22);; 
end

% Performance
fprintf('%7s | %15s  | %15s | %20s | %15s  \n', 'n_{ref}', 'Cells advanced', 'CPU time (s)',...
    'Normalized cells advanced', 'Normalized CPU time');

format = '%7d | %15d  | %15d | %1.4f | %1.4f | %50s \n';
for j=1:length(performance)
    normalizedCells = performance(j).ncells/performance(1).ncells;
    normalizedCPU = performance(j).runtime/performance(1).runtime;
    
   fprintf(format, performance(j).n_ref, performance(j).ncells, performance(j).runtime, ...
       normalizedCells, normalizedCPU, ...
       performance(j).folder);
end




end


function formatted = format_err_type(errType)

if errType(1) == 'L'
    formatted = ['$L_', errType(2), '$'];
else
    formatted = errType;
end

end

function formatted = nice_comp_name(compName)

if compName(1) == 'x' || compName(1) == 'y'
    formatted = compName;
    
    if strcmp(compName, 'xDarcyvelocity')
        formatted = '$x-$velocity';
    end
else
 formatted = compName;
end

end