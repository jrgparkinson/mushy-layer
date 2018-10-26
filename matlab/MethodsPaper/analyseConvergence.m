function analyseConvergence(folder)
close all;

%compName = 'Temperature';
%compName = 'xAdvectionvelocity';
%compName = 'xDarcyvelocity';
compName = 'Porosity';
%compName = 'Liquidconcentration';

errType = 'L2';
%errType = 'Max';

defaultFolder = 'MushyConvection';
defaultFolder = 'MushyConvectionLiquid';
defaultFolder = 'MushyConvectionLiquid-t0.1';

defaultFolder = 'MushyConvectionLiquid2-t0.1';
defaultFolder = 'MushyConvectionLiquid2-t0.2'; % longer time
defaultFolder = 'MushyConvectionLiquid2-t0.11'; % smaller dt
defaultFolder = 'MushyConvectionLiquid2-t0.12'; % smaller rayleigh
defaultFolder = 'MushyConvectionLiquid2-t0.13'; % bigger hole
defaultFolder = 'MushyConvectionLiquid2-t0.16'; % more diffusion
defaultFolder = 'MushyConvectionLiquid2-t0.08'; % all porous
defaultFolder = 'MushyConvectionLiquid2-t0.04'; % all porous, smaller dt

%defaultFolder = 'MushyConvectionLiquid2-t0.0'; % check initial conditions are second order
%defaultFolder = 'MushyConvectionLiquid2-t0.001'; % check 1 timestep
%defaultFolder = 'MushyConvectionLiquid2-t0.005'; % check 5 timesteps
%defaultFolder = 'MushyConvectionLiquid2-t0.05'; % check 50 timesteps
%defaultFolder = 'MushyConvectionLiquid2-t0.1'; % check 50 larger timesteps - works

defaultFolder = 'MushyConvectionLiquid2-t0.1'; % shorter - good
%defaultFolder = 'MushyConvectionLiquid2-t1.0'; % bad
defaultFolder = 'MushyConvectionLiquid2-t0.2'; % middle dist - good
defaultFolder = 'MushyConvectionLiquid2-t0.4'; % longer - good
defaultFolder = 'MushyConvectionLiquid2-t0.6'; % even longer - not great
defaultFolder = 'MushyConvectionLiquid2-t0.5'; % smaller dt than 0.6, seeing if this helps - good
defaultFolder = 'MushyConvectionLiquid2-t0.5-1proc'; % for timing
%defaultFolder = 'MushyConvectionLiquid2-t0.8'; % smaller dt than 0.5 + longer - trying to get regridding. Does it work ??



if nargin < 1
    %  folder = '/AMRConvergenceTest/DBVariablePorosityGaussian1proc-t1.6-v2/';
    %   folder = '/AMRConvergenceTest/MushyDB/';
    %   folder = '/AMRConvergenceTest/MushQuick/';
    
    folder = defaultFolder;
end

base_dir = getDataDir(['/AMRConvergenceTest/',folder,'/']);
folders = dir(base_dir);
runs = getRuns(base_dir);

CompareDir = [base_dir, 'ChomboCompare/'];

% Get the finest uniform resolution
fineNumCells = 1;
lowestNumCells = 1024;
for j = 1:length(folders)
    if strfind(folders(j).name, 'Uniform')
        thisRes = getRes(folders(j).name);
        fineNumCells = max(fineNumCells, thisRes);
        lowestNumCells = min(lowestNumCells, thisRes);
    end
end

% hard code for now
fineNumCells = 512;


% Restructure data for table
% consider 2^4 to 2^8
numCellsOffset = (log(lowestNumCells)/log(2) -1 );

for j=numCellsOffset+1:(log(fineNumCells)/log(2) -1 )
    thisErrStruct = struct();
    
    thisErrStruct.numCells = 2^j;
    thisErrStruct.singleLevel = NaN;
    thisErrStruct.richardson = NaN;
    thisErrStruct.rate = NaN;
    thisErrStruct.nref2 = NaN;
    thisErrStruct.nref4 = NaN;
    thisErrStruct.nref22 = NaN;
    
    errStruct(j-numCellsOffset) = thisErrStruct;
end


thisRunTime = NaN*ones(length(runs), 10);
thisError = thisRunTime;
thisNumCells = thisRunTime;

for i=1:length(runs)
    folders = dir([base_dir, runs{i}, '*']);
    
    for j=1:length(folders)
        folder_dir = [base_dir, folders(j).name];
        fprintf('%s \n', folders(j).name);
        
        %errorFile = [folder_dir, '/errorAnalysis.mat'];
        %if exist(errorFile, 'file') ~= 2
        %    continue
        %end
        
        try
            
            % Need to replace this with parsing the calculations in
            % ChomboCompare/Nxcoarse-Nxfine/
            %load(errorFile);
            
            % Need to get:
            % runtime, err.(compname).(errtype)
            runTime = getRuntime(folder_dir);
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
            
            thisRun = runs{i};
            shortRun = strrep(thisRun, [folder, '-'], '');
            
            chComparePoutFile = [CompareDir, shortRun, num2str(numCells), '-', num2str(fineNumCells), '/pout.0'];
            comparePoutFine = ChComparePout(chComparePoutFile);
            errFine = comparePoutFine.getErr();
            
            chComparePoutFileRichardson = [CompareDir, 'Richardson', num2str(numCells), '-', num2str(2*numCells), '/pout.0'];
            comparePoutRichardson = ChComparePout(chComparePoutFileRichardson);
            errRichardson = comparePoutRichardson.getErr();
            
            if numCells*totalRefinement == fineNumCells/2 && maxLev<2 ...
                    && (length(strfind(folders(j).name, 'AMR') ) > 0 || length(strfind(folders(j).name, 'Uniform') ) > 0)
                
                pout = Pout([folder_dir, '/pout.0']);
                
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
            
            %numCells = strrep(folders(j).name, runs{i}, '');
            %numCells = strrep(numCells, '-ref2--0', '');
            %numCells = strrep(numCells, '--0', '');
            %numCells = str2num(numCells);
            
            
            
            thisRunTime(i, j) = runTime;
            if isfield(errFine, compName)
                thisError(i, j) = errFine.(compName).(errType);
            end
            thisNumCells(i, j) = numCells;
            
            exponent = log(numCells)/log(2) - numCellsOffset;
            thisErrStruct = errStruct(exponent);
            
            if strfind(folders(j).name, 'Uniform')
                thisErrStruct.singleLevel = thisError(i, j);
                
                % Also get richardson error
                                
                thisErrStruct.richardson = errRichardson.(compName).(errType);
                
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


% Now do plots
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
    plot(log(thisNumCells(i, :))/log(2), log(thisRunTime(i, :))/log(2), marker{i});
end
hold off;
box on;
xlabel('log$_2$($N_x$)'); ylabel('log$_2$(Runtime)');

subplot(m, n, 2)

hold on;
for i=1:length(runs)
    plot(log(thisNumCells(i, :))/log(2), log(thisError(i, :))/log(2), marker{i});
end
hold off;
box on;
xlabel('log$_2$($N_x$)'); ylabel('log$_2$(Error)');



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


% figure();
% hold on;
% for i=1:length(runs)
%     ncell = log10(thisNumCells(i, :));
%     runt = log10(thisRunTime(i, :));
%     validVals = 1-isnan(ncell).*isnan(runt);
%     
%     
%     extrapncell = sort([ncell(validVals==1) max(ncell) max(ncell)*1.5 max(ncell)*2 max(ncell)*3]);
%     
%     P = polyfit(sort(ncell(validVals==1)), sort(runt(validVals==1)), 2);
%     
%     extraprun = P(3) + extrapncell.^2*P(1) +  extrapncell*P(2);
%     
%     %extrapncell = [ncell(validVals==1) max(ncell) max(ncell)*3];
%     
%     %extraprun = interp1(ncell(validVals==1), runt(validVals==1), extrapncell, 'linear', 'extrap');
%     
%     plot(log10(thisNumCells(i, :)), log10(thisRunTime(i, :)), marker{i});
%     plot(extrapncell, extraprun, '-');
%     
% end
% hold off;
% box on;
% xlabel('num cells'); ylabel('runtime');


% Compute single level rates


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print tables, first latex then human readable format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fineLevelDiff = sprintf('Single-level %d difference', fineNumCells); 

% Error
fprintf('%5s & %15s  & %5s  & %13s & %12s  & %12s  & %12s \\\\ \n \hline \n ', ...
    '$1/\Delta z$', 'Richardson', ...
    'Rate', fineLevelDiff, '$n_\text{ref} = 2$', '$n_\text{ref} = 4$', '$n_\text{ref} = (2,2)$');

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
fprintf('%7s & %15s  & %15s & %20s & %15s  \\\\ \n \hline \n', '$n_\text{ref}$', 'Cells advanced', 'CPU time (s)',...
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
fprintf('%5s | %15s  | %5s  | %13s | %12s  | %12s  | %12s \n', '1/dz', 'Richardson', 'Rate', fineLevelDiff, 'n_{ref} = 2', 'n_{ref} = 4', 'n_{ref} = (2,2)');

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




temp=0;


end





