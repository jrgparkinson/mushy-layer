close all;

getData = true;
if getData
    clear all;
    getData = true;
end

dim = 2; subcycled = true;

printQuality = '-r100'; % -r0 produces very large files

exactNu = 0.0;

% short_run_names = {'Uniform', ...
%     'AMRNoSubcycleNoReflux', 'AMRNoSubcycleReflux', ...
%     'AMRSubcycleNoReflux', 'AMRSubcycleReflux', ...
%     'AMRNoSubcycleRefluxFreestream', 'AMRSubcycleRefluxFreestream'};

% short_run_names = {'Uniform', ...
%     'AMRNoSubcycleNoReflux', 'AMRNoSubcycleReflux', ...
%     'AMRNoSubcycleRefluxFreestream0.5', ...
%     'AMRNoSubcycleRefluxFreestream0.53', ...
%     'AMRNoSubcycleRefluxFreestreamNoInit0.5'};

% errorPlotFormat = {'b-', 'red:', 'red--', 'k:', 'k--', 'g--', 'g-.'};


% short_run_names = {
%     'AMRNoSubcycleNoReflux-HRLOptShort', 'AMRNoSubcycleRefluxFreestream0.5-HRLOptShort', 'AMRNoSubcycleRefluxFreestreamNoInit0.5-HRLOptShort', ...
%     'AMRNoSubcycleNoReflux-HRLOptLong',   'AMRNoSubcycleRefluxFreestream0.5-HRLOptLong',  'AMRNoSubcycleRefluxFreestreamNoInit0.5-HRLOptLong'};
%
% errorPlotFormat = {'b--o', 'b:x', 'b-.d', 'red--o', 'red:x', 'red-.d'};

% short_run_names = {
%     'AMRNoSubcycleNoReflux-HRLOptShort', 'AMRNoSubcycleRefluxFreestream0.5-HRLOptShort', 'Uniform-HRLOptShort', ...
%     'AMRNoSubcycleNoReflux-HRLOptLong',   'AMRNoSubcycleRefluxFreestream0.5-HRLOptLong',  'Uniform-HRLOptLong'};
%
% errorPlotFormat = {'b--o', 'b:x', 'b-d', 'red--o', 'red:x', 'red-d'};
%

% HRL 
% short_run_names = {
%     'Uniform-HRLRa100', ...
%     'AMRSubcycleRefluxFreestream0.5-HRLRa100'};
% 
% errorPlotFormat = {'r-x', 'b-o'};
% baseDir = '/home/parkinsonjl/mnt/sharedStorage/AW002_PARKINSON_MUSH/AMRConvergenceTestHRL/';
% coarseNzs = [4, 8, 16];
% plotLabels = {'Uniform mesh', 'Adaptive mesh'};
% exactSolRes = '32';
% exactML = ChomboOutput(dim, 10000, [baseDir, 'Uniform-HRLRa100-',exactSolRes,'-0/'],...
%     ['AMRConvergenceTestHRL-Uniform-HRLRa100-',exactSolRes,'-']);
% exactSol = ChomboCompare(exactML);


%%%% Diffusive solidification
short_run_names = {
    'Uniform-diffusiveSolidification', ...  % previously had NoInit at the end
      'AMRSubcycleRefluxFreestream0.5Levels2-diffusiveSolidification', ... %'AMRSubcycleRefluxFreestream0.5-diffusiveSolidificationNoInit', ...
      'FixedSubcycleRefluxFreestream0.5-diffusiveSolidification', ...
      'AMRSubcycleRefluxFreestream0.5Levels3-diffusiveSolidification'
     };
 
short_run_names = {'Uniform-diffusiveSolidification'};
 
plotLabels = {'Uniform mesh', 'Adaptive Mesh (2 levels)', 'Variable Mesh (3 levels)', 'Adaptive Mesh (3 levels)'};
plotLabels = {'Uniform mesh'};
errorPlotFormat = {'r-x', 'b-o', 'm-d', 'k--x'};

coarseNzs = [8, 16, 32, 64, 128];
baseDir = '/home/parkinsonjl/mnt/sharedStorage/AW002_PARKINSON_MUSH/AMRConvergenceTestNoFlow/';
baseDir = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/run/AMRConvergenceTestNoFlowV2/';
plotPrefix = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/matlab/Tests/diffusiveSolidification';



% Convection with Darcy Brinkman
% short_run_names = {
%     'Uniform-convectionDarcyBrinkman', ...
%     'VariableMeshSubcycleRefluxFreestream0.45-convectionDarcyBrinkman', ...
%     'VariableMesh2SubcycleRefluxFreestream0.45-convectionDarcyBrinkman', ...
%     'VM3LevelsSubcycleRefluxFreestream0.45-convectionDarcyBrinkman'
%      };
%  
%   % 'AMRSubcycleRefluxFreestream0.5-convectionDarcyBrinkman', ...
%  
% plotLabels = {'Uniform mesh', '1/4 VM (2 levels)', '1/2 VM (2 levels)', '1/2 VM (3 levels)'};
% errorPlotFormat = {'r-x', 'b-o', 'm-d', 'k--x'};
% 
% coarseNzs = [8, 16, 32];
% baseDir = '/home/parkinsonjl/mnt/sharedStorage/AW002_PARKINSON_MUSH/AMRConvergenceTestConvectionDB/';
% plotPrefix = 'convectionDB';

% exactSolRes = '128';
% exactML = ChomboOutput(dim, 4600, [baseDir, 'Uniform-convectionDarcyBrinkman-',exactSolRes,'-0/'],...
%      ['AMRConvergenceTestConvectionDB-Uniform-convectionDarcyBrinkman-',exactSolRes,'-']);
% exactSol = ChomboCompare(exactML);
% exactDiags = getDiagnostics([baseDir, 'Uniform-convectionDarcyBrinkman-',exactSolRes,'-0/diagnostics.out']);
% exactNu = exactDiags.Nusselt(end);


% Darcy no subcycle
% short_run_names = {'Uniform-MushyLayerDarcy', ...
%     'AMRNoSubcycleNoReflux-MushyLayerDarcy', ...
%     'AMRNoSubcycleReflux-MushyLayerDarcy', ...
%     'AMRNoSubcycleRefluxFreestream0.5-MushyLayerDarcy'};
% errorPlotFormat = {'r-x', 'b--x', 'g:o', 'k--+'};

% Darcy subcycle
% short_run_names = {'Uniform-MushyLayerDarcy', ...
%     'AMRSubcycleNoReflux-MushyLayerDarcy', ...
%     'AMRSubcycleReflux-MushyLayerDarcy', ...
%     'AMRSubcycleRefluxFreestream0.5-MushyLayerDarcy', ...
%     'AMRSubcycleRefluxFreestream0.49-MushyLayerDarcy', ...
%     'AMRSubcycleRefluxFreestream0.51-MushyLayerDarcy'};
% errorPlotFormat = {'r-x', 'g:x', 'b--o', 'k-.+', 'c-.+', 'm-.+'};
% exactML = ChomboOutput(dim, 1000, [baseDir, 'Uniform-MushyLayerDarcy-',exactSolRes,'-0/'],...
%     ['AMRConvergenceTestMushyLayer-Uniform-MushyLayerDarcy-',exactSolRes,'-']);
% exactSol = ChomboCompare(exactML);

%short_run_names = {'Uniform'};
%'AMRNoSubcycleLinearReflux', 'AMRSubcycleLinearReflux', ...


%physicalProblem = 'diffusiveSolidificationNoOpt';
%physicalProblem = 'HRLOpt';

%

%errorPlotFormat = {'b-', 'red:', 'red--', 'k:', 'k--', 'g--', 'g-.'};
%errorPlotFormat = {'b-', 'red:', 'red--', 'k:', 'k--', 'g--', 'm--', 'c-.'};

runPlotFormat = errorPlotFormat; %{'bx-', 'redx:', 'redx--', 'kx:', 'kx--'};
%for i=1:length(errorPlotFormat)
%    runPlotFormat{end+1} = [errorPlotFormat{i}, 'x'];
%end



if getData
    
    err = NaN*ones(length(short_run_names), length(coarseNzs));
    nu_err = NaN*ones(length(short_run_names), length(coarseNzs));
    pointsUpdated = err;
    fluxConservation = err;
    
    % Make a figure to plot all the figures
    
    if (findstr(short_run_names{1}, 'diffusiveSolidification'))
        hErrPlots = figure('PaperPositionMode', 'auto');
        set(hErrPlots, 'Position', [100 100 1600 900]);
        subplot_i = 1;
        
        m = 1;
        n = length(coarseNzs)+1;

    end
    
    
    for nz_i = 1:length(coarseNzs)
        
        coarseNz = coarseNzs(nz_i);
        
        if (findstr(short_run_names{1}, 'diffusiveSolidification'))
            subplot(m, n, subplot_i);
            subplot_i = subplot_i + 1;
            hold on;
        else
            h  = figure('PaperPositionMode', 'auto');  set(gcf,'Visible', 'off');
            set(h, 'Position', [200 200 1600 900]);
            
            subplot_m = 3; subplot_n = 2;
            subplot_i = 1;
        end
        
        
        
        for r = 1:length(short_run_names)
            thisRun = [short_run_names{r}, '-', num2str(coarseNz)];
            
            output_dir = [baseDir, thisRun, '-0/'];
            
%             if (findstr(thisRun, 'MushyLayer'))
%                 plot_prefix = ['AMRConvergenceTestMushyLayer-', thisRun , '-'];
%             else
%                 
%                 plot_prefix = ['refluxTest-', thisRun , '-'];
%             end
            
           
            
            search = [output_dir, '*.2d.hdf5'];
            all_files = dir(search);
            files = [];
            num_files = length(all_files);
            f_i = 1;
            
            while f_i <= num_files
                name = all_files(f_i).name;
                if sum(name(1:3) == 'chk') == 3 % three, because 4 characters
                    all_files(f_i) = [];
                    num_files =num_files -1;
                else
                    f_i = f_i + 1;
                end
                
            end
            
            files = all_files;
            
            diags = getDiagnostics([output_dir, 'diagnostics.out']);
            
            pout = Pout([output_dir, 'pout.0']);
            
            
            
            if pout.pointsUpdated
                pointsUpdated(r, nz_i) = pout.pointsUpdated;
            end
            
            %             if pout.steadyState
            %                 if isfield(diags, 'Fh_abs_mismatch')
            %                     fluxConservation(r, nz_i) = diags.Fh_abs_mismatch(end);
            %                 else
            %                 fluxConservation(r, nz_i) = pout.finalHeatFluxMismatch;
            %                 end
            %             end
            
            if isfield(diags, 'Fh_abs_mismatch')
                fluxConservation(r, nz_i) = abs(diags.Fh_abs_mismatch(end));
            else
                fluxConservation(r, nz_i) = NaN;
            end
            
            frames = [];
            
            actual_plot_prefix = '';
            for file_i = 1:length(files)
                fname = strrep(files(file_i).name, '.2d.hdf5', '');
                parts = strsplit(fname, '-');
                frames(end+1) = str2num(parts{end});
                
                actual_plot_prefix = strrep(fname, parts{end}, '');
            end
            
                       
            frame = max(frames);
            %frame = 0; % testing
            
            thisML = MushyLayerOutput(dim, frame, output_dir, actual_plot_prefix, subcycled);
            
            if length(thisML.levelArray()) > 0
                
                if (findstr(short_run_names{1}, 'diffusiveSolidification'))
                    
                    %T_err = thisML.dataForComp(thisML.components.Terr);
                    %thisErr = max(max(abs(T_err)));
                    [L1, L2, Max, Sum] = AMRSum(thisML, thisML.components.Terr);
                    
                else
                    %                     try
                    if (findstr(thisRun, 'Long'))
                        MLdiff =  exactSolLong.diff(thisML, [thisML.components.Temperature], ...
                            [thisML.components.Temperature]);
                    elseif (findstr(thisRun, 'Short'))
                        MLdiff =  exactSolShort.diff(thisML, [thisML.components.Temperature], ...
                            [thisML.components.Temperature]);
                    else
                        MLdiff =  exactSol.diff(thisML, [thisML.components.Temperature], ...
                            [thisML.components.Temperature]);
                    end
                    
                    [L1, L2, Max, Sum] = AMRSum(MLdiff, thisML.components.Temperature);
                    
                    %                     catch
                    %                         warning('Error whilst trying to compute difference with analytic sol');
                    %                        L1 = NaN;
                    %                        L2 = NaN;
                    %                        Max = NaN;
                    %                        Sum = NaN;
                    %                     end
                    % L2 = NaN;
                end
                
                % Only record error if we've reached steady state
                if pout.steadyState
                    err(r, nz_i) = L2;
                    nu_err(r, nz_i) = abs(diags.Nusselt(end) - exactNu)/exactNu;
                end
                
                
                % Make plot
                
                [X, Y] = thisML.grid();  y = Y(:, 1); x = X(1, :);
                
                if (findstr(short_run_names{1}, 'diffusiveSolidification'))
                    dat = thisML.dataForComp(thisML.components.Terr);
                    errField = dat(1, :);
                    
                    formatStr = errorPlotFormat{r};
                    
                    plot(y, errField, formatStr);
                    
                    %set(gca, 'XTickLabel', []);
                    box on;
                    
                    % blank second line to avoid overlap with scientific notation
                    title({['Coarse $dz$ = ', num2str(coarseNz)], ''});
                    xlabel('$z$'); ylabel('Error');
                    
                    if nz_i == 3
                        legend(plotLabels, 'Location', 'southoutside');
                    end
                    %axis equal;
                    
                    
                else
                    
                    dat = MLdiff.dataForComp(MLdiff.components.Temperature);
                    %errField = dat(round(end/2), :);
                    
                    
                    
                    subplot(subplot_m, subplot_n, subplot_i);
                    H = pcolor(x, y, dat.');
                    set(H,'edgecolor','none');
                    caxis([-Max Max]);
                    colormap(bluewhitered);
                    
                    colorbar();
                    title(short_run_names{r});
                    subplot_i = subplot_i + 1;
                end
                
                
                
                
            end
            
            
            
        end
        
        hold off;
        
        set(gcf,'Units','Inches');
        pos = get(gcf,'Position');
        set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        
        print(gcf, '-dpdf', [plotPrefix, num2str(coarseNz), '.pdf'], printQuality);
        
    end
    
    
end % end if get data


% Second order error. Start from uniform coarsest error and then scale with
% dx^2
secondOrderErr = err(1, :);
for i=2:length(secondOrderErr)
    secondOrderErr(i) = ((coarseNzs(i)/coarseNzs(i-1))^(-2))*secondOrderErr(i-1);
end


% Plot error
h  = figure('PaperPositionMode', 'auto'); 
set(h, 'Position', [200 200 1600 800]);


m = 2;
n = 3;
subplot(m, n, 1);

hold on;
for r = 1:length(short_run_names)
    plot(log10(1./coarseNzs), log10(err(r, :)),  runPlotFormat{r});
end

plot(log10(1./coarseNzs), log10(secondOrderErr),':');
hold off;

%legend(plotLabels, 'Location', 'eastoutside');

box on; grid on;
xlabel('$\log_{10}(dz_{coarse})$');
ylabel('$\log_{10}$(error)');



subplot(m, n, [2, 3]);

hold on;
for r = 1:length(short_run_names)
    plot(log10(1./coarseNzs), log10(nu_err(r, :)),  runPlotFormat{r});
end

%plot(log10(1./coarseNzs), log10(secondOrderErr),':');
hold off;

legend(plotLabels, 'Location', 'eastoutside');

box on; grid on;
xlabel('$\log_{10}(dz_{coarse})$');
ylabel('$\log_{10}$(Nu err)');



subplot(m, n, 4);

hold on;
for r = 1:length(short_run_names)
    plot(log10(pointsUpdated(r, :)), log10(err(r, :)), runPlotFormat{r});
end
hold off;

box on; grid on;
xlabel('$\log_{10}$(grid points updated)');
ylabel('$\log_{10}$(error)');



subplot(m, n, 5);

hold on;
for r = 1:length(short_run_names)
    plot(log10(1./coarseNzs), log10(fluxConservation(r, :)), runPlotFormat{r});
end
hold off;

box on; grid on;
xlabel('$\log_{10}(dz_{coarse})$');
ylabel('$\log_{10}$(flux mismatch)');



set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf, '-dpdf', [plotPrefix, 'Summary.pdf'], printQuality);

