close all;

set(0, 'defaultlinelinewidth',2);
set(0, 'defaultaxeslinewidth',2);
set(0, 'defaultpatchlinewidth',2);
set(0, 'defaultAxesFontSize',22);


getData = true;
if getData
    clear all;
    getData = true;
end

dim = 2; subcycled = true;

printQuality = '-r100'; % -r0 produces very large files

exactNu = 0.0;

plotWidth = 0.6;
plotHeight = 0.38;
fluxPos = [0.15       0.12                        plotWidth, plotHeight];
errPos  = [fluxPos(1) fluxPos(2)+plotHeight+0.05 plotWidth plotHeight];


%%%% Diffusive solidification
 
short_run_names = {'Uniform-diffusiveSolidification', ...
    'AMRSubcycleRefluxFreestream0.45-diffusiveSolidification',  ...
    'VM3LevelsSubcycleRefluxFreestream0.45-diffusiveSolidification'};

 
%plotLabels = {'Uniform mesh', 'Adaptive Mesh (2 levels)', 'Variable Mesh (3 levels)', 'Adaptive Mesh (3 levels)'};
plotLabels = {'Uniform mesh', '$l_{max} = 1$', 'VM3'};
errorPlotFormat = {'r-x', 'b-o', 'm-d', 'k--x', 'c-x'};
errVar = '\theta';

coarseNzs = [8, 16, 32, 64, 128];
coarseNzs = [8, 16, 32, 64, 128];


%baseDir = '/home/parkinsonjl/mnt/sharedStorage/AMRConvergenceTestNoFlow/';
%baseDir = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/run/AMRConvergenceTestNoFlowV2/';
baseDir = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/run/AMRConvergenceTestNoFlowV3/';
plotPrefix = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/matlab/MethodsPaper/noFlowConvergence';

errplotLabels = {};

for i=1:length(coarseNzs)
   errplotLabels{end+1} = [num2str(coarseNzs(i))]; 
       
   
   errplotLabels{end+1} = ['(', num2str(coarseNzs(i)), ',', num2str(coarseNzs(i)*2), ')'];
   
   errplotLabels{end+1} = ['(', num2str(coarseNzs(i)), ',', num2str(coarseNzs(i)*2),',', num2str(coarseNzs(i)*4), ')'];
   %errplotLabels{end+1} = ['(', num2str(coarseNzs(i)), ',', num2str(coarseNzs(i)*2),',', num2str(coarseNzs(i)*4), ')'];
   %errplotLabels{end+1} = ['(', num2str(coarseNzs(i)), ',', num2str(coarseNzs(i)*2),',', num2str(coarseNzs(i)*4), ')'];
end


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
    
    L2err = NaN*ones(length(short_run_names), length(coarseNzs));
    L1err = L2err;
    Maxerr = L2err;
    nu_err = NaN*ones(length(short_run_names), length(coarseNzs));
    pointsUpdated = L2err;
    fluxConservation = L2err;
    
    % Make a figure to plot all the figures
    
    if (findstr(short_run_names{1}, 'diffusiveSolidification'))
        hErrPlots = figure('PaperPositionMode', 'auto');
        set(hErrPlots, 'Position', [100 100 1600 900]);
        %set(hErrPlots, 'visible', 'false');
      %  subplot_i = 1;
        
     %   m = 1;
     %   n = length(coarseNzs)+1;

    end
    
    
    for nz_i = 1:length(coarseNzs)
        
        coarseNz = coarseNzs(nz_i);
        
        if (findstr(short_run_names{1}, 'diffusiveSolidification'))
            %subplot(m, n, subplot_i);
            %subplot_i = subplot_i + 1;
            hold on;
        else
            h  = figure('PaperPositionMode', 'auto');  set(gcf,'Visible', 'off');
            set(h, 'Position', [200 200 1600 900]);
            
            subplot_m = 3; subplot_n = 2;
            subplot_i = 1;
        end
        
        
        
        for r = 1:length(short_run_names)
            thisRun = [short_run_names{r}, '-Nz', num2str(coarseNz)];
            
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
                    L1err(r, nz_i) = L1;
                    L2err(r, nz_i) = L2;
                    maxerr(r, nz_i) = Max;
                    nu_err(r, nz_i) = abs(diags.Nusselt(end) - exactNu)/exactNu;
                end
                
                
                % Make plot
                
                [X, Y] = thisML.grid();  y = Y(:, 1); x = X(1, :);
                
                if (findstr(short_run_names{1}, 'diffusiveSolidification'))
                    dat = thisML.dataForComp(thisML.components.Terr);
                    errField = dat(1, :);
                    
                  %  formatStr = errorPlotFormat{r};
                    
                    logErr = log10(abs(errField));
                    
                    
                    errplotFormat = {'-', '--', ':'};
                    
                    
                    plot(y, logErr, errplotFormat{length(thisML.levelArray())});
                    
                    %set(gca, 'XTickLabel', []);
                    box on;
                    
                    % blank second line to avoid overlap with scientific notation
                    %title({['Coarse $dz$ = ', num2str(coarseNz)], ''});
                    
                    ylabel('Error');
                    
                    
                    xlabel('$z$'); 
                    
                    %if nz_i == 3
                    %   
                    %end
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
        
    end % loop over Nzs
    
    if (findstr(short_run_names{1}, 'diffusiveSolidification'))
         legend(errplotLabels, 'Location', 'eastoutside');
         end
    
    
    
end % end if get data


% Second order error. Start from uniform coarsest error and then scale with
% dx^2
% secondOrderErr = err(1, :);
% for i=2:length(secondOrderErr)
%     secondOrderErr(i) = ((coarseNzs(i)/coarseNzs(i-1))^(-2))*secondOrderErr(i-1);
% end

% Start from finest error with an offset
secondOrderErr = L2err(1, :)*2.0;
for i=length(secondOrderErr)-1:-1:1
    secondOrderErr(i) = ((coarseNzs(i)/coarseNzs(i+1))^(-2))*secondOrderErr(i+1);
end



% Plot error
h  = figure('PaperPositionMode', 'auto'); 
set(h, 'Position', [200 200 900 800]);


m = 2;
n = 1;
subplot(m, n, 1);

hold on;
for r = 1:length(short_run_names)
    plot(log10(1./coarseNzs), log10(L2err(r, :)),  runPlotFormat{r});
end

plot(log10(1./coarseNzs), log10(secondOrderErr),':');
hold off;

plotLabels{end+1} = '2nd order';
legend(plotLabels, 'Location', 'eastoutside');

box on; grid on;
%xlabel('$\log_{10}(dz_{coarse})$');
ylabel(['$\log_{10}(',errVar,'_{err})$']);


axErr = gca;
%axErr.YTick = [-4,-3,-2];
%axErr.YTickLabels = {'-4','','-2'};
%ylim([-5 -2]);
axErr.XTick = [-2.5 -2 -1.5 -1];
axErr.XTickLabels = {'','','', ''};

axErr.Position = errPos;

subplot(m, n, 2);

hold on;
for r = 1:length(short_run_names)
    plot(log10(1./coarseNzs), log10(fluxConservation(r, :)), runPlotFormat{r});
end
hold off;

box on; grid on;
xlabel('$\log_{10}(dz_{coarse})$');
ylabel('$\log_{10}(\Delta \mathbf{F}_\theta)$');

axFluxCons = gca;
%axFluxCons.YTick = [-12,-11,-10];
%axFluxCons.YTickLabels = {'-12','','-10'};

%ylim([-12 -10]);

axFluxCons.Position = fluxPos;

% subplot(m, n, [2, 3]);
% 
% hold on;
% for r = 1:length(short_run_names)
%     plot(log10(1./coarseNzs), log10(nu_err(r, :)),  runPlotFormat{r});
% end
% 
% %plot(log10(1./coarseNzs), log10(secondOrderErr),':');
% hold off;
% 
% legend(plotLabels, 'Location', 'eastoutside');
% 
% box on; grid on;
% xlabel('$\log_{10}(dz_{coarse})$');
% ylabel('$\log_{10}$(Nu err)');

errFormat = '%12.4e';
rateFormat = '%12.4f';

fprintf('%18s&%12s&%12s&%12s&%12s \n', '$1/\Delta z^{0}$', 'Single level error', 'Rate', '$n_{ref}=2$','$n_{ref}=(2,2)$');
for i=1:length(coarseNzs)
    
    errorMetric = L2err;
    
    singleLevelErr =  errorMetric(1, i);
    
    se = size(errorMetric);
    if se(1) > 1
        TwolevErr = errorMetric(2, i);
    else
        TwolevErr = NaN;
    end
    
    if se(1) > 2
    ThreelevErr = errorMetric(3, i);
    else
        ThreelevErr= NaN;
    end
    
    
    if i == 1
        % No rate
        fprintf(['%18d&',errFormat,'&%12s&',errFormat,'&',errFormat,'\n'],...
            coarseNzs(i), singleLevelErr, '-', TwolevErr, ThreelevErr);
    else
        
        rate = (errorMetric(1,i-1)/errorMetric(1,i)) * (coarseNzs(i-1)/coarseNzs(i));
         
    fprintf(['%18d&',errFormat,'&',rateFormat,'&',errFormat,'&',errFormat,'\n'],...
        coarseNzs(i), singleLevelErr, rate, TwolevErr, ThreelevErr);
    end
end
%T = table(coarseNzs,err(1, :),'VariableNames',{'Coarse_dx','Single_level_err'})


% subplot(m, n, 4);
% 
% hold on;
% for r = 1:length(short_run_names)
%     plot(log10(pointsUpdated(r, :)), log10(err(r, :)), runPlotFormat{r});
% end
% hold off;
% 
% box on; grid on;
% xlabel('$\log_{10}$(grid points updated)');
% ylabel('$\log_{10}$(error)');
% 
% 

% Let's also produce a table of the results which we can fairly easily
% Copy into latex



set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf, '-dpdf', [plotPrefix, 'Summary.pdf'], printQuality);

