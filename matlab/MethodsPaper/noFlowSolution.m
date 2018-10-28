function noFlowSolution(dataFolder, figureName)

if nargin < 2
   %dataFolder = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/matlab/MethodsPaper/';
   dataFolder = getDataDir('Test/NoFlow/');
   figureName = 'benchmarkNoFlow.pdf';
end

close all;

set(0, 'defaultlinelinewidth',2);
set(0, 'defaultaxeslinewidth',2);
set(0, 'defaultpatchlinewidth',2);
set(0, 'defaultAxesFontSize',22);


savePlots = true;


uniformPrefix = 'Uniform-noFlow-';
gridRes = [8,16,32,64,128,256,512];
dz = 4./gridRes;
errUniform = NaN*ones(length(gridRes), 1);
errAMR = NaN*ones(length(gridRes), 1);

amrPrefix2lev = 'AMR-Subcycle-Reflux-Freestream0.95-MaxLevel2-ref2-noFlow-';
amrPrefix1lev = 'AMR-Subcycle-Reflux-Freestream0.95-MaxLevel1-ref2-noFlow-';


textFontSize = 16;
labelPos = [0.05 4.0];

plotWindowSize =  [200 200 1400 550];

plotHeight = 0.7;
plotWidth = 0.18;
plotBottom = 0.18;
horizSpacing = 0.03;

doColorbar = true;

fprintf('Loading 2 level solution \n');
output2lev        = getFinalPlotFile(fullfile(dataFolder, [amrPrefix2lev, '16--0']));

fprintf('Loading uniform fine solution \n');
outputUniformFine = getFinalPlotFile(fullfile(dataFolder, [uniformPrefix, '64--0']));
%output = MushyLayerOutput(2, 224, dataFolder, AMRplot_prefix, true);

perm = output2lev.dataForComp(output2lev.components.Porosity);
T = output2lev.dataForComp(output2lev.components.Temperature);
Terr = output2lev.dataForComp(output2lev.components.Terr);
streamfunction = output2lev.getStreamfunction(3000, 1);

Tuniform  = outputUniformFine.dataForComp(output2lev.components.Temperature);
chiUniform = outputUniformFine.dataForComp(output2lev.components.Porosity);

axPos(1,:) = [0.05 plotBottom 0.13 plotHeight];
axPos(2,:) = [axPos(1,1)+axPos(1,3)+horizSpacing plotBottom 0.17 plotHeight];
axPos(3,:) = [axPos(2,1)+axPos(2,3) plotBottom 0.1 plotHeight];
axPos(4,:) = [axPos(3,1)+axPos(3,3)+horizSpacing  plotBottom 0.15 plotHeight];
axPos(5,:) =  [axPos(4,1)+axPos(4,3)+horizSpacing*4 plotBottom plotWidth plotHeight];

%perm = perm(:, 40:(end-5));
%perm = log10(perm);

probDomain = output2lev.problemDomain;
dx = probDomain.dxCoarse;
numLevels = length(output2lev.levelArray);
dxFine = dx/2^(numLevels-1);

numx = 2^(numLevels-1)*((probDomain.domainExtent.hi_i - probDomain.domainExtent.lo_i) + 1);
numy = 2^(numLevels-1)*((probDomain.domainExtent.hi_j - probDomain.domainExtent.lo_j) + 1);
x = linspace(double(probDomain.domainExtent.lo_i), double(probDomain.domainExtent.hi_i)+1, numx)*dx;
y = linspace(double(probDomain.domainExtent.lo_j), double(probDomain.domainExtent.hi_j)+1, numy)*dx;

x = x - min(x);
x = x*1.02;
%y = y*1.02;

[X, Y] = meshgrid(x, y);

xlo = double(probDomain.domainExtent.lo_i)*dx;
xhi = double(probDomain.domainExtent.hi_i)*dx;
width = xhi-xlo;

h = figure();
set(h, 'Position', plotWindowSize);
%colormap bone;

m = 1;
n = 5;

subplot(m, n, 1);

axBCs = gca;
xlim([0 1]);
ylim([0 4]);
box on;

axBCs.YTick = [];
axBCs.XTick = [];
%axBCs.YTickLabels = {'-4','0'};
%axBCs.XTickLabels = {'0','1'};

%text(-0.2, 4.35, '$H=H(\theta_e=0, \chi_e)$,', 'FontSize', 16);
%text(0.2, 4.1, '$\Theta=-1$.', 'FontSize', 16);
text(0.0, 4.25,{'$H=H(\theta_e=0, \chi_e)$,', '$\Theta=-1$.'}, 'FontSize', 16);

%text(0.2, -0.2, '$\Theta=-1$,', 'FontSize', 16);
%text(-0.3, -0.42, '$H=H(\theta_{analytic}, \chi=1)$.', 'FontSize', 16);
text(0.0, -0.35, {'$H=H(\theta_{analytic}, \chi=1)$,','$\Theta=-1$.' }, 'FontSize', 16);

text(-0.1, 1.6, 'Periodic', 'FontSize', 16, 'Rotation', 90);
text(1.12, 1.6, 'Periodic', 'FontSize', 16, 'Rotation', 90);
text(0.52, 1.8, '$V$', 'FontSize', 16);
annotation('arrow',[0.1 0.1], [0.45 0.6]);

axBCs.Position = axPos(1,:);

text(-0.2,labelPos(2), '(a)', 'FontSize', textFontSize);

subplot(m, n, 2);

%colormap(makeColorMap([0.0039    0.0980    0.3216], [1 1 1]));
colormap(makeColorMap( [1 1 1], [0.0039    0.0980    0.3216]));
colormap(flipud(parula));


hcolor = pcolor(x, y, perm.');
set(hcolor,'edgecolor','none');

axLeft= gca;
set(axLeft, 'Layer', 'top'); % put box border above imagesc

axLeft.Position = axPos(2, :);

hold on
meshColor(1, :) = [1 1 1];
meshColor(2, :) = [0 0 0 ];
meshColor(3, :) = [0 0 0];

opacity(1) = 0.0;
opacity(2) = 0.0;
opacity(3) = 0.0;

edgeColor(1, :) = [1 1 1];
edgeColor(2, :) = [1 0 1];
edgeColor(3, :) = [0 1 1];

% Draw on the different meshes
for l = 2:length(output2lev.levelArray)
    lev = output2lev.levelArray(l);
    levDx = 0.5*lev.dx;
    lev_width = lev.xhi-lev.xlow+2*levDx;
    lev_height = lev.yhi-lev.ylow+2*levDx;
    lev_bottom = lev.ylow;
    fcl = meshColor(l, :);
    ecl = edgeColor(l, :);
    fcl(4) = opacity(l); % opacity
    rectangle('Position', [0.02 lev_bottom  0.99 3.98-lev_bottom],...
        'FaceColor', fcl, 'EdgeColor', ecl, 'LineWidth', 2.0);
    
    
end
hold off


if doColorbar
    c = colorbar('north');
    c.Label.String = '\chi';
    
    
    cbar_loc = 'top';
    if strcmp(cbar_loc, 'top')
        c.Label.Rotation = 0;
        
        thisaxPos = axLeft.Position;
        
        c.Position = [thisaxPos(1) + 0.04 thisaxPos(2) + thisaxPos(4) + 0.02 thisaxPos(3)*0.5 0.03];
        
        
        c.Label.Position = [0.9 1.6]; % this is relative to the colorbar position
        
    else
        c.Label.Rotation = 0;
        c.Label.Position = [2.2 0.6];
        c.Position = [0.8 0.22 0.06 0.5];
    end
    
    smallestPerm = min(min(perm));
    if smallestPerm < 0.1
        
        c.Ticks = [0.01 0.99];
        c.TickLabels = {'0', '1'};
    else
        if smallestPerm > 0.99
            smallestPerm = 0.0;
        end
        c.Ticks = [smallestPerm+0.01 0.99];
        c.TickLabels = {sprintf('%1.1f', smallestPerm), '1'};
        
    end
    
    c.TickLength = 0;
    
end

%hold off;



set(axLeft,'dataAspectRatio',[1 1 1])


%set(ax,'clim',[0 1]);

box on;

xlabel('$x$', 'Interpreter','latex');
ylabh = ylabel('$z$', 'Interpreter','latex');
set(ylabh, 'Rotation', 0);
%set(ylabh, 'Position', 0);


set(axLeft,'TickLength',[0 0])

axLeft.YTick = ([0 3.95]);
axLeft.YTickLabels = ({'-4', '0'});

axLeft.XTick = ([0 0.97]);
axLeft.XTickLabels = ({'0', '1'});

text(-0.7,labelPos(2), '(b)', 'FontSize', textFontSize); %'Color', [1 1 1]

% Now make the figure next to it
subplot(m, n, 3);


hold on;
% Save the colours for the final plot
plot(Tuniform(2, :), y, 'k-');
plot(chiUniform(2, :), y, 'k--');

box on;
hold off;
%xlabel('$\theta$');
axMiddle = gca;
axMiddle.Position = axPos(3,:);

xlim([0 1.2]);


axMiddle.YTick = [];

legend({'$\theta$','$\chi$'}, 'Location', 'southwest');

text(-0.3,labelPos(2), '(c)', 'FontSize', textFontSize);

% Now make the figure next to it
subplot(m, n, 4);



%semilogx(abs(Terr(2, :)), y);
plot(log10(abs(Terr(2, :))), y, '-');
hold on;
minErr = min(Terr(2, :));

gridResPlot = [16 32 64];
for f =1:length(gridResPlot)
    thisRes = gridResPlot(f);
    
    folder = [uniformPrefix, num2str(thisRes),'--0'];
    
    fprintf(['Loading 2 level solution with resolution ', num2str(thisRes), ' \n']);
    output2lev = getFinalPlotFile(fullfile(dataFolder, folder));
    
    %output = MushyLayerOutput(2, uniformframes(f), dataFolder, plot_prefix, true);
    if length(output2lev.levelArray) > 0
        
        Terr = output2lev.dataForComp(output2lev.components.Terr);
        
        probDomain = output2lev.problemDomain;
        dx = probDomain.dxCoarse;
        numLevels = length(output2lev.levelArray);
        dxFine = dx/2^(numLevels-1);
        
        numx = 2^(numLevels-1)*((probDomain.domainExtent.hi_i - probDomain.domainExtent.lo_i) + 1);
        numy = 2^(numLevels-1)*((probDomain.domainExtent.hi_j - probDomain.domainExtent.lo_j) + 1);
        x = linspace(double(probDomain.domainExtent.lo_i), double(probDomain.domainExtent.hi_i)+1, numx)*dxFine;
        y = linspace(double(probDomain.domainExtent.lo_j), double(probDomain.domainExtent.hi_j)+1, numy)*dxFine;
        
        %semilogx(abs(Terr(2, :)), y);
        Te = squeeze(Terr(2, :));
        plot(log10(abs(Te)), y, '--');
       
%         minErr = min(minErr, min(Te));
%         
%         errUniform(f) = mean(abs(Te));
        
    end
    
%     plot_prefix = [amrPrefix1lev, num2str(gridRes(f)),'--0'];
%     output = getFinalPlotFile([dataFolder, plot_prefix]);
%     if length(output.levelArray) > 0
%         TerrAMR = output.dataForComp(output.components.Terr);
%         Te = squeeze(TerrAMR(2,:));
%         errAMR(f) = mean(abs(Te));
%     end
    
end
box on;
hold off;
xlabel('log$_{10} (|\theta_{err}|)$');


forceRecalculate = false;
err2lev = getErr(dataFolder, amrPrefix2lev, gridRes, forceRecalculate);
errAMR = getErr(dataFolder, amrPrefix1lev, gridRes, forceRecalculate);
errUniform = getErr(dataFolder, uniformPrefix, gridRes);

timeUniform = getTimes(dataFolder, uniformPrefix, gridRes);
timeAMR = getTimes(dataFolder, amrPrefix1lev, gridRes);
time2lev = getTimes(dataFolder, amrPrefix2lev, gridRes);

%Print out times
fprintf('Nx,t_uniform,t_amr,t_amr2lev \n');
for i=1:length(gridRes)
   fprintf('%d,%1.5f,%1.5f,%1.5f \n', gridRes(i), timeUniform(i), timeAMR(i),time2lev(i));
    
end

err2lev(end-1:end) = NaN;
errAMR(end) = NaN;

%xlim([-1.4e-4 2e-3]);

axRight = gca;
axRight.Position = axPos(4,:);

axRight.YTick = [];
%axRight.XTick = [0  1e-4, 1e-3];
%axRight.XTickLabels = {'0', '10^{-4}', '10^{-3}'};

%axRight.XTick = [0  1e-4, 1e-3];
%axRight.XTickLabels = {'0', '10^{-4}', '10^{-3}'};

% legend({'16-32-64', '16', '32', '64'}, 'Location', 'northwest');
axRightPos = axRight.Position;
legend({'(16,32,64)', '16', '32', '64'}, 'Position',...
    [axRightPos(1)-0.08 axRightPos(2)+axRightPos(4)+0.05 0.2 0.06], ...
    'Orientation', 'horizontal');

text(-11.5,labelPos(2),'(d)','FontSize', textFontSize);

subplot(m, n, 5);
%Error plot
ax = gca;
box on;

plotUniform = plot(log10(dz), log10(errUniform), 'x-');
hold on;
plotAMR = plot(log10(dz), log10(errAMR), 'x-');
plotAMR2lev = plot(log10(dz), log10(err2lev), 'x-');

plotComparison = plot([-2.1 -0.1], [-5 -1], '--');

hold off;

xlim([-4 0]);

ax.Position = axPos(5,:);
xlabel('log$_{10}(\Delta z)$');
ylabel('log$_{10}(\bar{|\theta_{err}|})$');
legend([plotUniform, plotAMR, plotAMR2lev, plotComparison], {'Uniform', '$n_{ref}=2$', ...
    '$n_{ref}=(2,2)$', '2nd order'}, 'Location', 'northwest');

text(-5.2,-1,'(e)','FontSize', textFontSize);


h =gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

if savePlots
    %print(h,[dataFolder, 'noFlowSolution.png'],'-dpng','-r300')
    fprintf('Saved to %s \n', figureName);
    print(h,figureName,'-dpdf','-r50')
    
end

end



function err = getErr(dataFolder, prefix, gridRes, forceRecalculate)
if nargin < 4
    forceRecalculate  = false;
end

err = NaN*ones(length(gridRes), 1);

for f =1:length(gridRes)
    
    thisRes = gridRes(f);
    
    folder_name = [prefix, num2str(gridRes(f)),'--0'];
    thisFolder = fullfile(dataFolder, folder_name);
    errFile = [thisFolder, '/err.mat'];
    
    if exist(errFile, 'file') == 2 && ~forceRecalculate
        load(errFile)
    else
   
        fprintf('Loading final plot file from %s \n', thisFolder);
        output = getFinalPlotFile(thisFolder);

        if length(output.levelArray) > 0

            Terr = output.dataForComp(output.components.Terr);
            Te = squeeze(Terr(2,:));
            e.meanTerr =  mean(abs(Te));
            save(errFile, 'e');

        end


    end
    
    err(f) = e.meanTerr;
    
end


end


function times = getTimes(dataFolder, prefix, gridRes)

times = NaN*ones(length(gridRes), 1);

for f =1:length(gridRes)
    
thisRes = gridRes(f);
    
    folder_name = [prefix, num2str(thisRes),'--0'];
    thisFolder = fullfile(dataFolder, folder_name);
    
    times(f)  = getRuntime( thisFolder );
end

end

