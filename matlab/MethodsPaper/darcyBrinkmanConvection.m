% Plot a typical AMR solution and the error convergence
function darcyBrinkmanConvection

clear all; close all;

gridRes = [8 16 32 64];
dataFolder = getDataDir('AMRConvergenceTestConvectionDB-chi0.4-Da0.01-Ra1e5/');
figureName = [dataFolder, 'benchmark2Convergence'];
saveFigure = true;
%plotPrefixUniform = 'Uniform-convectionDarcyBrinkman-';
%plotPrefixVariable = 'VariableMeshSubcycleRefluxFreestream0.45-convectionDarcyBrinkman-';
%plotPrefixAMR = @N

dz = 1./gridRes;

highRes = getFinalPlotFile([dataFolder, plotPrefixUniform(gridRes(end)*2)]);
highResML = ChomboCompare(highRes);

errUniform = getErr(dataFolder, @plotPrefixUniform, gridRes, highResML);
errAMR1lev = getErr(dataFolder, @plotPrefixAMR, gridRes, highResML);


h = figure();
h.Position = [200 200 1200 500];
axHeight = 0.55;
axBottom = 0.25;

m = 1; n = 3;
subplot(m, n, 1)
axBCs = gca;
axBCs.Position = [0.05 axBottom 0.25 axHeight];


xlim([0 1]);
ylim([0 1]);
box on;

axBCs.YTick = [];
axBCs.XTick = [];
daspect([1 1 1]);



%text(0.0, 1.13 ,{'$H=H(\theta=1, \chi=1)$,', '$\Theta=0, \mathbf{U} = 0.$'}, 'FontSize', 16);
%text(0.0, -0.15, {'$H=H(\theta=0, \chi=1)$,','$\Theta=0, \mathbf{U} = 0.$' }, 'FontSize', 16);


text(0.1, 0.5, {'Eq''s 3-6 with:', '$\chi=0.4, St=0, Pr=1$.' }, 'FontSize', 16);

insulatingStr = '$\partial H/\partial x = \partial \Theta / \partial x = 0, \mathbf{U} = 0.$';
%text(-0.1, 0.0, insulatingStr, 'FontSize', 16, 'Rotation', 90);
%text(1.12 , 0.0, insulatingStr, 'FontSize', 16, 'Rotation', 90);


text(-0.1, 0.0, '$H=H(\theta=1, \chi=1), \Theta=0, \mathbf{U} = 0.$', 'FontSize', 16, 'Rotation', 90);
text(1.12 , 0.0, '$H=H(\theta=0, \chi=1), \Theta=0, \mathbf{U} = 0.$', 'FontSize', 16, 'Rotation', 90);

text(0.0, 1.13 ,insulatingStr, 'FontSize', 16);
text(0.0, -0.15, insulatingStr, 'FontSize', 16);

%axis arrows
annotation('arrow',[0.03 0.07], [0.06 0.06]);
annotation('arrow',[0.03 0.03], [0.06 0.16]);
text(-0.05, -0.4, '$x$', 'FontSize', 16);
text(-0.23 , -0.2   , '$z$', 'FontSize', 16);


text(0.04, 0.93, '(a)', 'FontSize', 16);

subplot(m, n, 2)

%AMRsol = getFinalPlotFile([dataFolder, plotPrefixUniform(32)]);
AMRsol = getFinalPlotFile([dataFolder, plotPrefixAMR(32)]);
UniformFineSol = getFinalPlotFile([dataFolder, plotPrefixUniform(64)]);

T = AMRsol.dataForComp(AMRsol.components.Temperature).';
psi = AMRsol.dataForComp(AMRsol.components.streamfunction).';
psiUniform = UniformFineSol.dataForComp(UniformFineSol.components.streamfunction).';

[X,Y] = AMRsol.grid();
[XUniform,YUniform] = UniformFineSol.grid();
%Create smooth psi



% Temperature
axTemperature = gca;

plot_T = T*2-1;
pcolor(X,Y,T*2-1); daspect([1 1 1]);
colormap(bluewhitered);
c = colorbar('Location', 'northoutside');
c.Ticks = [min(min(plot_T))*0.99 max(max(plot_T))*0.99];
c.TickLabels = {'0', '1'};
c.Label.String = '\theta';
oldPos = c.Label.Position;
c.Label.Position = [oldPos(1) oldPos(2)+1.2];
thetaColorbarPos = c.Position;

box on;
set(axTemperature, 'Layer', 'top');

% Psi
% Plot psi uniform because it's smoother
axPsi = axes;
contour(XUniform,YUniform,psiUniform); daspect([1 1 1]);
colormap(axPsi,flipud(parula));
cPsi = colorbar(axPsi,'Location', 'southoutside');
oldPos = cPsi.Position;

linkaxes([axTemperature,axPsi]);
axPsi.Position = [0.35 axBottom 0.3 axHeight];
axTemperature.Position = [0.35 axBottom 0.3 axHeight];



% Plot meshes
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
for l = 2:length(AMRsol.levelArray)
    lev = AMRsol.levelArray(l);
    levDx = 0.5*lev.dx;
    
    for b = 1:length(lev.boxArray)
        thisBox = lev.boxArray(b);
    
        %lev_width = lev.xhi-lev.xlow+2*levDx;
        %lev_height = lev.yhi-lev.ylow+2*levDx;
        %lev_bottom = lev.ylow;
        fcl = meshColor(l, :);
        ecl = edgeColor(l, :);
        fcl(4) = opacity(l); % opacity
        %rectangle('Position', [thisBox.ylow thisBox.xlow  thisBox.yhi thisBox.xhi ],...
        rectangle('Position', [thisBox.xlow thisBox.ylow  thisBox.xhi thisBox.yhi ],...
            'FaceColor', fcl, 'EdgeColor', ecl, 'LineWidth', 2.0);
    
    end
    
end
hold off






axPsi.Visible = 'off';


cPsi.Label.String = '\psi';
maxPsi = max(max(psi));
%minPsi = min(min(psi));
%deltaPsi = maxPsi-minPsi;
%caxis([minPsi+0.01*deltaPsi maxPsi-0.01*deltaPsi]);
currentTicks = cPsi.Ticks;
cPsi.Ticks = [currentTicks(1) currentTicks(end)];
format = '%1.1f';
cPsi.TickLabels = {sprintf(format,currentTicks(1)), sprintf(format,currentTicks(end))};

cPsi.Position = [oldPos(1) oldPos(2)-0.05 oldPos(3)-0.035 oldPos(4)];

oldPos = cPsi.Label.Position;
cPsi.Label.Position = [oldPos(1) oldPos(2) + 1.2];

axPsi.XTick = [0 1];
axPsi.YTick = [0 1];


axTemperature.XTick = [0.02 0.98];
axTemperature.YTick = [0.02 0.98];
axTemperature.XTickLabel = {'0', '1'};
axTemperature.YTickLabel = {'0', '1'};


text(0.04, 0.92, '(b)', 'FontSize', 16, 'Color', [1 1 1]);
%text(-0.2, 1.2, '(b)', 'FontSize', 16);

subplot(m, n, 3);
hold on
plot(log10(dz), log10(errUniform), '-x');
plot(log10(dz), log10(errAMR1lev), '-x');
%plot(log10(dz), log10(errUniform), '-x');
plot([-2 -1], [-3.5 -1.5], '--');
hold off

%xlim([-6 0]);
%ylim([-6 2]);

legend({'Uniform', '$n_{ref}=2$'}, 'Location', 'northwest');

xlabel('$log_{10}(\Delta z)$');
ylabel('Error');

box on; 

axErr = gca;
axErr.Position =  [0.76 axBottom 0.2 axHeight];


text(axErr.XLim(2)-0.5, axErr.YLim(1)+0.3, '(c)', 'FontSize', 16);


set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
set(h, 'InvertHardCopy', 'off'); % keep white text white

if saveFigure
    fprintf('Saved to %s \n', figureName);
    print(h,[figureName, '.pdf'],'-dpdf','-r50')
    print(h,[figureName, '.png'],'-dpng','-r1000')
    
end

end


function err = getErr(dataFolder, prefixFunc, gridRes, highResML)
err = NaN*ones(length(gridRes), 1);

exactNu = 1.41008;

for f =1:length(gridRes)
    thisRes = gridRes(f);
    thisFolder = [dataFolder, prefixFunc(thisRes)];
    
    thisML = getFinalPlotFile([dataFolder, prefixFunc(thisRes)]);
    
    if length(thisML.levelArray) > 0
    MLdiff =  highResML.diff(thisML, [thisML.components.Temperature], ...
                            [thisML.components.Temperature]);
             
                    
    [L1, L2, Max, Sum] = AMRSum(MLdiff, thisML.components.Temperature);
        
    err(f) = L1;
    end
    
    
    %diags = getDiagnostics([thisFolder, '/diagnostics.out']);
    %fprintf('Loading diags %s \n', thisFolder);
    %if isfield(diags,'Nusselt')
    %Nuerr = abs(diags.Nusselt(end) - exactNu)/exactNu;   
    %err(f) = Nuerr;
    %end
    
    
   %errFile = [thisFolder, '/err.mat'];
   % if exist(errFile, 'file') == 2
   %     load(errFile)
   % else
        %output = getFinalPlotFile(thisFolder);
       
        %if length(output.levelArray) > 0
        %    Terr = output.dataForComp(output.components.Terr);
        %    Te = squeeze(Terr(2,:));
        %    e.meanTerr =  mean(abs(Te));
        %    save(errFile, 'e');
        %end
   % end
    
   % err(f) = e.meanTerr;
end

end

function f = plotPrefixAMR(N)
f =  ['AMR-Subcycle-Reflux-Freestream0.45-MaxLevel1-convectionDarcyBrinkman-',num2str(N),'-ref2--0'];
end

function f = plotPrefixUniform(N)
f =  ['Uniform-convectionDarcyBrinkman-',num2str(N),'--0'];
end


