% Plot  typical AMR solutions for a mushy layer test problem
function MushySimplePlot

clear all; close all;

%gridRes = [8 16 32 64];
coarseResolution = 16;


dataFolderAMR = getDataDir(sprintf('AMRConvergenceTest/MushyDB/AMR-Subcycle-Reflux-Freestream0.99-MaxLevel1-ref2-MushyDB-%d--0/', coarseResolution));
dataFolderUniform = getDataDir(sprintf('AMRConvergenceTest/MushyDB/Uniform-MushyDB-%d--0/', coarseResolution*2));

thisFilename = mfilename('fullpath');
thisFolder = strrep(thisFilename, 'MushySimplePlot', '');
figureName =  [thisFolder, 'simpleMushFigure']; %[dataFolderNu, 'benchmark2Convergence'];
saveFigure = false;


uniformPrefix = sprintf('Uniform-MushyDB-%d-',coarseResolution*2);
amrPrefix = sprintf('AMR-Subcycle-Reflux-Freestream0.99-MaxLevel1-ref2-MushyDB-%d-',coarseResolution);


earlyFrameAMR = 200;
%highRes = getFinalPlotFile([dataFolderNu, plotPrefixUniform(gridRes(end)*2)]);


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


text(0.1, 0.55, {'Eq''s 3-6 and 8'}, 'FontSize', 16);

insulatingStr = '$\partial H/\partial x = \partial \Theta / \partial x = 0, \mathbf{U} = 0.$';
text(-0.1, 0.0, insulatingStr, 'FontSize', 16, 'Rotation', 90);
text(1.12 , 0.0, insulatingStr, 'FontSize', 16, 'Rotation', 90);


%text(-0.1, 0.0, '$H=H(\theta=1, \chi=1), \Theta=0, \mathbf{U} = 0.$', 'FontSize', 16, 'Rotation', 90);
%text(1.12 , 0.0, '$H=H(\theta=0, \chi=1), \Theta=0, \mathbf{U} = 0.$', 'FontSize', 16, 'Rotation', 90);

%text(0.0, 1.13 ,insulatingStr, 'FontSize', 16);
%text(0.0, -0.15, insulatingStr, 'FontSize', 16);

text(0.0, 1.06, '$H=0, \Theta=-1, \mathbf{U} = 0.$', 'FontSize', 16);
text(0.0, -0.15, {'$H=H(\theta=1.3, \chi=1), \Theta=-1,$','$\mathbf{U}$: Inflow/outflow'}, 'FontSize', 16);

%axis arrows
annotation('arrow',[0.03 0.07], [0.06 0.06]);
annotation('arrow',[0.03 0.03], [0.06 0.16]);
text(-0.05, -0.4, '$x$', 'FontSize', 16);
text(-0.23 , -0.2   , '$z$', 'FontSize', 16);


text(0.04, 0.93, '(a)', 'FontSize', 16);


subplot(m, n, 2)

%AMRsol = getFinalPlotFile([dataFolder, plotPrefixUniform(32)]);
AMRsol = MushyLayerOutput(2, earlyFrameAMR, dataFolderAMR, amrPrefix, true); %getFinalPlotFile(dataFolderAMR);
UniformFineSol = MushyLayerOutput(2, earlyFrameAMR*2, dataFolderUniform, uniformPrefix, true); %getFinalPlotFile(dataFolderUniform);
makeSubplot(AMRsol, UniformFineSol, [0.35 axBottom 0.3 axHeight])

midTime = sprintf('t=%1.2f', UniformFineSol.t);

%text(0.04, 0.92, '(b)', 'FontSize', 16, 'Color', [0 0 0]);
text(0.00, 0.54, '(b)', 'FontSize', 16, 'Color', [0 0 0]);
text(0.20, 0.54, midTime, 'FontSize', 20, 'Color', [0 0 0]);

%text(-0.2, 1.2, '(b)', 'FontSize', 16);

subplot(m, n, 3);

AMRsol = getFinalPlotFile(dataFolderAMR);
UniformFineSol = getFinalPlotFile(dataFolderUniform);
makeSubplot(AMRsol, UniformFineSol, [0.35 axBottom 0.3 axHeight])


% text(axErr.XLim(2)-0.5, axErr.YLim(1)+0.3, '(c)', 'FontSize', 16);
text(0.04, 0.92, '(c)', 'FontSize', 16, 'Color', [0 0 0]);

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%set(h, 'InvertHardCopy', 'off'); % keep white text white

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
f =  ['VariableMesh2SubcycleRefluxFreestream0.45-convectionDB-',num2str(N),'-ref2--0'];
end

function f = plotPrefixUniform(N)
f =  ['Uniform-convectionDB-',num2str(N),'--0'];
end



function makeSubplot(AMRsol,UniformFineSol, axisExtent)

T = AMRsol.dataForComp(AMRsol.components.Temperature).';
psi = AMRsol.dataForComp(AMRsol.components.streamfunction).';
porosity = AMRsol.dataForComp(AMRsol.components.Porosity).';


psiUniform = UniformFineSol.dataForComp(UniformFineSol.components.streamfunction).';
Tuniform = UniformFineSol.dataForComp(UniformFineSol.components.Temperature).';
Porosityuniform = UniformFineSol.dataForComp(UniformFineSol.components.Porosity).';

if max(max(T)) > 1.5
   T = T/2; 
   Tuniform =Tuniform/2;
end

[X,Y] = AMRsol.grid();
[XUniform,YUniform] = UniformFineSol.grid();
%Create smooth psi


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Porosity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot psi uniform because it's smoother
axPorosity = gca;
%Ncontours = min(10, round(max(max(psiUniform))/1.5e-5));
%v =  linspace(4e-6,max(max(psiUniform)), Ncontours) ;

%[Cpsi, hpsi] = contour(XUniform,YUniform,psiUniform, v); 
%[Cpsi, hpsi] = contour(XUniform,YUniform,plot_T); 
pcolor(X,Y,porosity); 

daspect([1 1 1]);

cmap = flipud(plasma(100));
cmap = flipud(viridis);
cmap = blues();
colormap(axPorosity,cmap);

box on;

cPorosity = colorbar(axPorosity,'Location', 'southoutside');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Psi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%psiUniform = max(psiUniform, 0);
%psi = max(psi, 0);
% Plot psi uniform because it's smoother
axPsi = axes;
%Ncontours = min(10, round(max(max(psiUniform))/1.5e-5));
%v =  linspace(4e-6,max(max(psiUniform)), Ncontours) ;
maxPsi = max(max(psi));
Ncontours = 11; %9 or 17;
minContour = maxPsi/(Ncontours/2);
minContour = round(minContour, 1, 'significant');

%maxPsiRounded = round(maxPsi, 1, 'significant');
%v = linspace(-maxPsiRounded, maxPsiRounded, Ncontours);

newMaxPsi = minContour*((Ncontours-1)/2);
v = linspace(-newMaxPsi, newMaxPsi, Ncontours);

% Remove 0 contour
vnew = v;
vnew(round(Ncontours/2):end-1) = v(round(Ncontours/2)+1:end);
vnew = vnew(1:end-1);
v = vnew;

[Cpsi, hpsi] = contour(XUniform,YUniform,psiUniform, v); 
%[Cpsi, hpsi] = contour(XUniform,YUniform,plot_T); 
%pcolor(X,Y,psi); 
hpsi.LineColor = [0 0 0];
hpsi.LineColor = [1 1 1];
hpsi.LineWidth = 2.0;
%hpsi.LineStyle = '--';

starti = floor(Ncontours/2-1);
endi = ceil(Ncontours/2+1);

%clabel(Cpsi, hpsi, v(starti:endi), 'Color', 'white', 'FontWeight', 'bold');
clabel(Cpsi, hpsi, v(starti:endi), 'LabelSpacing', 2000,  'Color', 'white', 'FontWeight', 'bold');

daspect([1 1 1]);

%cmap = flipud(plasma(100));
%cmap = flipud(viridis);
%colormap(axPsi,cmap);

box on;

%cPsi = colorbar(axPsi,'Location', 'southoutside');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Temperature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%axTemperature = axes;
% 
% plot_T = T*2-1;
% plot_TUniform = Tuniform*2 - 1.3;
% [CTemp, hTemp] = contour(XUniform, YUniform, plot_TUniform);
% %hTemp.LineWidth = 2.0;
% 
% daspect([1 1 1]);
% colormap(axTemperature, bluewhitered);
% cTemp = colorbar(axTemperature, 'Location', 'northoutside');
% cTemp.Ticks = cTemp.Limits;
% cTemp.TickLabels = {'0', '1.3'};
% cTemp.Label.String = '\theta';
% oldPos = cTemp.Label.Position;
% cTemp.Label.Position = [oldPos(1) oldPos(2)-1.2];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combine axes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
linkaxes([axPsi, axPorosity]); %axTemperature
axPsi.Position = axisExtent;
%axTemperature.Position = axisExtent;
axPorosity.Position = axisExtent;


set(axPorosity, 'Layer', 'top');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot meshes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        
       % rectangle('Position', [thisBox.xlow+levDx    thisBox.ylow  ...
       %     thisBox.xhi-thisBox.xlow-levDx thisBox.yhi-thisBox.ylow ],...
       %     'FaceColor', fcl, 'EdgeColor', ecl, 'LineWidth', 2.0);
       
       xPoints = [thisBox.xlow-levDx thisBox.xlow-levDx thisBox.xhi+levDx thisBox.xhi+levDx];
       yPoints = [thisBox.ylow-levDx thisBox.yhi+levDx thisBox.yhi+levDx thisBox.ylow-levDx];
       
       thisPoly = polyshape(xPoints, yPoints);
       %polygons(end+1) = thisPoly;
       
       if b==1
           mergedPolyshape = thisPoly;
       else
           mergedPolyshape = union(mergedPolyshape, thisPoly);
       end
       
       %plot(thisPoly);
    %
    end
    
end

%mergedPolyshape = polygons(1);%

%for i=2:length(polygons)
%    mergedPolyshape = union(mergedPolyshape, polygons(i));
%end

fcl = meshColor(l, :);
        ecl = edgeColor(l, :);
        fcl(4) = opacity(l); % opacity
        
       % rectangle('Position', [thisBox.xlow+levDx    thisBox.ylow  ...
       %     thisBox.xhi-thisBox.xlow-levDx thisBox.yhi-thisBox.ylow ],...
       %     'FaceColor', fcl, 'EdgeColor', ecl, 'LineWidth', 2.0);

plot(mergedPolyshape, 'EdgeColor', ecl, 'FaceColor', fcl, 'LineWidth', 2.0, 'FaceAlpha',0.0);

hold off

axPsi.Visible = 'off';

%axTemperature.Visible = 'off';

% Setup Psi colorbar
% cPsi.Label.String = '\psi';
% maxPsi = max(max(psi));
% currentTicks = cPsi.Ticks;
% newTicks = currentTicks;
% 
% newTicks = cPsi.Limits;
% cPsi.Ticks = [newTicks(1) newTicks(end)];
% 
% if maxPsi < 1e-3
%     format = '%1.1e';
% else
%     format = '%1.1f';
% end
% cPsi.TickLabels = {sprintf(format,newTicks(1)), sprintf(format,newTicks(end))};
% 
% oldPos = cPsi.Position;
% thetaColorbarPos = cTemp.Position;
% %cPsi.Position = [oldPos(1) oldPos(2)-0.05 oldPos(3)-0.035 oldPos(4)];
% cPsi.Position = [thetaColorbarPos(1) axisExtent(2)-0.15 thetaColorbarPos(3) oldPos(4)];
% 
% oldPos = cPsi.Label.Position;
% cPsi.Label.Position = [oldPos(1) oldPos(2) + 1.2];

% Setup porosity colobar
cPorosity.Label.String = '\chi';
currentTicks = cPorosity.Ticks;
newTicks = currentTicks;

newTicks = cPorosity.Limits;
cPorosity.Ticks = [newTicks(1) newTicks(end)];

format = '%1.2f';

cPorosity.TickLabels = {sprintf(format,newTicks(1)), sprintf(format,newTicks(end))};

oldPos = cPorosity.Position;
%thetaColorbarPos = cTemp.Position;
%cPsi.Position = [oldPos(1) oldPos(2)-0.05 oldPos(3)-0.035 oldPos(4)];
%cPorosity.Position = [thetaColorbarPos(1) axisExtent(2)-0.15 thetaColorbarPos(3) oldPos(4)];

oldPos = cPorosity.Label.Position;
cPorosity.Label.Position = [oldPos(1) oldPos(2) + 1.2];

%axPsi.XTick = [0 1];
%axPsi.YTick = [0 1];

axPorosity.XTick = [0.0 0.5];
axPorosity.YTick = [0.0 0.5];
axPorosity.XTickLabel = {'0', '0.5'};
axPorosity.YTickLabel = {'0', '0.5'};

% Setup temperature axes
%axTemperature.XTick = [0.02 0.98];
%axTemperature.YTick = [0.02 0.98];
%axTemperature.XTickLabel = {'0', '1'};
%axTemperature.YTickLabel = {'0', '1'};

end