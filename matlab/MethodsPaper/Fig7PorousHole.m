% Plot  typical AMR solutions for darcy brinkman problems 
function Fig7PorousHole(dataFolder, figureDirectory)

if nargin < 1
   % dataFolder = getDataDir('AMRConvergenceTest/MushyConvectionLiquid2-t0.5/');
   base_dir = '/home/parkinsonjl/mnt/sharedStorage/TestFinal/';
   dataFolder = fullfile(base_dir, '/PorousMushyHole-t0.00015-hole0.04/');
end

if nargin < 2
    figureDirectory = '/home/parkinsonjl/mnt/sharedStorage/TestDiffusiveTimescale/';
end

close all;

figureName =  fullfile(figureDirectory, 'Fig7PorousHoleFigure'); % thisFilename; %[dataFolderNu, 'benchmark2Convergence'];
saveFigure = true;

%highRes = getFinalPlotFile([dataFolderNu, plotPrefixUniform(gridRes(end)*2)]);

h = figure();
%h.Position = [200 200 1200 500];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(h,'Units','Inches');
h.Position = [2.0 2.0 6.5 2.7];
% Caption font size is 8pt
textFontSize = 8;
legendFontSize = 8;
domainFontSize = 8;

set(0, 'defaultlinelinewidth',1);
set(0, 'defaultaxeslinewidth',1);
set(0, 'defaultpatchlinewidth',1);
set(0, 'defaultAxesFontSize', textFontSize);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axHeight = 0.55;
axBottom = 0.25;

m = 1; n = 3;
subplot(m, n, 1)
axBCs = gca;
axBCs.Position = [0.05 axBottom 0.25 axHeight];


xlim([0 1]);
ylim([0 1]);
box on;

axBCs.YTick = [0 1];
axBCs.XTick = [0 1];
daspect([1 1 1]);


text(0.04, 0.5, {'Eqn''s 3-6.', 'All boundaries are periodic.'}, 'FontSize', domainFontSize);

periodicStr = {'Periodic'};

xl = xlabel('$x$');
yl = ylabel('$z$');

text(-0.2, 1.2, '(a)', 'FontSize', textFontSize);


subplot(m, n, 2)

%AMRsol = getFinalPlotFile([dataFolder, plotPrefixUniform(32)]);
%AMRsol = getFinalPlotFile([dataFolderNu, plotPrefixAMR(32)]);
UniformFineSol = getInitialPlotFile(fullfile(dataFolder, plotPrefixUniform(64)));
%makeSubplot(AMRsol, UniformFineSol, [0.35 axBottom 0.3 axHeight])
makeSubplot(UniformFineSol, UniformFineSol, [0.35 axBottom 0.3 axHeight])


text(0.16, 0.85, '(b)', 'FontSize', textFontSize);
%text(-0.2, 1.2, '(b)', 'FontSize', 16);

subplot(m, n, 3);

AMRsol = getFinalPlotFile(fullfile(dataFolder, plotPrefixAMR(64, dataFolder)));
UniformFineSol = getFinalPlotFile(fullfile(dataFolder,  plotPrefixUniform(64)));

makeSubplot(AMRsol, UniformFineSol, [0.7 axBottom 0.3 axHeight])

axErr = gca;

 text(0.16, 0.85, '(c)', 'FontSize', textFontSize, 'Color', [0 0 0]);

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%set(h, 'InvertHardCopy', 'off'); % keep white text white

if saveFigure
    fprintf('Saved to %s \n', figureName);
   % set(gcf,'renderer','Painters');
   % This exports better as a pdf
    %print(h,[figureName, '.eps'],'-depsc','-r50')
    print(h,[figureName, '.pdf'],'-dpdf','-r50')
    %  print(h,[figureName, '.png'],'-dpng','-r1000')
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

function f = plotPrefixAMR(N, dataFolder)
% Use dataFolder to determine the freestream value

folders = dir(fullfile(dataFolder, 'AMR-Subcycle-Reflux-Freestream*'));

for i=1:length(folders)
    fname = folders(i).name;
    if contains(fname, ['-',num2str(N),'--0'])
       f = fname;
       return;
    end
end
%f =  ['AMR-Subcycle-Reflux-Freestream0.95-MaxLevel1-ref2-PorousMushyHole-',num2str(N),'--0'];
end

function f = plotPrefixUniform(N)
f =  ['Uniform-PorousMushyHole-',num2str(N),'--0'];
end



function makeSubplot(AMRsol,UniformFineSol, axisExtent)

onLaptop = false;
[~, name] = system('hostname'); 
name = strtrim(name);
if strcmp(name, 'atmlxlap005')
onLaptop = true;
end
% Get matlab version
[v, d] = version;
year = str2num(d(end-4:end));
if year > 2016
    recentMatlab = true;
else
    recentMatlab = false;
end

%T = AMRsol.dataForComp(AMRsol.components.Temperature).';

%psi = AMRsol.dataForComp(AMRsol.components.streamfunction).';
%psiUniform = UniformFineSol.dataForComp(UniformFineSol.components.streamfunction).';

porosity = AMRsol.dataForComp(AMRsol.components.Porosity).';
porosityUniform = UniformFineSol.dataForComp(UniformFineSol.components.Porosity).';

%Tuniform = UniformFineSol.dataForComp(UniformFineSol.components.Temperature).';
Sl = AMRsol.dataForComp(AMRsol.components.Liquidconcentration).';
SlUniform = UniformFineSol.dataForComp(AMRsol.components.Liquidconcentration).';


Uuni = UniformFineSol.dataForComp(AMRsol.components.xDarcyvelocity).';
Vuni = UniformFineSol.dataForComp(AMRsol.components.yDarcyvelocity).';

% if max(max(T)) > 1.5
%    T = T/2; 
%    Tuniform =Tuniform/2;
% end

[X,Y] = AMRsol.grid();
[XUniform,YUniform] = UniformFineSol.grid();
%Create smooth psi

Xbig = XUniform; Ybig = YUniform; Vbig = Vuni; Ubig =Uuni;

% Zoom region:
refinement = 2;
X = zoomInto(X, refinement);
Y = zoomInto(Y, refinement);
XUniform = zoomInto(XUniform, refinement);
YUniform = zoomInto(YUniform, refinement);
porosity = zoomInto(porosity, refinement);
SlUniform = zoomInto(SlUniform, refinement);
Sl = zoomInto(Sl, refinement);

Uuni = zoomInto(Uuni, refinement);
Vuni = zoomInto(Vuni, refinement);

porosityUniform = zoomInto(porosityUniform, refinement);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Liquid Salinity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axSl = gca;
plot_Sl = SlUniform;

minSl = min(min(plot_Sl))*0.98;
maxSl = max(max(plot_Sl))*1.01;

%plot_TUniform = Tuniform*2 - 1;
v = linspace(minSl, maxSl,5);

%[CTemp, hTemp] = contour(XUniform, YUniform, plot_Sl, v);
%hTemp.LineWidth = 2.0;
pcolor(X,Y,Sl); 
daspect([1 1 1]);

%pcolor(X,Y,T*2-1); 
%pcolor(X,Y,psiUniform); 

%daspect([1 1 1]);

if recentMatlab
    SlLabOffset = -1.2;
else
    SlLabOffset = +1.0; 
end

colormap(axSl, flipud(bluewhitered(257)));

cSl = colorbar(axSl, 'Location', 'northoutside');
cSl.Ticks = [minSl*0.98 maxSl*1.01];
cSl.TickLabels = {sprintf('%1.1f', minSl), sprintf('%1.1f', maxSl)};
cSl.Label.String = '\Theta_l';

%cSl.Position(4) = cSl.Position(4)*SlHeightScale;
 
oldPos = cSl.Label.Position;
cSl.Label.Position = [(minSl+maxSl)/2 oldPos(2)+SlLabOffset];

box on;
%set(axTemperature, 'Layer', 'top');

xl = xlabel('$x$');
yl = ylabel('$z$');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Porosity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axPorosity = axes;
 
daspect([1 1 1]);

v = linspace(0.7, 1.0, 6);
[Cchi, hChi] = contour(XUniform, YUniform, porosityUniform, v);
hChi.LineWidth = 1.0;

daspect([1 1 1]);

caxis([0.7 1.0]);

box on;

cPorosity = colorbar(axPorosity,'Location', 'southoutside');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axVel = axes;

xlimits = [0.25 0.75];
ylimits = [0.25,0.75];

xlim(xlimits);
ylim(ylimits);

daspect([1 1 1]);
axVel.Visible = 'off';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combine axes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
linkaxes([axSl,axPorosity]);
linkaxes([axSl, axVel]);

axPorosity.Position = axisExtent;
axSl.Position = axisExtent;
axVel.Position = axisExtent;

% Set pcolor on top
set(axSl, 'Layer', 'top');


cSl.Position(4) = 0.05;
cPorosity.Position(4) = 0.05;


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
% for l = 2:length(AMRsol.levelArray)
%     lev = AMRsol.levelArray(l);
%     levDx = 0.5*lev.dx;
%     
%     for b = 1:length(lev.boxArray)
%         thisBox = lev.boxArray(b);
%     
%         %lev_width = lev.xhi-lev.xlow+2*levDx;
%         %lev_height = lev.yhi-lev.ylow+2*levDx;
%         %lev_bottom = lev.ylow;
%         fcl = meshColor(l, :);
%         ecl = edgeColor(l, :);
%         fcl(4) = opacity(l); % opacity
%         %rectangle('Position', [thisBox.ylow thisBox.xlow  thisBox.yhi thisBox.xhi ],...
%         rectangle('Position', [thisBox.xlow+levDx    thisBox.ylow  ...
%             thisBox.xhi-thisBox.xlow-levDx thisBox.yhi-thisBox.ylow ],...
%             'FaceColor', fcl, 'EdgeColor', ecl, 'LineWidth', 2.0);
%     
%     end
%     
% end

% New method
if exist('polyshape')
   lev1 = AMRsol.levelArray(1);
   lev1Dx = 0.5*lev1.dx;
    
    domainBox = polyshape([xlimits(1)+lev1Dx xlimits(1)+lev1Dx xlimits(end)-lev1Dx xlimits(end)-lev1Dx], ...
        [ylimits(1)+lev1Dx ylimits(end)-lev1Dx ylimits(end)-lev1Dx ylimits(1)+lev1Dx]);
    
    meshes = AMRsol.getMeshes(domainBox);
    
    % Plot meshes
    if length(meshes) > 1
        pshape = meshes(2);
        mesh2 = plot(pshape); % meshes on chombo level 1

        %plot(domainBox);
        mesh2.FaceAlpha = 0;
        mesh2.EdgeColor = [1 0 1]; % magenta
        mesh2.LineWidth = 2;
    end
    
    if length(meshes) > 2
        pshape = meshes(3);
                
        mesh3 = plot(pshape); % meshes on chombo level 2
        mesh3.FaceAlpha = 0;
        mesh3.EdgeColor = [0 1 1]; %cyan
        mesh3.LineWidth = 2;
    end
    
else
    fprintf('Not plotting level outlines as polyshape cannot be found \n');
end


hold off

pcAxis = axSl;
contAxis = axPorosity;

%axPsi.Visible = 'off';
%axSl.Visible = 'off';
contAxis.Visible = 'off';

cPorosity.Label.String = '\chi';
%maxPsi = max(max(psi));
%minPsi = min(min(psi));
%deltaPsi = maxPsi-minPsi;
%caxis([minPsi+0.01*deltaPsi maxPsi-0.01*deltaPsi]);
currentTicks = cPorosity.Ticks;
newTicks = currentTicks;

newTicks = cPorosity.Limits;
cPorosity.Ticks = [newTicks(1) newTicks(end)];


box on;


% if maxPsi < 1e-3
%     formatMax = '%1.1e';
% else
%     formatMax = '%1.1f';
% end
% 
% if newTicks(1) == 0
%     formatMin = '%1.1f';
% else
%     formatMin = formatMax;
% end

if onLaptop
    PorosityLabOffset = 1.2;
    PorosityColorbarShift = -0.15;
else
    PorosityLabOffset = 1.2;
    PorosityColorbarShift = -0.15;
end


formatMin = '%1.1f';
formatMax = '%1.1f';

cPorosity.TickLabels = {sprintf(formatMin,newTicks(1)), sprintf(formatMax,newTicks(end))};

oldPos = cPorosity.Position;
thetaColorbarPos = cSl.Position;
%cPsi.Position = [oldPos(1) oldPos(2)-0.05 oldPos(3)-0.035 oldPos(4)];
cPorosity.Position = [thetaColorbarPos(1) axisExtent(2)+PorosityColorbarShift thetaColorbarPos(3) oldPos(4)];

oldPos = cPorosity.Label.Position;


cPorosity.Label.Position = [oldPos(1) oldPos(2) + PorosityLabOffset];

minX = min(min(X));
maxX = max(max(X));

minY = min(min(Y));
maxY = max(max(Y));

pcAxis.XTick = [minX maxX];
pcAxis.YTick = [minY maxY];


%axPorosity.XTickLabel = {sprintf('%1.2f', minX), sprintf('%1.2f', maxX)};
%axPorosity.YTickLabel = {sprintf('%1.2f', minY), sprintf('%1.2f', maxY)};

pcAxis.XTickLabel = {'0.25', '0.75'};
pcAxis.YTickLabel = {'0.25', '0.75'};



xl.Position = [0.5 0.24];
yl.Position = [0.24 0.5];

%%axSl.XTick = [0.02 0.98];
%axSl.YTick = [0.02 0.98];
%axSl.XTickLabel = {'0', '1'};
%axSl.YTickLabel = {'0', '1'};

end


function X = zoomInto(X, refinement)

s = size(X);
newS = round(s/refinement);

ax = 0.5*s(1)*(1-1/refinement);
ay = 0.5*s(2)*(1-1/refinement);

xl = 1 + ax:1:(ax + newS(1));
yl = 1 + ay:1:(ay + newS(2));

X = X(xl, yl);

temp = 0;

end