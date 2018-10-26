%clear all;
function viewWideDomainTimeSeriesEGU

close all;

set(0, 'defaultlinelinewidth',2);
set(0, 'defaultaxeslinewidth',2);
set(0, 'defaultpatchlinewidth',2);
set(0, 'defaultAxesFontSize',18);
tic;

%dataFolder = getDataDir('middleton/'); %'/media/parkinsonjl/FREECOM HDD/';
%dataFolder = '/home/parkinsonjl/mnt/raymaster/convection-in-sea-ice/MushyLayer/matlab/EGUPoster/';
dataFolder = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/matlab/EGUPoster/';

plotHeight = 400;
plotWidth = 1000;

h = figure();
set(h, 'Position', [100 100 plotWidth plotHeight]);

Slmax = 0.085; %55g/kg

axPos = [0.5 0.3 0.35 0.6]
bcPos = [0.05 0.21 0.35 0.7];


% BC diagram
axBCs = axes;

yTop = 2; yBottom = 0;
xBottom = 0; xTop = 2.5;
xlim([xBottom xTop]);
ylim([yBottom yTop]);
daspect([1 1 1]);
text(0.1,yTop+0.13, '$\mathbf{U}\cdot \mathbf{n} = 0, \, T = T_a, \, S=35$g/kg', 'FontSize', 16);

%text(0.1,yBottom+0.15, '$\mathbf{U}\cdot \mathbf{n} = 0, \, T = -2^\circ$ C$, \, \partial S / \partial z = 0$', 'FontSize', 16);

text(0.1,yBottom+0.45, {'Inflow:', '$T = -2^\circ$C', '$S=35$g/kg', '$\partial \mathbf{U}/\partial z = 0$' }, 'FontSize', 16);
%text(xTop-0.8,yBottom+0.45, {'Outflow:', '$\partial  T / \partial  z = 0$', '$\partial  S / \partial  z = 0$', '$\partial \mathbf{U}/\partial z = 0$' }, 'FontSize', 16);
text(xTop-1.2,yBottom+0.45, {'Outflow:', '$T:$ extrapolation', '$S:$ extrapolation', '$\partial \mathbf{U}/\partial z = 0$' }, 'FontSize', 16);

arrowY = [0.13 0.23];
annotation('arrow', [0.12 0.12], arrowY);
annotation('arrow', [0.3 0.3], fliplr(arrowY));


text(0.15,2*(yBottom+yTop)/3, 'Periodic', 'Rotation', 90, 'FontSize', 16);
text(xTop-0.15,2*(yBottom+yTop)/3, 'Periodic', 'Rotation', 90, 'FontSize', 16);

text((xTop+xBottom)/5, (yBottom+yTop)/1.5, '$\textit{Axes not to scale}$', 'FontSize', 16);


axBCs.XTick = [xBottom xTop];
axBCs.YTick = [yBottom yTop];
axBCs.YTickLabel = {'-8', '0'};
axBCs.XTickLabel = {'0', '2'};
xlab = xlabel('$x$ (cm)');
ylab = ylabel('$z$ (cm)');

xlabPos = xlab.Position; ylabPos = ylab.Position;
xlab.Position = [xlabPos(1) xlabPos(2) +0.1];
ylab.Position = [ylabPos(1) + 0.08 ylabPos(2)];

box on;



axBCs.Position = bcPos;



% Diffusive solution diagram

[axPorosity, axSl] = doPlot([dataFolder, 'T-15/toPlot1200.mat'], 1.5, true, axPos, Slmax);

axPorosity.XTick = [0.05 1.95];
axPorosity.XTickLabels = {'0', '2'};
axPorosity.XLabel.String = '$x$ (cm)';

cbSl = colorbar(axSl, 'Location', 'eastoutside');
cbChi = colorbar(axPorosity, 'Location', 'eastoutside');

cbSlPos = cbSl.Position; cbChiPos = cbChi.Position;
cbSl.Position  = [cbSlPos(1)+0.08   cbSlPos(2)-0.18 cbSlPos(3)  0.35];
cbChi.Position = [cbChiPos(1)+0.08 cbChiPos(2)+0.28 cbChiPos(3) 0.35];

cbSl.Label.String = 'S_l (g/kg)';
cbChi.Label.String = '\chi';


cbSl.Ticks = [0 Slmax];
cbSl.TickLabels = {'35', num2str(round(35+Slmax*233))};

cbChi.Ticks = [0.0 0.98];
cbChi.TickLabels = {'0', '1'};

%title('t = 2 mins');



 h =gcf;
 h.InvertHardcopy = 'off';
 h.Color = 'white';
  set(h,'Units','Inches');
  pos = get(h,'Position');
  set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
  filename = ['/home/parkinsonjl/convection-in-sea-ice/figures', '/EGUBCs'];
  print(h,[filename, '.pdf'],'-dpdf','-r0')
  print(h,[filename, '.png'],'-dpng','-r400')

end


function [axPorosity, axSl] = doPlot(dataFile, depth, doLabel, axPos, SlMax)

blues = makeColorMap([8 48 107]/255, [104, 172, 212]/255, [1 1 1]); blues = flipud(blues);

%slMap = makeColorMap([  0.9763    0.9831    0.0538] , [ 0.0488    0.5772    0.8228], [8 48 107]/255); 
slMap = makeColorMap([  0.9763    0.9831    0.0538] , [ 0.6488    0.7772    0.2228], [8 48 107]/255); 

axPorosity = axes;

load(dataFile);

zmin = 900;
zmax=1024;

zmin = round(1024*(1-depth/8));


Xplot = X(zmin:zmax, :);
Zplot = Z(zmin:zmax, :);
chiPlot = porosity(:, zmin:zmax).';
SlPlot = Sl(:, zmin:zmax).';

pcolor(Xplot,Zplot,chiPlot);

colormap(axPorosity, (blues));

axPorosity.XTick = [];
if doLabel
    axPorosity.YTick = [min(min(Z(zmin:zmax, :))) max(max(Z)) ];
    bottomLabel = sprintf('%1.1f', -depth);
    axPorosity.YTickLabels = { bottomLabel, '0.0'};
    ylab = ylabel('$z$ (cm)');
    oldPos = ylab.Position;
    %ylab.Position = [oldPos(1)+0.08 oldPos(2)];
else
    axPorosity.YTick = [];
end

daspect([1 1 1]);

axPorosity.Position = axPos;

caxis([0 1]);

%Sl axes
axSl = axes;

smoothTransition = false;

%SlMax = 0.16;
SlDiff  =  min(SlPlot+(1-SlMax), 0);

if smoothTransition
SlDiff(chiPlot<0.8) = NaN;
else
    SlDiff(chiPlot<1.0) = NaN;
end

SlDiff = -min(min(SlDiff))+SlDiff;



colormap(axSl, flipud(slMap));
%minSl = -1;
avPorosity = mean(chiPlot.^6, 1);

chiFilter = repmat(avPorosity,256,1); %porosity.^2;

if smoothTransition
SlDiff =SlDiff.*chiFilter;
end

%min(min(SlDiff(:, zmin:zmax)))

pcolor(Xplot, Zplot, SlDiff);
%caxis([min(min(SlDiff)) max(max(SlDiff))]);
caxis([0 SlMax]);
%colorbar('Location', 'eastoutside');

daspect([1 1 1]);

axis(axSl, 'off');
axSl.Visible = 'off';

axSl.Position = axPos;  

linkaxes([axPorosity axSl])

timeStr = sprintf('$t=%1.1f$ mins', t*800/70);
text(0.05, 8-depth+0.05,timeStr, 'Color', [1 1 1], 'FontSize', 24);


end
