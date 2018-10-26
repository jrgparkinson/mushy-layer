%clear all;
close all;

set(0, 'defaultlinelinewidth',2);
set(0, 'defaultaxeslinewidth',2);
set(0, 'defaultpatchlinewidth',2);
set(0, 'defaultAxesFontSize',22);

dataFolder = '/home/parkinsonjl/mnt/sharedStorage/optimalStates-highRes-new/';
doColorbar = false;


plot_prefix = 'CR1.1RaC200Le200ChiCubedPermeabilitypts96-0';
plot_prefix = 'CR1.1RaC150Le200ChiCubedPermeabilitypts64-0';
plot_prefix = 'CR1.15RaC150Le200ChiCubedPermeabilitypts64-0';
plot_prefix = 'CR1.15RaC200Le200ChiCubedPermeabilitypts60-0';
plot_prefix = 'CR1.1RaC100Le200ChiCubedPermeabilitypts64-0';

output = getFinalPlotFile([dataFolder, plot_prefix]);

%output = MushyLayerOutput(dim, frame, output_dir, plot_prefix, subcycled);

perm = output.dataForComp(output.components.Porosity);
streamfunction = output.getStreamfunction(3000, 1); 


plotHeight = 0.7;


%perm = perm(:, 40:(end-5));
%perm = log10(perm);

probDomain = output.problemDomain;
dx = probDomain.dxCoarse;
numx = (probDomain.domainExtent.hi_i - probDomain.domainExtent.lo_i) + 1;
numy = (probDomain.domainExtent.hi_j - probDomain.domainExtent.lo_j) + 1;
x = linspace(double(probDomain.domainExtent.lo_i), double(probDomain.domainExtent.hi_i), numx)*dx;
y = linspace(double(probDomain.domainExtent.lo_j), double(probDomain.domainExtent.hi_j), numy)*dx;

x = x - min(x);
x = x*1.02;
y = y*1.02;

xlo = double(probDomain.domainExtent.lo_i)*dx;
xhi = double(probDomain.domainExtent.hi_i)*dx;
width = xhi-xlo;

h = figure();
set(h, 'Position', [200 200 400 900]);
%colormap bone;

%colormap(makeColorMap([0.0039    0.0980    0.3216], [1 1 1]));
colormap(makeColorMap( [1 1 1], [0.0039    0.0980    0.3216]));
colormap(flipud(parula));

hold on;
hcolor = pcolor(x, y, perm.');
set(hcolor,'edgecolor','none');


if doColorbar
c = colorbar();
c.Label.String = '\chi';
c.Label.Rotation = 0;
c.Label.Position = [2.2 0.6];
c.Position = [0.8 0.22 0.06 0.5];
c.Ticks = [0.01 0.99];
c.TickLabels = {'0', '1'};
c.TickLength = 0;
end

 [X, Y] = meshgrid(x, y);
   offset = 0.002; numContours = 20; maxVal = 0.2;
    contourValsNeg = linspace(-maxVal, -offset, numContours);
    contourValsPos = linspace(offset, maxVal, numContours);
    
 [CconNeg, hconNeg] = contour(X, Y, streamfunction, contourValsNeg);
 [CconPos, hconPos] = contour(X, Y, streamfunction, contourValsPos);
 
  hconNeg.LineColor = 'k';
    hconNeg.LineWidth = 2;
    hconPos.LineColor = 'k';
    hconPos.LineWidth = 2;

hold off;

set(gca,'dataAspectRatio',[1 1 1])

ax = gca;

set(ax,'clim',[0 1]);
 set(gca, 'Layer', 'top'); % put box border above imagesc
 box on;

%ax.YTick = [0,0.05, 0.1, 0.15, 0.195];
%ax.YTickLabel ={'0','', '0.1','', '0.2'};
%ax.XTick = [xlo, xhi];
%ax.XTickLabel ={'0.05', ''};



xlabel('$x$', 'Interpreter','latex');
ylabel('$y$', 'Interpreter','latex');
% 
% set(get(gca,'YLabel'),'Rotation',0, ...
%    'Units', 'Normalized', 'Position', [-1.5, 0.5, 0])

%set(gca, 'position', [0.4 0.22 0.4  plotHeight]);


set(gca,'TickLength',[0 0])

h =gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%print(h,[plot_prefix, '.pdf'],'-dpdf','-r100')
filename = plot_prefix;

if doColorbar
    filename = [filename, '-colorbar'];
end
print(h,[filename, '.png'],'-dpng','-r300')

