%clear all;
close all;

set(0, 'defaultlinelinewidth',2);
set(0, 'defaultaxeslinewidth',2);
set(0, 'defaultpatchlinewidth',2);
set(0, 'defaultAxesFontSize',22);
% Plot upper and lower branches

% Length scale (metres)
L = 2;

% /home/parkinsonjl/mnt/sharedStorage/AW002_PARKINSON_MUSH/optimalStates-highRes-newSoluteFlux/
% CR1.15RaC200Le200ChiCubedPermeabilitypts64-0/
%mushyLayerLowC-upperBranch-insulating-CR1.15RaC200Le200ChiCubedPermeabilitypts64-341843.2d.hdf5

dataFolder = '/home/parkinsonjl/mnt/sharedStorage/AW002_PARKINSON_MUSH/';


plot_prefix = '/mushyLayerLowC-upperBranch-insulating-CR1.15RaC200Le200ChiCubedPermeabilitypts64-';

output_dir = [ dataFolder, ...
    '/optimalStates-highRes-newSoluteFlux/CR1.15RaC200Le200ChiCubedPermeabilitypts64-0/'];

frame = 341843;
dim = 2; subcycled = true;

output = MushyLayerOutput(dim, frame, output_dir, plot_prefix, subcycled);

perm = output.dataForComp(output.components.Permeability);


plotHeight = 0.7;


%perm = perm(:, 40:(end-5));
%perm = log10(perm);

probDomain = output.problemDomain;
dx = probDomain.dxCoarse;
numx = (probDomain.domainExtent.hi_i - probDomain.domainExtent.lo_i) + 1;
numy = (probDomain.domainExtent.hi_j - probDomain.domainExtent.lo_j) + 1;
x = linspace(double(probDomain.domainExtent.lo_i), double(probDomain.domainExtent.hi_i), numx)*dx;
y = linspace(double(probDomain.domainExtent.lo_j), double(probDomain.domainExtent.hi_j), numy)*dx;


xlo = double(probDomain.domainExtent.lo_i)*dx;
xhi = double(probDomain.domainExtent.hi_i)*dx;
width = xhi-xlo;



zoom_hi_i = 8;
zoom_hi_x = 8*dx;
zoom_lo_j = 44;
zoom_hi_j = 76;
zoom_lo_z = zoom_lo_j*dx;
zoom_hi_z = zoom_hi_j*dx;
chimneyZoomPerm = perm((end-zoom_hi_i):end, zoom_lo_j:zoom_hi_j);
xzoom = x((end-zoom_hi_i):end);
yzoom = y(zoom_lo_j: zoom_hi_j);

h = figure();
set(h, 'Position', [200 200 600 400]);
%colormap bone;
colormap(makeColorMap([0.0039    0.0980    0.3216], [1 1 1]));

subplot(1, 2, 1)
imagesc(xzoom, yzoom, chimneyZoomPerm.');

ax = gca;
ax.YTick = [zoom_lo_z, zoom_hi_z*0.99];

n = -2;

disp_lo_z = round(zoom_lo_z*10^(-n))/(10^(-n));
disp_hi_z = round(zoom_hi_z*10^(-n))/(10^(-n));
disp_hi_x = round(zoom_hi_x*10^(-n))/(10^(-n));

ax.YTickLabel ={num2str(disp_lo_z),num2str(disp_hi_z)};
ax.XTick = [min(xzoom), max(xzoom)];
ax.XTickLabel ={num2str(disp_hi_x), '0'};

set(ax,'dataAspectRatio',[1 1 1])
set(ax,'Ydir','normal');
set(ax,'Xdir','reverse');
%set(gca, 'position', [-0.05 0.5 0.9 0.73]);
ylabel('$z$', 'Interpreter','latex');
xlabel('$x$', 'Interpreter','latex');

set(ax,'TickLength',[0 0])



set(get(gca,'YLabel'),'Rotation',0, ...
   'Units', 'Normalized', 'Position', [-0.8, 0.5, 0])

set(gca, 'position', [0.02 0.22 0.4 plotHeight]);



subplot(1, 2, 2)

imagesc(x, y, perm.');
set(gca,'Ydir','normal');
set(gca,'Xdir','reverse');
c = colorbar();
c.Label.String = '\Pi';
c.Label.Rotation = 0;
c.Label.Position = [2.0 0.6];
c.Position = [0.8 0.22 0.05 0.7];
c.Ticks = [0.01 0.99];
c.TickLabels = {'0', '1'};
c.TickLength = 0;

set(gca,'dataAspectRatio',[1 1 1])

ax = gca;
ax.YTick = [0,0.05, 0.1, 0.15, 0.195];
ax.YTickLabel ={'0','', '0.1','', '0.2'};
ax.XTick = [xlo, xhi];
ax.XTickLabel ={'0.05', ''};


xlabel('$x$', 'Interpreter','latex');

%set(get(gca,'YLabel'),'Rotation',0, ...
%    'Units', 'Normalized', 'Position', [-1.5, 0.5, 0])

set(gca, 'position', [0.4 0.22 0.4  plotHeight]);


set(gca,'TickLength',[0 0])

h =gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'snapshot.pdf','-dpdf','-r0')

