%clear all;
close all;

set(0, 'defaultlinelinewidth',2);
set(0, 'defaultaxeslinewidth',2);
set(0, 'defaultpatchlinewidth',2);
set(0, 'defaultAxesFontSize',22);
% Plot upper and lower branches

% Length scale (metres)
L = 2;

baseDir = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/';

dataFolder = '/media/parkinsonjl/FREECOM HDD/';

%/media/parkinsonjl/FREECOM HDD/mushyLayer-periodic/CR1.25RaC100.0Le200ChiCubedPermeabilitypts512-0v1/.2d.hdf5
output_dir = [ dataFolder, 'mushyLayer-periodic/CR1.25RaC100.0Le200ChiCubedPermeabilitypts512-0v1/'];
plot_prefix = ['mushyLayer-periodic-CR1.25RaC100.0Le200ChiCubedPermeabilitypts512-'];

frame = 026200;
dim = 2; subcycled = true;

output = MushyLayerOutput(dim, frame, output_dir, plot_prefix, subcycled);

porosity = output.dataForComp(output.components.Porosity);
%streamfunction = output.getStreamfunction(100, 1);

xvel = output.dataForComp(output.components.xAdvectionvelocity);
yvel = output.dataForComp(output.components.yAdvectionvelocity);

Sl = output.dataForComp(output.components.Liquidconcentration);

probDomain = output.problemDomain;
dx = probDomain.dxCoarse;

numx = (probDomain.domainExtent.hi_i - probDomain.domainExtent.lo_i) + 1;
numy = (probDomain.domainExtent.hi_j - probDomain.domainExtent.lo_j) + 1;

x = linspace(double(probDomain.domainExtent.lo_i), double(probDomain.domainExtent.hi_i), numx)*dx;
y = linspace(double(probDomain.domainExtent.lo_j), double(probDomain.domainExtent.hi_j), numy)*dx;


xlo = double(probDomain.domainExtent.lo_i)*dx;
xhi = double(probDomain.domainExtent.hi_i)*dx;
width = xhi-xlo;



%zoom_lo_z = zoom_lo_j*dx;
%zoom_hi_z = zoom_hi_j*dx;
%streamfunction = streamfunction.';

chi = porosity;
U = xvel.';
V = yvel.';
%psi = streamfunction;
xzoom = x;
yzoom = y;

%Make it four times wider
stretch = 2;
originalWidth = length(xzoom);
for i=1:(stretch-1)
   xzoom(end+1:end+originalWidth) = xzoom(1:originalWidth);
   chi(end+1:end+originalWidth, :) = chi(1:originalWidth, :);
end

[X, Y] = meshgrid(xzoom, yzoom);

refinement = 8;
chiRef = resizem(chi, refinement, 'bicubic');
XRef = resizem(X, refinement, 'bilinear');
YRef = resizem(Y, refinement, 'bilinear');

figureWidth = 1000;
figureHeight = round(figureWidth/(4*stretch));

% Mush liquid boundary
%[mush_x, mush_y] = find(chi<1.0);
%idx_liquid = 1-idx_mush;
%refinement = 16;
%chiRef = resizem(chi, refinement, 'bicubic');
%Xref = resizem(X, refinement, 'bilinear');
%Yref = resizem(Y, refinement, 'bilinear');
%chiRef = smoothn(chiRef, 'robust');

% Velocity
%offset = 5.0; numContours = 6; maxVal = 90;
%contourValsNeg = linspace(-maxVal, -offset, numContours);
%contourValsPos = linspace(offset, maxVal, numContours);


h = figure();
set(h, 'Position', [200 200 figureWidth figureHeight]);



% left bottom width height
axisExtent = [0.0 0.0 1.0 1.0];

axPorosity = axes;
% This should be the 'blues' colormap from visit
colormap(axPorosity, makeColorMap([8 48 107]/255, [104, 172, 212]/255, [1 1 1]));
caxis(axPorosity, [0 1]);
imagesc(XRef(1, :), YRef(:, 1), chiRef.');
set(axPorosity,'dataAspectRatio',[1 stretch 1])
set(axPorosity,'Ydir','normal');
axPorosity.XTick = [];
axPorosity.YTick = [];
set(axPorosity, 'position', axisExtent);
set(axPorosity,'TickLength',[0 0])

box off;
axis off;



h =gcf;
 set(h,'Units','Inches');
 pos = get(h,'Position');
 set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,['background.pdf'],'-dpdf','-r0');

