function viewMushyOutputWideDomain
%clear all;
close all;

set(0, 'defaultlinelinewidth',2);
set(0, 'defaultaxeslinewidth',2);
set(0, 'defaultpatchlinewidth',2);
set(0, 'defaultAxesFontSize',22);
% Plot upper and lower branches

getData = true;
saveImage = true;
plotPsi = false;
plotSl = false;

%baseDir = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/';
%baseDir = getDataDir('channelSpacing');

%dataFolder = getDataDir('directionalSolidification/'); %'/media/parkinsonjl/FREECOM HDD/';
%plot_prefix = 'fixedChill-periodic-CR10.000RaC50Le100KozenyPermeabilityDa1.0e-02R1.0e-03pts256-St5-domWidth12.0-np4-';
%frames = [1000];
%output_dir = fullfile(dataFolder, 'CR10.000RaC50Le100KozenyPermeabilityDa1.0e-02R1.0e-03pts256-St5-domWidth12.0-np4-0');

% Higher res:
dataFolder = getDataDir('channelSpacingV2/'); %'/media/parkinsonjl/FREECOM HDD/';
plot_prefix = 'plt';
frames = [100000];
output_dir = fullfile(dataFolder, 'CR10.000RaC200Le100KozenyPermeabilityDa1.0e-03R1.0e-03pts512-width8.0-aspect8.0-0');


for f_i = 1:length(frames)
    
frame = frames(f_i);
dim = 2; subcycled = true;

if getData
output = MushyLayerOutput(dim, frame, output_dir, plot_prefix, subcycled);

porosity = output.dataForComp(output.components.Porosity);
streamfunction = output.getStreamfunction(100, 1);

xvel = output.dataForComp(output.components.xAdvectionvelocity);
yvel = output.dataForComp(output.components.yAdvectionvelocity);

Sl = output.dataForComp(output.components.Liquidconcentration);

probDomain = output.problemDomain;
dx = probDomain.dxCoarse;

numx = (probDomain.domainExtent.hi_i - probDomain.domainExtent.lo_i) + 1;
numy = (probDomain.domainExtent.hi_j - probDomain.domainExtent.lo_j) + 1;

xlo = double(probDomain.domainExtent.lo_i)*dx;
xhi = double(probDomain.domainExtent.hi_i)*dx;
width = xhi-xlo;

ylo = double(probDomain.domainExtent.lo_j)*dx;
yhi = double(probDomain.domainExtent.hi_j)*dx;

% if limits are almost an integer, round them up
xlo = correctLimit(xlo, dx);
xhi = correctLimit(xhi, dx);
ylo = correctLimit(ylo, dx);
yhi = correctLimit(yhi, dx);

%x = linspace(double(probDomain.domainExtent.lo_i), double(probDomain.domainExtent.hi_i), numx)*dx;
%y = linspace(double(probDomain.domainExtent.lo_j), double(probDomain.domainExtent.hi_j), numy)*dx;
% actually, rescale
x = linspace(double((xlo)), double((xhi)), numx);
y = linspace(double((ylo)), double((yhi)), numy);



%zoom_lo_z = zoom_lo_j*dx;
%zoom_hi_z = zoom_hi_j*dx;
streamfunction = streamfunction.';

chi = porosity;
U = xvel.';
V = yvel.';
%psi = streamfunction;
xzoom = x;
yzoom = y;

[X, Y] = meshgrid(xzoom, yzoom);



end

aspectRatio = (max(xzoom)-min(xzoom))/(max(yzoom)-min(yzoom));

figureWidth = 1400;
figureHeight = round(figureWidth/(0.85*aspectRatio));


contourValsSl = linspace(-0.9, 0, 7);

h = figure();
set(h, 'Position', [200 200 figureWidth figureHeight]);



%axisExtent = [0.02 0.02 0.96 0.96];
axisExtent = [0.01 0.01 0.9 0.97];

axPorosity = axes;
% This should be the 'blues' colormap from visit
%colormap(axPorosity, makeColorMap([8 48 107]/255, [104, 172, 212]/255, [1 1 1]));
colormap(axPorosity, flipud(parula));
caxis(axPorosity, [0 1]);
pcolor(xzoom, yzoom, chi.');
set(axPorosity,'dataAspectRatio',[1 1 1])
set(axPorosity,'Ydir','normal');


axPorosity.XTick = [];
axPorosity.YTick = [];
%axPorosity.YTick = linspace(round(ylo), round(yhi);

%xlabel('$x$'); ylabel('$z$');


c = colorbar();
c.Label.String='Porosity, \chi';

set(axPorosity, 'position', axisExtent);
set(axPorosity,'TickLength',[0 0])

box on;


if plotSl

axSl = axes;

%colormap(axSl, makeColorMap([0 0 0] , [230, 175, 0]/255, [200 0 0]/255));
salinityColormap = makeColorMap([0 0 0] ,  [200 0 0]/255);
salinityColormap = makeColorMap([1 1 1] ,  [0.5 0.5 0.5], [200 0 0]/255);
salinityColormap = makeColorMap([200 0 0]/255, [200 0 0]/255);

salinityColormap = makeColorMap([1 1 1] ,  [1 1 1], [1 1 1]);

colormap(axSl, salinityColormap);
caxis(axSl, [-1 0]);


[CSl, hSl] = contour(X, Y, Sl.', contourValsSl);


hSl.LineWidth = 1.0;

set(axSl,'dataAspectRatio',[1 1 1])
set(axSl, 'position', axisExtent);
set(axSl,'TickLength',[0 0])

axis(axSl, 'off');
axSl.Visible = 'off';

linkaxes([axPorosity axSl]);

end

if plotPsi
axPsi = axes;
maxPsi = 0.2; %max(max(abs(streamfunction)));
psiVals = linspace(-maxPsi,maxPsi,6);
psiColormap = makeColorMap([1 1 1], [1 1 1]);
colormap(axPsi, psiColormap);
%[cPsi,hPsi] = contour(X,Y,streamfunction.', psiVals);
%hPsi.LineWidth = 2.0;
contourValsNeg = linspace(-maxPsi,-maxPsi/5,3);
contourValsPos = linspace(maxPsi/5,maxPsi,3);
hold on;
[CconNeg, hconNeg] = contour(X, Y, streamfunction.', contourValsNeg);
[CconPos, hconPos] = contour(X, Y, streamfunction.', contourValsPos);

hconNeg.LineWidth = 2.0;
hconPos.LineWidth = 2.0;
hconPos.LineStyle = '--';
%Make arrows!
%plotArrows(CconNeg, hconNeg, 'w', 0.6);
%plotArrows(CconPos, hconPos, 'w', 0.6);
hold off;

%hPsi.LineStyle = '--';
set(axPsi,'dataAspectRatio',[1 1 1])
set(axPsi, 'position', axisExtent);
set(axPsi,'TickLength',[0 0])

axis(axPsi, 'off');
axPsi.Visible = 'off';

linkaxes([axPorosity axPsi]);
end

%axis([xzoom(1)+dx xzoom(end)-dx yzoom(1)+dx yzoom(end)-dx]);
%ax = gca;
%set(gca,'linewidth',2); %make box border bold to hide chimney outline
%set(ax,'clim',[0 1]);
%set(gca, 'Layer', 'top'); % put box border above imagesc



if saveImage
     h =gcf;
     set(h,'Units','Inches');
     pos = get(h,'Position');
     set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
     saveFile = fullfile(output_dir, ['t',num2str(output.t)]);
     
     %print(h,[saveFile, '.svg'],'-dsvg','-r50')
     print(h,[saveFile, '.png'],'-dpng','-r250')
     
     fprintf('Saved to %s \n', saveFile);
     
end

end

end

function output = correctLimit(lim, dx)

output = lim;

if abs(round(lim)-lim) < 2*dx
    output = round(output);
end


end
