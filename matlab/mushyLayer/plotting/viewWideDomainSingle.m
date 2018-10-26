%clear all;
close all;

set(0, 'defaultlinelinewidth',2);
set(0, 'defaultaxeslinewidth',2);
set(0, 'defaultpatchlinewidth',2);
set(0, 'defaultAxesFontSize',22);
% Plot upper and lower branches

getData = false;

%baseDir = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/';
%baseDir = getDataDir('channelSpacing');

dataFolder = getDataDir('channelSpacing/'); %'/media/parkinsonjl/FREECOM HDD/';

plot_prefix = 'CR1.25RaC250Le200ChiCubedPermeabilitypts1024-';
output_dir = [dataFolder, 'CR1.25RaC250Le200ChiCubedPermeabilitypts1024-0/'];
frame = 064200;
dim = 2; subcycled = true;

if getData
output = MushyLayerOutput(dim, frame, output_dir, plot_prefix, subcycled);
end

figureWidth = 900;
figureHeight = round(figureWidth/3.9);

h = figure();
set(h, 'Position', [200 200 figureWidth figureHeight]);

axisExtent = [0.02 0.02 0.96 0.96];

plotWideDomain(output, axisExtent);

h =gcf;
 set(h,'Units','Inches');
 pos = get(h,'Position');
 set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
 print(h,[output_dir, '/t',num2str(output.t),'.pdf'],'-dpdf','-r0')
