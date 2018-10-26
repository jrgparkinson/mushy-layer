%clear all;
close all;

set(0, 'defaultlinelinewidth',2);
set(0, 'defaultaxeslinewidth',2);
set(0, 'defaultpatchlinewidth',2);
set(0, 'defaultAxesFontSize',18);
% Plot upper and lower branches

getData = true;

%baseDir = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/';
%baseDir = getDataDir('channelSpacing');

dataFolder = getDataDir('channelSpacing/'); %'/media/parkinsonjl/FREECOM HDD/';

%plot_prefix = 'CR1.25RaC200Le200ChiCubedPermeabilitypts2048-';
%output_dir = [dataFolder, 'CR1.25RaC200Le200ChiCubedPermeabilitypts2048-0/'];
%frames = [0 4400 5850 8700 16150 27000 36000, 38000];


plot_prefix = 'fixedChill-periodic-CR1.25RaC200Le200ChiCubedPermeabilitypts256-';
output_dir = [dataFolder, 'CR1.25RaC200Le200ChiCubedPermeabilitypts256-steady/'];
frames = [500 1000 2500 4000 8500];
depths = [1 1 1.5 2 2.5];
totalDepth = 3.2;

dim = 2; subcycled = true;

figureAspectRatio = 7.8;

indFigHeight = 0.94/(length(frames));
figAxisHeight = indFigHeight*0.95;
indFigWidth = 0.83; %Including space for the colorbar, was 0.96

figureWidth = 900;
figureHeight = round(length(frames)*figureWidth/figureAspectRatio);
plotWidth = figureWidth + 100; % Extra space for a colorbar

h = figure();
set(h, 'Position', [50 50 plotWidth figureHeight]);


for f_i =1:length(frames)
    frame = frames(f_i);
    
if getData
output = MushyLayerOutput(dim, frame, output_dir, plot_prefix, subcycled);
end

if f_i == 1
    toffset = output.t;
end



%axisExtent = [0.02 0.02 0.96 0.96];
axisExtent = [0.06 0.06+indFigHeight*(length(frames)-f_i) indFigWidth figAxisHeight];

colorbar = (f_i==2);
axesLabels = (f_i==length(frames));

 drawTime = false;
 drawSl = true;
 maxPsi = 0.35;
 Npsi = 5;
 
plotWideDomain(output, axisExtent, colorbar, axesLabels,toffset, ...
     drawTime, drawSl, maxPsi, Npsi);

end

h =gcf;
 set(h,'Units','Inches');
 pos = get(h,'Position');
 set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
 print(h,[output_dir, '/timeSeries2.pdf'],'-dpdf','-r0')
