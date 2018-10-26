clear all;
close all;

set(0, 'defaultlinelinewidth',2);
set(0, 'defaultaxeslinewidth',2);
set(0, 'defaultpatchlinewidth',2);
set(0, 'defaultAxesFontSize',16);
% Plot upper and lower branches

getData = true;
savePlot = true;
%baseDir = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/';
%baseDir = getDataDir('channelSpacing');

dataFolder = getDataDir('optimalStates-highRes-new/'); %'/media/parkinsonjl/FREECOM HDD/';


%Process series of files:
files(1).CR = '1.15'; files(1).Ra = '400'; files(1).Nx = '64';
permeability = 'ChiCubed';
%files(1).CR = '1.04'; files(1).Ra = '600'; files(1).Nx = '40';
%files(1).CR = '3.0'; files(1).Ra = '35'; files(1).Nx = '168';
%files(1).CR = '3.0'; files(1).Ra = '200'; files(1).Nx = '76';
extras = '';
axisExtent = [ 0.18 0.1 0.54 0.83];
figureWidth = 300;
aspectRatio = 0.8;

%dataFolder = getDataDir('optimalStates-Brinkman/'); %'/media/parkinsonjl/FREECOM HDD/';
%files(1).CR = '1.15'; files(1).Ra = '30000'; files(1).Nx = '128';
%permeability = 'Kozeny';
%extras = 'Da5e-3R1e-10';
%axisExtent = [ 0.07 0.06 0.83 0.9];
%figureWidth = 700; aspectRatio = 1.0;
for f_i = 1:length(files)
    CR = files(f_i).CR;
    Ra = files(f_i).Ra;
    Nx = files(f_i).Nx;
    
    plot_prefix = ['CR',CR,'RaC',Ra,'Le200',permeability,'Permeability',extras,'pts',Nx,'-'];
    %CR1.15RaC400Le200ChiCubedPermeabilitypts64-steady
    
    output_dir = [dataFolder, plot_prefix, 'steady/'];
    %frame = 064200;
    %dim = 2; subcycled = true;
    
    if getData
        %output = MushyLayerOutput(dim, frame, output_dir, plot_prefix, subcycled);
        output = getFinalPlotFile(output_dir, plot_prefix);
    end
    
    
    
    h = figure();
    
    %drawLabels = false;
    %axisExtent = [0.02 0.02 0.96 0.96];
    
    drawLabels = true;
    
    
    drawColorbar = true;
    toffset = 0;
    drawTime = false;
    drawSl = true;
    maxPsi = 0.15; %0.4;
    Npsi = 10;
    mushyZoom = true;
    
    axisLabelPrecision = 2;
    setOriginToZero = [true true];
        
    [xzoom, yzoom] = plotSingleChannelReflected(output, axisExtent, drawColorbar, drawLabels, toffset, ...
        drawTime, drawSl, maxPsi, Npsi, mushyZoom, axisLabelPrecision, setOriginToZero);
    
    
    
    figureHeight = round(aspectRatio*figureWidth*length(yzoom)/length(xzoom))
    set(h, 'Position', [300 400 figureWidth figureHeight]);
    
    if savePlot
    h =gcf;
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    saveTo = [dataFolder, '/',plot_prefix,'t',num2str(output.t)];
    saveTo = ['/home/parkinsonjl/convection-in-sea-ice/figures/',plot_prefix,'t',num2str(output.t)];
    fprintf('Save to %s \n', saveTo);
    print(h,[saveTo,'.pdf'],'-dpdf','-r20')
     print(h,[saveTo,'.png'],'-dpng','-r500')
    end
end