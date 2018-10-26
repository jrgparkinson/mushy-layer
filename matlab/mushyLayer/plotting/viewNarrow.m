%clear all;
close all;

set(0, 'defaultlinelinewidth',2);
set(0, 'defaultaxeslinewidth',2);
set(0, 'defaultpatchlinewidth',2);
set(0, 'defaultAxesFontSize',16);
% Plot upper and lower branches

getData = true;

%baseDir = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/';
%baseDir = getDataDir('channelSpacing');

dataFolder = getDataDir('optimalStates-highRes-new/'); %'/media/parkinsonjl/FREECOM HDD/';


%Process series of files:
files(1).CR = '1.04'; files(1).Ra = '250'; files(1).Nx = '84';
%files(1).CR = '1.04'; files(1).Ra = '600'; files(1).Nx = '40';
%files(1).CR = '3.0'; files(1).Ra = '35'; files(1).Nx = '168';
%files(1).CR = '3.0'; files(1).Ra = '200'; files(1).Nx = '76';

for f_i = 1:length(files)
    CR = files(f_i).CR;
    Ra = files(f_i).Ra;
    Nx = files(f_i).Nx;
    
    plot_prefix = ['CR',CR,'RaC',Ra,'Le200ChiCubedPermeabilitypts',Nx,'-'];
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
    axisExtent = [ 0.18 0.12 0.7 0.8];
    
    drawColorbar = false;
    toffset = 0;
    drawTime = false;
    drawSl = true;
    maxPsi = 0.05;
    Npsi = 10;
    mushyZoom = true;
    
    axisLabelPrecision = 2;
    setOriginToZero = [true true];
    
    [xzoom, yzoom] = plotWideDomain(output, axisExtent, drawColorbar, drawLabels, toffset, ...
        drawTime, drawSl, maxPsi, Npsi, mushyZoom, axisLabelPrecision, setOriginToZero);
    
    
    figureWidth = 300;
    figureHeight = round(0.8*figureWidth*length(yzoom)/length(xzoom));
    set(h, 'Position', [300 400 figureWidth figureHeight]);
    
    h =gcf;
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(h,[output_dir, '/',plot_prefix,'t',num2str(output.t),'.pdf'],'-dpdf','-r0')
    
end