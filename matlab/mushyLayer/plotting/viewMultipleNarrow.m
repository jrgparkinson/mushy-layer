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
files(1).CR = '1.04'; files(1).Ra = '250'; files(1).Nx = '84'; files(1).title = 'a';
files(2).CR = '1.04'; files(2).Ra = '600'; files(2).Nx = '40'; files(2).title = 'b';
files(3).CR = '3.0'; files(3).Ra = '35'; files(3).Nx = '168'; files(3).title = 'c';
files(4).CR = '3.0'; files(4).Ra = '200'; files(4).Nx = '76'; files(4).title = 'd';
plotfilename = 'CRRaExampleSol2';

h = figure();
figureWidth = 900;
figureHeight = 1000;
set(h, 'Position', [100 100 figureWidth figureHeight]);
  
    
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
    
    
    
   
    %drawLabels = false;
    %axisExtent = [0.02 0.02 0.96 0.96];
    
    drawLabels = true;
    
    %axisExtent = [ 0.18 0.12 0.7 0.8];
    axisExtent(3) = 0.4;
    axisExtent(4) = 0.38;
    
    % Layout is
    %  1   2
    %  3   4
    
    if f_i == 1 || f_i == 2
        axisExtent(2) = 0.58;
    else
        axisExtent(2) = 0.1;
    end
    
    if f_i == 1 || f_i == 3
       axisExtent(1) = 0.1; 
    else
        axisExtent(1) = 0.6;
    end
    
    if f_i == 2
        drawColorbar = true;
    else
        drawColorbar = false;
    end
    
    
    toffset = 0;
    drawTime = false;
    drawSl = true;
    maxPsi = 0.45;
    Npsi = 20;
    mushyZoom = true;
    
    axisLabelPrecision = 2;
    setOriginToZero = [true true];
    
    [xzoom, yzoom] = plotWideDomain(output, axisExtent, drawColorbar, drawLabels, toffset, ...
        drawTime, drawSl, maxPsi, Npsi, mushyZoom, axisLabelPrecision, setOriginToZero);
    
    title(['(', files(f_i).title, ') $\mathcal{C}=',num2str(str2num(CR)-1),', Rm_S=',Ra,'$']);
   
end


   
    h =gcf;
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    saveTo = ['/home/parkinsonjl/convection-in-sea-ice/MushyLayer/matlab/ConfirmationOfStatus/',plotfilename,'.pdf'];
    print(h,saveTo,'-dpdf','-r0')
    
    fprintf('Saved to %s \n', saveTo);