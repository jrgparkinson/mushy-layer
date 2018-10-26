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

dataFolderDarcy = getDataDir('optimalStates-highRes-new/'); %'/media/parkinsonjl/FREECOM HDD/';
dataFolderDarcy = getDataDir('narrowDarcySameBCsAsBrinkman/'); %'/media/parkinsonjl/FREECOM HDD/';

dataFolderBrinkman = getDataDir('optimalStates-Brinkman/'); %'/media/parkinsonjl/FREECOM HDD/';

Rm = '35';
CR = '6.0';
Da = 5e-3;

Ra = '7000';

%Process series of files:
%files(1).CR = '1.04'; files(1).Ra = '250'; files(1).Nx = '84'; files(1).title = 'a';
%files(2).CR = '1.04'; files(2).Ra = '600'; files(2).Nx = '40'; files(2).title = 'b';
%files(3).CR = '3.0'; files(3).Ra = '35'; files(3).Nx = '168'; files(3).title = 'c';
%files(4).CR = '3.0'; files(4).Ra = '200'; files(4).Nx = '76'; files(4).title = 'd';
plotfilename = 'CompareDarcyBrinkman';

h = figure();
figureWidth = 900;
figureHeight = 400;
set(h, 'Position', [100 300 figureWidth figureHeight]);
  

    plot_prefixDarcy = ['CR',CR,'RaC',Rm,'Le200ChiCubedPermeabilitypts168-'];
    output_dirDarcy = [dataFolderDarcy, 'CR6.0RaC35Le200ChiCubedPermeabilityH1.0168-0/'];
    outputDarcy = getFinalPlotFile(output_dirDarcy, plot_prefixDarcy);
    
    plot_prefixDB = ['CR',CR,'RaC',Ra,'Le200KozenyPermeabilitypts96-'];
    output_dirDB = [dataFolderBrinkman, 'CR6.0RaC7000Le200KozenyPermeabilityDa5e-3R1e-10pts168-stoppedEarly2/'];
    outputDB = getFinalPlotFile(output_dirDB);
   
   
    %drawLabels = false;
    %axisExtent = [0.02 0.02 0.96 0.96];
    
    drawLabels = true;
    
    %axisExtent = [ 0.18 0.12 0.7 0.8];
    widthDarcy = 0.35;
    widthDB = widthDarcy*1.03; % Determined empirically
    
    heightDB = 0.7;
    heightDarcy = heightDB*1.03;
    
    axisExtentDarcy = [0.08 0.83-heightDarcy widthDarcy heightDarcy];
    axisExtentDB = [0.6 0.83-heightDB widthDB heightDB];
    
    toffset = 0;
    drawTime = false;
    drawSl = true;
    maxPsi = 4.5;
    Npsi = 20;
    mushyZoom = true;
    drawColorbar = true;
    
    axisLabelPrecision = 2;
    setOriginToZero = [true true];
    
    [xzoomDarcy, yzoomDarcy] = plotWideDomain(outputDarcy, axisExtentDarcy, drawColorbar, drawLabels, toffset, ...
        drawTime, drawSl, maxPsi, Npsi, mushyZoom, axisLabelPrecision, setOriginToZero);
    
    %title(['(', files(f_i).title, ', $\mathcal{C}=',num2str(str2num(CR)-1),', Rm_S=',Ra,'$)']);
    CRprint = num2str(str2num(CR)-1);
    %title({['(a) $\mathcal{C} = ',CRprint,', Rm_S = ',Rm, '$'], ''});
    title('(a)');
    
    axDarcy = gca;
   
    drawColorbar = false;
    [xzoomDB, yzoomDB] = plotWideDomain(outputDB, axisExtentDB, drawColorbar, drawLabels, toffset, ...
        drawTime, drawSl, maxPsi, Npsi, mushyZoom, axisLabelPrecision, setOriginToZero);
    RmDB = num2str(str2num(Ra)*Da);
%title({['(b) $\mathcal{C} = ',CRprint,', Rm_S = ',RmDB, '$'],''});
 title('(b)');
 
axDB = gca;


ratio = (xzoomDB(end)-xzoomDB(1))/(xzoomDarcy(end)-xzoomDarcy(1));

%pbaspect(axDarcy, 

   
    h =gcf;
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    outputFolder = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/matlab/ConfirmationOfStatus/';
    print(h,[outputFolder,plotfilename,'.pdf'],'-dpdf','-r0')
    
    fprintf('Saved to %s \n', [outputFolder, '/',plotfilename,'.pdf']);