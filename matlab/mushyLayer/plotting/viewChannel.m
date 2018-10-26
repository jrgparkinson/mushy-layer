clear all;
close all;

set(0, 'defaultlinelinewidth',2);
set(0, 'defaultaxeslinewidth',2);
set(0, 'defaultpatchlinewidth',2);
set(0, 'defaultAxesFontSize',16);
% Plot upper and lower branches

getData = true;

%baseDir = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/';
%baseDir = getDataDir('channelSpacing');

dataFolder = getDataDir('saltDiffusionChannelWidth/'); %'/media/parkinsonjl/FREECOM HDD/';


%Process series of files:
%Da = 5e-3;
% Ra = '4000'; permeability = 'KozenyPermeability';
% Conc = '6.0';
% files(1).CR = Conc; files(1).Ra = Ra; files(1).Le = '1000.0'; files(1).perm=permeability; files(1).Nx = '176';
% files(2).CR = Conc; files(2).Ra = Ra; files(2).Le = '100.0'; files(2).perm=permeability; files(2).Nx = '176';
% files(3).CR = Conc; files(3).Ra = Ra; files(3).Le = '10.0'; files(3).perm=permeability; files(3).Nx = '176';
% plotName = 'RaC4e3Le1000_100_10';
% 
% 
% Lewis = '500.0';
% files(1).CR = Conc; files(1).Ra = '2000.0'; files(1).Le = Lewis; files(1).perm=permeability; files(1).Nx = '176';
% files(2).CR = Conc; files(2).Ra = '4000'; files(2).Le = Lewis; files(2).perm=permeability; files(2).Nx = '176';
% files(3).CR = Conc; files(3).Ra = '10000.0'; files(3).Le = Lewis; files(3).perm=permeability; files(3).Nx = '176';
% plotName = 'RaC2e3_4e3_10e3Le500';


Ra = '4000'; permeability = 'KozenyPermeability';
Conc = '6.0';
files(1).CR = Conc; files(1).Ra = Ra; files(1).Le = '1000.0'; files(1).perm=permeability; files(1).Nx = '176';
files(2).CR = Conc; files(2).Ra = Ra; files(2).Le = '100.0'; files(2).perm=permeability; files(2).Nx = '176';
files(3).CR = Conc; files(3).Ra = Ra; files(3).Le = '10.0'; files(3).perm=permeability; files(3).Nx = '176';

Lewis = '500.0';
files(4).CR = Conc; files(4).Ra = '2000.0'; files(4).Le = Lewis; files(4).perm=permeability; files(4).Nx = '176';
files(5).CR = Conc; files(5).Ra = '4000'; files(5).Le = Lewis; files(5).perm=permeability; files(5).Nx = '176';
files(6).CR = Conc; files(6).Ra = '10000.0'; files(6).Le = Lewis; files(6).perm=permeability; files(6).Nx = '176';

plotName = 'All';


% Use same lims for all plots
lims = [-0.1 0.12; -0.4 0.65];


filesWide = 3;
filesHigh = 2;

axisXOffset = 0.0;

plotWidth = (1-axisXOffset)/filesWide;
plotHeight = (1-0.35)/filesHigh;

h = figure();

figureWidth = 200*filesWide;
figureHeight = 400*filesHigh;
set(h, 'Position', [100 100 figureWidth figureHeight]);

titles = {'a','b','c','d','e','f'};

for f_i = 1:length(files)
    CR = files(f_i).CR;
    Ra = files(f_i).Ra;
    Nx = files(f_i).Nx;
    Le = files(f_i).Le;
    perm = files(f_i).perm;
    letterCaption = titles{f_i};
    
    plot_prefix = ['CR',CR,'RaC',Ra,'Le',Le,perm,'pts',Nx,'-'];
    output_dir = [dataFolder, plot_prefix, 'steady/'];
    
    if getData
        %output = MushyLayerOutput(dim, frame, output_dir, plot_prefix, subcycled);
        
        output = getFinalPlotFile(output_dir, plot_prefix);
        
        if length(output.levelArray) == 0
            output_dir = [dataFolder, plot_prefix, '0/'];
            output = getFinalPlotFile(output_dir, plot_prefix);
        end
    end
    
    %drawLabels = false;
    %axisExtent = [0.02 0.02 0.96 0.96];
    
    drawLabels = true;
    
    if f_i > 3
        bottom = 0.08;
    else
        bottom = 0.6;
    end
    axisExtent = [axisXOffset + plotWidth*(mod(f_i-1,filesWide)) bottom plotWidth plotHeight];
    
    drawColorbar = false;
    toffset = 0;
    drawTime = false;
    drawSl = true;
    maxPsi = 6.0;
    Npsi = 15;
    mushyZoom = true;
    
    axisLabelPrecision = 2;
    setOriginToZero = [true true];
    
    
    [xzoom, yzoom] = plotChannel(output, axisExtent, drawColorbar, drawLabels, toffset, ...
        drawTime, drawSl, maxPsi, Npsi, mushyZoom, axisLabelPrecision, setOriginToZero, lims);
    
    titleLine1 = sprintf('(%s) $Le = %d$,', letterCaption, str2num(Le));
    titleLine2 = sprintf('$Ra_s = %d$', str2num(Ra));
    title({titleLine1, titleLine2});
    
end

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
outputFile = [dataFolder, '/channelFields',plotName,'.pdf'];
print(h,outputFile,'-dpdf','-r0')
fprintf('Saved as %s', outputFile);