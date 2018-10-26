%clear all;
close all;

set(0, 'defaultlinelinewidth',2);
set(0, 'defaultaxeslinewidth',2);
set(0, 'defaultpatchlinewidth',2);
set(0, 'defaultAxesFontSize',18);
tic;

% Plot upper and lower branches

getData = true;

%baseDir = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/';
%baseDir = getDataDir('channelSpacing');

dataFolder = getDataDir('middleton/'); %'/media/parkinsonjl/FREECOM HDD/';

%plot_prefix = 'CR1.25RaC200Le200ChiCubedPermeabilitypts2048-';
%output_dir = [dataFolder, 'CR1.25RaC200Le200ChiCubedPermeabilitypts2048-0/'];
%frames = [0 4400 5850 8700 16150 27000 36000, 38000];


plot_prefix = 'plt';
output_dir = [dataFolder, 'CR1.179RaC800Le200KozenyPermeabilitypts128-S3.5-TB-20.0-R0.013-domWidth1.0-0/'];
frames = [5475];

output_dir = [dataFolder, 'CR1.179RaC80Le200KozenyPermeabilitypts256-S3.5-TB-20.0-R0.13-0/'];

out_dir = {[dataFolder, 'CR1.179RaC80Le200KozenyPermeabilitypts256-S3.5-TB-10.0-R0.13-0/'], ...
    [dataFolder, 'CR1.179RaC80Le200KozenyPermeabilitypts256-S3.5-TB-15.0-R0.13-0/'], ...
    [dataFolder, 'CR1.179RaC80Le200KozenyPermeabilitypts256-S3.5-TB-20.0-R0.13-0/']    };

t = 0.13;
frames = [1450 1450 1500];


depths = [0.5];
totalDepth = 3.0;

dim = 2; subcycled = true;

figureAspectRatio = 1.0;

indFigHeight = 0.95;
figAxisHeight = indFigHeight*0.95;
indFigWidth = 0.25;


plotHeight = 200;
plotWidth = 1200;

h = figure();
set(h, 'Position', [100 300 plotWidth plotHeight]);



        toffset = 0;
    
    
    
    %axisExtent = [0.02 0.02 0.96 0.96];
    
    % Plot three axes
    for axi = 1:3
        
        % frame = frames(f_i);
    frame = frames(axi);
    output_dir = out_dir{axi};
    
    if getData
        output = MushyLayerOutput(dim, frame, output_dir, plot_prefix, subcycled);
 porosity = output.dataForComp(output.components.Porosity);
         Sl = output.dataForComp(output.components.Liquidconcentration);
         [X,Z] = output.grid();
         t = output.t;
         save([output_dir, '/toPlot',num2str(frame),'.mat'], 'porosity', 'Sl','X','Z','t');
    end
    
    
    axisExtent = [0.06+indFigWidth*1.1*(axi-1) 0.06 indFigWidth figAxisHeight];
    
    colorbar = false;
    %axesLabels = (f_i==length(frames));
    if axi == 1
        axesLabels  = [false true];
    else
        axesLabels = [false false];
    end
    
    drawTime = false;
    drawSl = false;
    maxPsi = 15; %0.35;
    Npsi = 5;
    
   % plotWideDomainEGU(output, axisExtent, colorbar, axesLabels,toffset, ...
   %     drawTime, drawSl, maxPsi, Npsi, true);
    
    %t = toc;
    %fprintf('t=%1.5f, \n', t);
    fprintf('t=%1.5f \n', output.t);
    
    end
    


 h =gcf;
  set(h,'Units','Inches');
  pos = get(h,'Position');
  set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
  filename = ['/home/parkinsonjl/convection-in-sea-ice/figures', '/timeSeries-t',num2str(t)];
  filename = [dataFolder, '/timeSeries-t',num2str(t)];
  print(h,[filename, '.pdf'],'-dpdf','-r0')
  print(h,[filename, '.png'],'-dpng','-r400')
