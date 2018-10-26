%clear all;
function viewWideDomainTimeSeriesEGU

close all;

set(0, 'defaultlinelinewidth',2);
set(0, 'defaultaxeslinewidth',2);
set(0, 'defaultpatchlinewidth',2);
fontsize = 14;
set(0, 'defaultAxesFontSize',fontsize);

savePlots = false;

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

dataFolder = '/home/parkinsonjl/mnt/raymaster/convection-in-sea-ice/MushyLayer/matlab/EGUPoster/';
dataFolder = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/matlab/EGUPoster/';

Temperature = [-10 -15 -20];
out_dir = {[dataFolder, 'T-10/'], ...
    [dataFolder, 'T-15/'], ...
    [dataFolder, 'T-20/']    };

% frames = [1400 1400 1400;
%     2000 2000 2000;
%     11000 11000 114500];

% minutes: 2, 15, 40, 50,   80 (OR LARGEST possible)

frames = [2000 2000 2000; %2 mins
    %11000 11000 8000; % 12 mins 
    %11000 11000 11000; % 15 mins (need to change -15 and -20)
   % 15000 15000 15000; % 23 mins
   % 29000 29000 29500; % 47 mins
     12500 13000 13000; % 20 mins - previously used this
     19000 19000 19300; % 30 mins
   % 25000 25000 26000; % 40 mins
    32000 31500 32000; % 50 mins
    %42500 36000 37000; % 58 mins
     44000 37500 38000; % 60 mins
   % 52000 43500 175000; %70 mins5
   % 72000 58000 196000 %rubbish
    %%50000 50000 50000;
    %75000 156500 114500
    ];


dim = 2; subcycled = true;


plotHeight = 900;
plotWidth = 550;

h = figure();
set(h, 'Position', [100 100 plotWidth plotHeight]);

% First do the diffusive growth plot + BCs
w =  0.25; h=0.15;

Slmax = 0.085; %55g/kg


depths = [1.0 1.5 2.5 2.5 2.5];

baseX = 0.16;

yPos1 = 0.84;
axExtent(1,1,:) = [baseX         yPos1 w h*depths(1)];
axExtent(1,2,:) = [baseX+w*1.1   yPos1 w h*depths(1)];
axExtent(1,3,:) = [baseX+2*w*1.1 yPos1 w h*depths(1)];

yPos2 = 0.675;
axExtent(2,1,:) = [baseX         yPos2 w h*depths(2)];
axExtent(2,2,:) = [baseX+w*1.1   yPos2 w h*depths(2)];
axExtent(2,3,:) = [baseX+2*w*1.1 yPos2 w h*depths(2)];

yPos3 = 0.4;
axExtent(3,1,:) = [baseX         yPos3 w h*depths(3)];
axExtent(3,2,:) = [baseX+w*1.1   yPos3 w h*depths(3)];
axExtent(3,3,:) = [baseX+2*w*1.1 yPos3 w h*depths(3)];

yPos4 = 0.17;
axExtent(4,1,:) = [baseX         yPos4 w h*depths(4)];
axExtent(4,2,:) = [baseX+w*1.1   yPos4 w h*depths(4)];
axExtent(4,3,:) = [baseX+2*w*1.1 yPos4 w h*depths(4)];
    
yPos5 = -0.06;
axExtent(5,1,:) = [baseX         yPos5 w h*depths(5)];
axExtent(5,2,:) = [baseX+w*1.1   yPos5 w h*depths(5)];
axExtent(5,3,:) = [baseX+2*w*1.1 yPos5 w h*depths(5)];







for time = 1:length(depths)
    
    
    depth = depths(time);
    
    %axisExtent = [0.02 0.02 0.96 0.96];
    
    % Plot three axes
    for Ti = 1:length(Temperature)
        
        % frame = frames(f_i);
        frame = frames(time, Ti);
        output_dir = out_dir{Ti};
        
        
        
        
        %    colorbar = (f_i==2);
        %axesLabels = (f_i==length(frames));
        
        
        dataFile = [output_dir, 'toPlot', num2str(frame),'.mat'];
        if Ti==1
            doLabel = true;
        else
            doLabel = false;
        end
        
        
        
        axisExtent = axExtent(time,Ti,:);
        
        axPorosity = makeHeleShawPlot(dataFile, depth, doLabel, axisExtent, Slmax);
        
        if time == 1
            
           title(['$T_a=',num2str(Temperature(Ti)),'^\circ$ C'], 'FontSize', fontsize); 
        end
        
        %thisAx.Position = axisExtent;
        
    end
    
end


if savePlots

 h =gcf;
 h.InvertHardcopy = 'off';
 h.Color = 'white';
  set(h,'Units','Inches');
  pos = get(h,'Position');
  set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
  filename = ['/home/parkinsonjl/convection-in-sea-ice/figures', '/timeSeries-t'];
  print(h,[filename, '.pdf'],'-dpdf','-r0')
  print(h,[filename, '.png'],'-dpng','-r400')
  
end

end

