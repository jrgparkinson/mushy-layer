function makeSingleVideo(folder, redoFrames, makeVideo)
close all;

thisDir = strrep(mfilename('fullpath'), 'makeSingleVideo', '');

%local_base_dir = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/matlab/EGUPoster/';
%remote_base = '/home/parkinsonjl/mnt/sharedStorage/middleton/';

local_base_dir = thisDir; %'/home/parkinsonj/convection-in-sea-ice/MushyLayer/matlab/EGUPoster/';
remote_base = '/network/group/aopp/oceans/AW002_PARKINSON_MUSH/middleton/';


TESTING = false;

if nargin < 1
    folder = 'T-10';
end

if nargin < 2
    redoFrames = false;
end

if nargin < 3
    makeVideo = true;
end

output_dir = [local_base_dir, folder, '/'];

% Download data
downloadData(output_dir, folder, remote_base);

% Diagnostics for this simulation:
load([output_dir, 'diagConcat.mat']);

% Make frames from plotdata
files = getFiles(output_dir, 'toPlot');
frames = getFiles(output_dir, 'videoFrame');

numFiles = length(files);


if TESTING
numFiles = 2;
else
  set(0, 'DefaultFigureVisible', 'off');  
end

for i=1:numFiles
    file = files{i};
    %  frame = strrep(file, '.mat', '.png');
    %  frame = strrep(frame, 'toPlot', 'videoFrame');
    
    frameNum = strrep(file, 'toPlot', '');
    frameNum = str2num(strrep(frameNum, '.mat', ''));
    
    frame = sprintf('videoFrame%06d.png', frameNum);
    
    fprintf('Processing %s \n', file);
    
    if ~any(strcmp(frames,frame)) || redoFrames
        % Make frame
        h = makeFrame(output_dir, file, concat_time, concat_Si, concat_Vi);
        
        % Save frame
        h.InvertHardcopy = 'off';
        h.Color = 'white';
        % set(h,'Units','Inches');
        % pos = get(h,'Position');
        % set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        filename = [output_dir, frame];
        %print(h,[filename],'-dpdf','-r0')
        print(h,[filename],'-dpng','-r100')
        
        img = imread(filename);
        fprintf('Image size = %d x %d \n', size(img, 2), size(img, 1));
        
        %if size(img,1) ~= 947 || size(img, 2) ~= 2105
        %    fprintf('WARNING - inconsistent image size');
        %end
        
        %h2 = figure();
        %h2.Position = [0 0 2105 947];
       % imshow(img);
        
       % temp = 0;
    end
    
    close all;
    
end

set(0, 'DefaultFigureVisible', 'on');
if ~TESTING
close all;
end

% Now make the video
if makeVideo
    
    % First, create a list of the frames we want to go into the movie
    frameList = [];
    
     % How long a frame should last for. 
    % If the time between this and the next frame is less than half of
    % this, we will duplicate the current frame.
    % Equally, if the next frame comes too quickly we will skip it
    
    frame_dt = 0.01; 

    frames = getFiles(output_dir, 'videoFrame');
    
    for j=1:length(frames)
        frame = frames{j};   
    
        % Find the frame of this video
        frameNumber = strrep(frame, 'videoFrame', '');
        frameNumber = str2num(strrep(frameNumber, '.png', ''));
        
        frameList(end+1) = frameNumber;
        
         frame_i = find(concat_frames==frameNumber);
         frame_t = concat_time(frame_i);
         
         frame_iNext = frame_i+1;
         frame_tNext = concat_time(frame_iNext);
         
         dt = frame_tNext - frame_t;
         fprintf('dt = %1.5f \n', dt);
%         
%         while (frame_tNext-frame_t) < frame_dt;
%             frame_iNext = frame_iNext+1;
%             frame_tNext = concat_times(frame_iNext);
%         end
%         
%         % Find how many times we need to print this frame
%         numThisFrame = floor((frame_tNext-frame_t)/frame_dt);
%         
%         if length(imgSize) == 0
%             
%             % Size of first frame sets size of video
%             imgSize = size(img);
%         elseif sum(size(img) == imgSize) < 3
%             
%             % Skip frames which don't match video size
%             % Should recalculate these in future
%             fprintf('!!Skipping - wrong size frame!! \n');
%             continue            
%         end
        
    end
    
    
   
    outputVideo = VideoWriter([output_dir,'video.avi']);
    outputVideo.FrameRate = 32;
    open(outputVideo);
    
    firstImage = imread([output_dir, sprintf('videoFrame%06d.png', frameList(1))]);
    imgSize = size(firstImage);
    
    for i = 1:length(frameList)
        fprintf('Making video, frame %d/%d \n', i, length(frameList));
        frame = frameList(i);
        frameFile = sprintf('videoFrame%06d.png', frame);
        img = imread([output_dir, frameFile]);
        
        if sum(size(img) == imgSize) < 3
            fprintf('!!Skipping - wrong size frame!! \n');
            continue
        end;
        
        
        writeVideo(outputVideo,img)
        
    end
    
    close(outputVideo)
    
end


end

function frames = getFiles(d, prefix)
files = dir([d, prefix, '*']);

frames = {};

for i=1:length(files)
    file = files(i);
    frames{end+1} = file.name;
end



%frames = files;
end

function h = makeFrame(d, frame, concat_time, concat_Si, concat_Vi)

h = figure();
scale = 1.053;
h.Position = [300 300 854/scale 480/scale];

dataFile = [d, frame];
depth = 3.0;
doLabel = true;
axPos = [0.06 0.18 0.3 0.7];
SlMax = 0.085; %55g/kg

blues = makeColorMap([8 48 107]/255, [104, 172, 212]/255, [1 1 1]); blues = flipud(blues);

%slMap = makeColorMap([  0.9763    0.9831    0.0538] , [ 0.0488    0.5772    0.8228], [8 48 107]/255);
slMap = makeColorMap([  0.9763    0.9831    0.0538] , [ 0.6488    0.7772    0.2228], [8 48 107]/255);

axPorosity = axes;

load(dataFile);

zmin = 900;
zmax=1024;

zmin = round(1024*(1-depth/8));


Xplot = X(zmin:zmax, :);
Zplot = Z(zmin:zmax, :);
chiPlot = porosity(:, zmin:zmax).';
SlPlot = Sl(:, zmin:zmax).';

pcolor(Xplot,Zplot,chiPlot);

colormap(axPorosity, (blues));

axPorosity.XTick = [0.05 1.95];
axPorosity.XTickLabels = {'0', '2'};
xlabel('$x$ (cm)');

if doLabel
    axPorosity.YTick = [min(min(Z(zmin:zmax, :))) max(max(Z)) ];
    bottomLabel = sprintf('%1.1f', -depth);
    axPorosity.YTickLabels = { bottomLabel, '0.0'};
    ylab = ylabel('$z$ (cm)');
    oldPos = ylab.Position;
    ylab.Position = [oldPos(1)+0.12 oldPos(2)];
else
    axPorosity.YTick = [];
end

daspect([1 1 1]);

axPorosity.Position = axPos;

caxis([0 1]);

%Sl axes
axSl = axes;

smoothTransition = false;

%SlMax = 0.16;
SlDiff  =  min(SlPlot+(1-SlMax), 0);

if smoothTransition
    SlDiff(chiPlot<0.8) = NaN;
else
    SlDiff(chiPlot<1.0) = NaN;
end

SlDiff = -min(min(SlDiff))+SlDiff;



colormap(axSl, flipud(slMap));
%minSl = -1;
avPorosity = mean(chiPlot.^6, 1);

chiFilter = repmat(avPorosity,256,1); %porosity.^2;

if smoothTransition
    SlDiff =SlDiff.*chiFilter;
end

%min(min(SlDiff(:, zmin:zmax)))

pcolor(Xplot, Zplot, SlDiff);
%caxis([min(min(SlDiff)) max(max(SlDiff))]);
caxis([0 SlMax]);
%colorbar('Location', 'eastoutside');

daspect([1 1 1]);

axis(axSl, 'off');
axSl.Visible = 'off';

axSl.Position = axPos;

linkaxes([axPorosity axSl])

dimensionalTime = t*800/60;
if dimensionalTime > 48 && dimensionalTime < 51
    dimensionalTime = 50;
end
timeStr = sprintf('$t=%1.0f$ mins', dimensionalTime);
text(0.08, 8-depth+0.1,timeStr, 'Color', [1 1 1], 'FontSize', 14);


title(['$T_a=',d(end-3:end-1),'^\circ$ C']);


cbSl = colorbar(axSl, 'Location', 'eastoutside');
cbChi = colorbar(axPorosity, 'Location', 'eastoutside');

cbSlPos = cbSl.Position; cbChiPos = cbChi.Position;
cbSl.Position  = [cbSlPos(1)+0.08   cbSlPos(2)-0.1 cbSlPos(3)  0.3];
cbChi.Position = [cbChiPos(1)+0.08 cbChiPos(2)+0.3 cbChiPos(3) 0.3];

cbSl.Label.String = 'S_l (g/kg)';
cbChi.Label.String = '\chi';


cbSl.Ticks = [0 SlMax];
cbSl.TickLabels = {'35', num2str(round(35+SlMax*233))};

cbChi.Ticks = [0.0 0.98];
cbChi.TickLabels = {'0', '1'};



% Now add info
descr = {'Sea ice grown in a Hele-Shaw cell (2cm x 8cm x 1mm).';
    ['Horizontally periodic, open bottom boundary.'];
    ['Initial $S=35$ g/kg, $T=-2^\circ$ C.'];
    ['Colour indicates liquid fraction $\chi$ in the ice,'];
    ['and salinity $S_l$ in the liquid.'];
    ['Sea ice permeability $K=K_0 \chi^3/(1-\chi)^2$; $K_0=10^{-9}$.']
    ['Below: - salt flux $F_S$ (kg m$^{-2}$ s$^{-1}$),'];
    ['$\hspace{8ex}$- total salt release $\Delta S$ (kg m$^{-2}$).']
    
    };

%descr = 'Test';
text(3.2, 8-0.21, descr, 'FontSize', 12);




% Finally, plot diagnostics

t = getTime(concat_time);
mi = concat_Vi*(0.08)*1000; % mass of water/m^2
deltaS =  (35- getSalinity(concat_Si)).*mi; %g/m^2
saltRelease = deltaS/1000;
saltFlux = gradient(saltRelease,t);


doPlot = (~isnan(saltRelease)).*(concat_time ~= 0);
saltFlux = min(saltFlux, 0.018);


diagFluxPos = [0.63 0.4 0.33 0.2];
diagSaltPos = [0.63 0.2 0.33 0.2];

maxi = find(t < dimensionalTime, 1, 'last' );


maxTime = max(t);
% Fs axes
axDiagsFlux = axes;
axDiagsFlux.Position = diagFluxPos;
plot(t(1:maxi), saltFlux(1:maxi));
axDiagsFlux.XTickLabels = {};

axDiagsFlux.YLim = [-0.001 0.02];


ylabel('$F_S$');
axDiagsFlux.YTick = [0 0.01 0.02 0.03 0.04];
axDiagsFlux.YTickLabels = {'0', '','0.02','','0.04'}; 

axDiagsFlux.XLim = [0 maxTime];

% Salt axes
axDiagsSalt = axes;
axDiagsSalt.Position = diagSaltPos;
plot(t(1:maxi), saltRelease(1:maxi));

xlabel('$t$ (mins)');

axDiagsSalt.YTick = [0 0.1 0.2 0.3];
axDiagsSalt.YTickLabels = {'0', '', '0.2', ''};
axDiagsSalt.YLim = [0 0.35];

ylabel('$\Delta S$');


axDiagsSalt.XLim = [0 maxTime];


scale = 1.042;
h.Position = [300 300 854/scale 480/scale];


end


function S = getSalinity(dimlessSalinity)
Se = 233;
deltaS = 233-35;

S = dimlessSalinity*deltaS + Se;
end

function t = getTime(dimlessTime)
timescale = 800; %seconds
t = dimlessTime*timescale/60; % minutes
%t = dimlessTime;
end


function downloadData(local_dir, folder, remote_base)

endings = {'0', 'restart', 'restart2'};
%remote_base = '/home/parkinsonjl/mnt/sharedStorage/middleton/';

T = folder(2:end);

remote_dirs = {};
for i =1:length(endings)  
   remote_dirs{i} = [remote_base, 'CR1.179RaC80Le200KozenyPermeabilitypts256-S3.5-TB',T,'.0-R0.13-',endings{i},'/'];
end

remote_dirs{end+1} = [remote_base, folder, 'Outflow/'];


for i =1:length(remote_dirs)  
   remote_dir = remote_dirs{i};
    
   files = dir([remote_dir, 'toPlot*.mat']);
   for j=1:length(files)
       name = files(j).name;
       
       frame = strrep(name, 'toPlot', '');
       frame = str2num(strrep(frame, '.mat',''));
       
       if length(frame) == 0
           continue
       end
       
       localFileExists = (exist([local_dir, name], 'file') == 2);
       getFile = (mod(frame, 200) == 0);
       if  getFile && ~localFileExists
           
          fprintf('Copying %s \n', name);
          copyfile([remote_dir, name], local_dir);
          
       end
   end
end


end