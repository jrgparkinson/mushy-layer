function Fig9MakeVideo(folder)

set(0, 'defaultlinelinewidth',2);
set(0, 'defaultaxeslinewidth',2);
set(0, 'defaultpatchlinewidth',2);
set(0, 'defaultAxesFontSize',20);
set(groot, 'DefaultTextInterpreter', 'latex')
set(groot, 'DefaultLegendInterpreter', 'latex')
set(groot,'DefaultSurfaceEdgeColor','none') % removes edges from pcolor plots

font = 'Latin Modern Math';
set(groot, 'defaultAxesFontName', font); 
set(groot, 'defaultTextFontName', font); 

if nargin < 1
    
    % Remote
    baseFolder = '/network/group/aopp/oceans/AW002_PARKINSON_MUSH/';
    
    % local
    baseFolder = '/home/parkinsonjl/mnt/sharedStorage/'; 
    folder= [baseFolder, 'TestDiffusiveTimescale/FixedChill-t5.0e-02-Ra1e+06-Da4.8e-04-C2.00-Rel1.0e-04-0/'];
end

framePrefix = 'frame';
plot_prefix = 'plt';

videoFolder = fullfile(folder, 'video');
if exist(videoFolder, 'dir') ~= 7
    fprintf('Making directory %s \n', videoFolder);
    mkdir(videoFolder);
end

testing = false;

options.SlMax = 0.1;
options.axPos = [0.15 0.15 0.65 0.75];
options.topFraction = 0.5;
options.pcolorShadingInterp = true;
options.fluidVelScale = 200;
options.numFrames = 1;
options.pcolorField = 'S'; %porosity, Sl
options.subPlotLabel = ''; % no a/b/c/ label
options.timeTitle = true;
options.cbarPos = [];


close all;


% 1) make frames if they don't exist
plotFiles = dir(fullfile(folder, [plot_prefix, '*']));

filesToDo = 1:length(plotFiles);
%filesToDo = 23:23; %testing
if testing
    filesToDo = 1:20;
end

frameList = [];

% Start timer
tic;

for i=filesToDo

try
    
     plotFile = plotFiles(i).name;
    
    regexStr = ['[^\d]+(\d+)\.2d\.hdf5'];
    [tokens,matches] = regexp(plotFile,regexStr,'tokens','match');
    
    thisFrame = str2num(tokens{1}{1});
    frameList(i) = thisFrame;
    
    close all;
    
    frame  = -1;
    prefix = '';
    
    saveFile = fullfile(videoFolder, [framePrefix, num2str(frameList(i)) ]);
    
    if exist([saveFile, '.png'], 'file') == 2
        fprintf('frame %d exists \n', i);
        %todo - uncomment line below
        continue
    else
        fprintf('Making frame %d \n', i);
    end
    
    h = figure();
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(h,'Units','Inches');
    figPos = [2.0 2.0 4.5 2.0] ;

    h.Position = figPos;
    textFontSize = 9;
    legendFontSize = 8;
    domainFontSize = 8;

    set(0, 'defaultlinelinewidth',1);
    set(0, 'defaultaxeslinewidth',1);
    set(0, 'defaultpatchlinewidth',1);
    set(0, 'defaultAxesFontSize', textFontSize);
    set(0, 'defaultAxesFontName', 'times');
    set(0, 'defaultTextFontName', 'times');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    set(h, 'Position', figPos);
    fprintf('Setting figure position: [%d %d %d %d]\n', figPos(1), figPos(2), figPos(3), figPos(4));
      
    %success = makeFrameFunction(folder, plotFile, options);
    
    success = Fig9PlotFrame(options, folder, plot_prefix, thisFrame, 1); 
    
    % pause(0.5);
    if success
        % Set position again to be sure
        set(h, 'Position', figPos);
        fprintf('Setting figure position: [%d %d %d %d]\n', figPos(1), figPos(2), figPos(3), figPos(4));
           
        % Make sure we specify the paper position
        paperPos =  [0 0 figPos(3) figPos(4)];
        set(h, 'PaperUnits', 'inches', 'PaperPosition', paperPos);
        h.InvertHardcopy = 'off';
        h.Color = 'white';
        print(h,[saveFile, '.png'],'-dpng','-r200')
        
        fprintf('Saved to %s \n', saveFile);
        
    end
    
    elapsedTime = toc;
    timePerRun = elapsedTime/i;
    timeLeft = (length(filesToDo)-i)/timePerRun;
    fprintf('Elapsed time: %f, predicted remaining time: %f \n', elapsedTime, timeLeft);
    
   catch e
   fprintf('Exception \n ');
   fprintf('%s', e.message);
   end
   
end


% Now compile frames into one video
%frameList=filesToDo;

frameSpacing = 2;
fps = 32;
video_name = sprintf('video-%dfps-%dframeSkip.avi', fps, frameSpacing);
outputVideo = VideoWriter(fullfile(videoFolder, video_name));
outputVideo.FrameRate = fps; 
open(outputVideo);



firstImage = imread(fullfile(videoFolder, sprintf('frame%d.png', frameList(1) ) ) );
imgSize = size(firstImage);


% Skip first 10 frames
framesToMakeVideo = 1:frameSpacing:length(frameList);
%if testing
%framesToMakeVideo = 400:410;
%end

for i = framesToMakeVideo
    fprintf('Making video, frame %d/%d \n', i, length(frameList));
    
    frame = frameList(i);
    frameFile = sprintf('frame%d.png', frame);
    fullFilename = fullfile(videoFolder, frameFile);
    if exist(fullFilename, 'file') == 2
        
        img = imread(fullFilename);
        
        % Trim size:
        %img = img(1:880, 1:3100, :);
        imgSize = size(img);
        
        fprintf('Image size: %dx%d \n', imgSize(1), imgSize(2));
        
        % if sum(size(img) == imgSize) < 3
        %     fprintf('!!Skipping - wrong size frame!! \n');
        %     continue
        % end;
        
        
        writeVideo(outputVideo, img)
        
    else
        fprintf('Frame does not exist');
    end
    
end

fprintf('Writing video \n');
close(outputVideo)

 elapsedTime = toc;
 
fprintf('Wrote video, script took %1.1f seconds.', elapsedTime );




end




