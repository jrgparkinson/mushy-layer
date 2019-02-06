% This scripts makes use of the Streamfunction Color function:
% https://uk.mathworks.com/matlabcentral/fileexchange/69268-streamfunction-color

function Fig9FixedChill(output_dir, frames, figure_output_dir)

if nargin < 1
    
    % Remote
    %baseDir = '/network/group/aopp/oceans/AW002_PARKINSON_MUSH/';
    
    % Local
    baseDir = '/home/parkinsonjl/mnt/sharedStorage/'; 
    
    output_dir = [baseDir, 'TestDiffusiveTimescale/FixedChill-t5.0e-02-Ra1e+06-Da4.8e-04-C2.00-Rel1.0e-04-0/'];
    
end

if nargin < 2
    frames = [-1]; % just process the last frame by default
    frames = [3000, 4800, 17000]; % Testing
    frames = [2288, 3874, 4576]; % Testing
    %frames = [2286, 3872, 4574]; % Testing
end

if nargin < 3
    figure_output_dir = '/home/parkinsonjl/mnt/sharedStorage/TestDiffusiveTimescale/';
end


processSpecificFolder(output_dir, frames, figure_output_dir);
 
end


function processSpecificFolder(output_dir, frames, figure_output_dir)



savePNG = true;

close all; 

domWidth = 1.0;
domHeight = 1.0;
plotHeight = 0.5; % plot top 0.5 of domain

plotAspectRatio = domWidth/plotHeight;

plotScreenWidth = 500;
plotScreenHeight = plotScreenWidth/plotAspectRatio;

% Start plotting stuff
h = figure();
%h.Position = [200 200 plotScreenWidth+70 plotScreenHeight*length(frames)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(h,'Units','Inches');
width = 3.5;
height = (plotScreenHeight*length(frames)/(plotScreenWidth+70)) * width;
h.Position = [2.0 2.0 width height];
textFontSize = 8;
legendFontSize = 8;
domainFontSize = 8;

set(0, 'defaultlinelinewidth',1);
set(0, 'defaultaxeslinewidth',1);
set(0, 'defaultpatchlinewidth',1);
set(0, 'defaultAxesFontSize', textFontSize);

font = 'Latin Modern Math';
set(0, 'defaultAxesFontName', font); 
set(0, 'defaultTextFontName', font); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Common for all plots
options.pcolorShadingInterp = true;
options.topFraction = plotHeight;
options.fluidVelScale = 200;
options.pcolorField = 'S'; %porosity, Sl
options.timeTitle = false;

% Setup axes based on number of plotfiles
axPositions = {};

options.numFrames = length(frames);
pltWidth = 0.75;
pltHeight = 0.8/options.numFrames;

for i=1:length(frames)
   axPositions{end+1} = [0.07, 0.06 + (options.numFrames-i)*(pltHeight+0.06), pltWidth, pltHeight];
end

% Get all frames in folder
[actual_plot_prefix, allFrames] = getFrames(output_dir, '/');

% Load plotfile
for frame_i=1:length(frames)
    options.axPos = axPositions{frame_i}; %[0.08 0.08 0.8 0.75];
    options.maxSearchFrame = 100000;
    thisFrame = frames(frame_i);
    
   
    
    textLabels = {'(a)', '(b)', '(c)'};
    options.subPlotLabel = textLabels{frame_i};
    
    Fig9PlotFrame(options, output_dir, actual_plot_prefix, thisFrame, frame_i);
    
end

% h = gcf;
h.InvertHardcopy = 'off';
h.Color = 'white';
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
filename = fullfile(figure_output_dir, 'Fig9FixedChillSimulation');

if options.pcolorShadingInterp
    filename = [filename, '-shadingInterp'];
else
    filename = [filename, '-noShadingInterp'];
end

%%if savePDF
% print(h,[filename, '.pdf'],'-dpdf','-r0')
% end

if savePNG
    print(h,[filename, '.png'],'-dpng','-r800')
end

% this is abot 40MB
%print(h,[filename, '.eps'],'-depsc','-r50')



    
end




