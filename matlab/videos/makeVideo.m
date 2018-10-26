function makeVideo(run_name)

if nargin < 1
    run_name = 'CR10.000RaC200Le100KozenyPermeabilityDa1.0e-03R1.0e-03pts512-width8.0-aspect8.0-0';
    run_name = 'CR10.000RaC200Le100KozenyPermeabilityDa1.0e-03R1.0e-03pts1024-width8.0-aspect4.0-0';
end

close all;

base_dir = getDataDir('channelSpacingV2');

% Testing:
%base_dir = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/matlab/videos/';
%run_name = 'testing';

folder = fullfile(base_dir, run_name);


% 1) make frames if they don't exist
plotFiles = dir(fullfile(folder, 'plt*'));

filesToDo = 1:length(plotFiles);
%filesToDo = 1:3; %testing

% Start timer
tic;

for i=filesToDo
    
    frame=-1;
    prefix='';
    
    saveFile = fullfile(folder, ['frame',num2str(i)]);
    
    if exist([saveFile, '.png'], 'file') == 2
        fprintf('frame %d exists \n', i);
        continue
    else
        fprintf('Making frame %d \n', i);
    end
    
    plotFile = plotFiles(i).name;
    
    folders_regex = '(plt.*)(\d+).2d.hdf5';
    [~,matches,~]  = regexp(plotFile, folders_regex, 'match', 'tokens', 'tokenExtents');
    
    if isempty(matches)
        continue;
    else
        match = matches{1};
        prefix = match{1};
        frame = str2double(match{2});
    end
    
    
    % Load data
    output = MushyLayerOutput(2, frame, folder, prefix, true);
    
    
    porosity = output.dataForComp(output.components.Porosity);
    streamfunction = output.getStreamfunction(100, 1);
    
    xvel = output.dataForComp(output.components.xAdvectionvelocity);
    yvel = output.dataForComp(output.components.yAdvectionvelocity);
    
    Sl = output.dataForComp(output.components.Liquidconcentration);
    
    probDomain = output.problemDomain;
    dx = probDomain.dxCoarse;
    
    numx = (probDomain.domainExtent.hi_i - probDomain.domainExtent.lo_i) + 1;
    numy = (probDomain.domainExtent.hi_j - probDomain.domainExtent.lo_j) + 1;
    
    xlo = double(probDomain.domainExtent.lo_i)*dx;
    xhi = double(probDomain.domainExtent.hi_i)*dx;
    width = xhi-xlo;
    
    ylo = double(probDomain.domainExtent.lo_j)*dx;
    yhi = double(probDomain.domainExtent.hi_j)*dx;
    
    % if limits are almost an integer, round them up
    xlo = correctLimit(xlo, dx);
    xhi = correctLimit(xhi, dx);
    ylo = correctLimit(ylo, dx);
    yhi = correctLimit(yhi, dx);
    
    %x = linspace(double(probDomain.domainExtent.lo_i), double(probDomain.domainExtent.hi_i), numx)*dx;
    %y = linspace(double(probDomain.domainExtent.lo_j), double(probDomain.domainExtent.hi_j), numy)*dx;
    % actually, rescale
    x = linspace(double((xlo)), double((xhi)), numx);
    y = linspace(double((ylo)), double((yhi)), numy);
    
    
    
    %zoom_lo_z = zoom_lo_j*dx;
    %zoom_hi_z = zoom_hi_j*dx;
    streamfunction = streamfunction.';
    
    chi = porosity;
    U = xvel.';
    V = yvel.';
    %psi = streamfunction;
    xzoom = x;
    yzoom = y;
    
    [X, Y] = meshgrid(xzoom, yzoom);
    
    % Plot data
    aspectRatio = (max(xzoom)-min(xzoom))/(max(yzoom)-min(yzoom));
    
    figureWidth = 1500;
    figureHeight = round(figureWidth/6.0);
    figureHeight = round(figureWidth/(aspectRatio*0.9));
    
    h = figure();
    figPos = [200 200 figureWidth figureHeight];
    
    set(h, 'Position', figPos);
    
    
    
    %axisExtent = [0.02 0.02 0.96 0.96];
    axisExtent = [0.03 0.05 0.8 0.95];
    
    axSl = axes;
    % This should be the 'blues' colormap from visit
    %colormap(axPorosity, makeColorMap([8 48 107]/255, [104, 172, 212]/255, [1 1 1]));
    %colormap(axSl, flipud(parula));
    colormap(axSl, parula);
    caxis(axSl, [-1 0]);
    pcolor(xzoom, yzoom, Sl.');
    set(axSl,'dataAspectRatio',[1 1 1])
    set(axSl,'Ydir','normal');
    
    caxis(axSl, [-1 0]);
    
    % xlabel('$x$'); ylabel('$z$');
    axSl.XTick = [];
    axSl.YTick = [];
    
    c = colorbar();
    c.Label.String='Salinity';
    c.Ticks = [-0.99 -0.01];
    c.TickLabels = {'0', '1'};
    
    cPos = c.Position;
    
    c.Position = [0.85 0.2 cPos(3) cPos(4)];
    
    set(axSl, 'position', axisExtent);
    set(axSl,'TickLength',[0 0])
    
    axSl.LineWidth=3.0;
    
    box on;
    
    
    % Now plot mush-liquid interface
    
    axPorosity = axes;
    
    %colormap(axSl, makeColorMap([0 0 0] , [230, 175, 0]/255, [200 0 0]/255));
    %salinityColormap = makeColorMap([0 0 0] ,  [200 0 0]/255);
    %salinityColormap = makeColorMap([1 1 1] ,  [0.5 0.5 0.5], [200 0 0]/255);
    %salinityColormap = makeColorMap([200 0 0]/255, [200 0 0]/255);
    
    %salinityColormap = makeColorMap([1 1 1] ,  [1 1 1], [1 1 1]);
    % salinityColormap = makeColorMap([0 0 0] ,  [0 0 0]);
    
    %colormap(axPorosity, salinityColormap);
    caxis(axPorosity, [0 1]);
    
    [CSl, hSl] = contour(X, Y, chi.', [-0.1 0.95]);
    
    
    hSl.LineWidth = 4.0;
    hSl.LineStyle = '--';
    hSl.LineColor='r';
    
    set(axPorosity,'dataAspectRatio',[1 1 1])
    set(axPorosity, 'position', axisExtent);
    set(axPorosity,'TickLength',[0 0])
    
    axis(axPorosity, 'off');
    axPorosity.Visible = 'off';
    
    linkaxes([ axSl axPorosity ]);
    
    
    box on;
    
    % Todo: add ice/liquid labels
    fs = 26;
    text(0.5, max(yzoom)-0.13, 'Porous Ice', 'FontSize', fs)
    text(0.5, 0.15, 'Liquid', 'FontSize', fs, 'Color', 'white')
    
    text(mean(xzoom)*0.8, -0.2, 'Warm inflow/outflow', 'FontSize', fs, 'Color', 'red')
    text(mean(xzoom)*0.8, max(yzoom)+0.15, 'Cold boundary', 'FontSize', fs, 'Color', 'blue')
    
    text(max(xzoom)*0.97, min(yzoom)-0.3, '- - Ice-liquid interface', 'FontSize', fs*0.8, 'Color', 'red')
    
    
    
    
    
    % Trying to ensure all figures are same size
    set(axPorosity, 'position', axisExtent);
    set(axSl, 'position', axisExtent);
    
    
    
    pause(0.5);
    set(h, 'Position', figPos);
    fprintf('Setting figure position: [%d %d %d %d]\n', figPos(1), figPos(2), figPos(3), figPos(4));
    pause(0.5);
    
    h.InvertHardcopy = 'off';
    h.Color = 'white';
    print(h,[saveFile, '.png'],'-dpng','-r200')
    
    fprintf('Saved to %s \n', saveFile);
    
    elapsedTime = toc;
    timePerRun = elapsedTime/i;
    timeLeft = (length(filesToDo)-i)/timePerRun;
    fprintf('Elapsed time: %f, predicted remaining time: %f \n', elapsedTime, timeLeft);
    
end


% Now compile frames into one video
frameList=filesToDo;

outputVideo = VideoWriter(fullfile(folder, 'video.avi'));
outputVideo.FrameRate = 16;
open(outputVideo);

firstImage = imread(fullfile(folder, sprintf('frame%d.png', frameList(1) ) ) );
imgSize = size(firstImage);


% Skip first 10 frames
for i = 10:length(frameList)
    fprintf('Making video, frame %d/%d \n', i, length(frameList));
    
    frame = frameList(i);
    frameFile = sprintf('frame%d.png', frame);
    img = imread(fullfile(folder, frameFile));
    
    % Trim size:
    %img = img(1:880, 1:3100, :);
    imgSize = size(img);
    
    fprintf('Image size: %dx%d \n', imgSize(1), imgSize(2));
    
    % if sum(size(img) == imgSize) < 3
    %     fprintf('!!Skipping - wrong size frame!! \n');
    %     continue
    % end;
    
    
    writeVideo(outputVideo, img)
    
end

close(outputVideo)



end



function output = correctLimit(lim, dx)

output = lim;

if abs(round(lim)-lim) < 2*dx
    output = round(output);
end


end