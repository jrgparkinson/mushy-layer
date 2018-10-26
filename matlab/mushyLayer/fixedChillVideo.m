%clear all;
close all;

set(0, 'defaultlinelinewidth',2);
set(0, 'defaultaxeslinewidth',2);
set(0, 'defaultpatchlinewidth',2);
set(0, 'defaultAxesFontSize',22);


baseDir = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/';
dataFolder = '/media/parkinsonjl/FREECOM HDD/';
%output_dir = [dataFolder, 'mushyLayer-periodic/CR1.25RaC250.0Le200ChiCubedPermeabilitypts512-0/'];
%output_dir = '/media/parkinsonjl/FREECOM HDD/fixedChill-insulating/KozenyPermeabilitypts512-1/';
output_dir = '/media/parkinsonjl/FREECOM HDD/fixedChill-insulating/KozenyPermeabilitypts256-0-stopped/';
makeFrames = false;
plotDiagnostics = true;
hideFigures = true;
videoName = 'fixedChill';
makeVideo = true;

%Let's try and work everything else out from just the output dir
files = dir([output_dir, '*.2d.hdf5']);
plotFiles = {};
frames = [];
pattern = '(.*\-)(\d+)';
for i=1:length(files)
    fileName = files(i).name;
    TF = strfind(fileName, 'chk');
    if length(TF) == 0
        plotFiles{end+1} = fileName;
        [~,tok,~]  = regexp(fileName, pattern, 'match', ...
            'tokens', 'tokenExtents');
        plot_prefix = tok{1}{1};
        frame = tok{1}{2};
        frames(end+1) = str2num(frame);
    end
end

initialFrame = frames(1);

%tempFrames = frames;
%frames = [tempFrames(60)];

% plot_prefix = 'mushyLayer-periodic-CR1.25RaC250.0Le200ChiCubedPermeabilitypts512-';
%
% RaC = 250;
% % Get diagnostic info to plot on graph next to video
% pout = Pout([output_dir, 'pout.0']);
% % Get frames
% frames = 8000:400:54800; % should be 400
%
% % Length scale (metres)
% L = 2;

% already loaded these
pout = Pout([output_dir, 'pout.0']);
inputs = readInputs([output_dir, 'inputs']);
L = -1;

if isfield(inputs, 'height')
    L = str2num(inputs.height);
end

timescale = -1;
if L > -1
    timescale = pout.timescale;
end

RaC = -1;
if isfield(inputs, 'rayleighComp')
    RaC = str2num(inputs.rayleighComp);
end

numCellsStr = strsplit(inputs.num_cells, ' ');
numCells = [];
for i = 1:length(numCellsStr)
    numCells(i) = str2num(numCellsStr{i});
end
aspectRatio = numCells(2)/numCells(1);


% For converting salinities
S_e = str2double(inputs.eutecticComposition);
S_inf = str2double(inputs.initialComposition);

% S = S_e + Theta*(S_e - S_inf)

if makeFrames
    
    for frame_i = 1:length(frames)
        frame = frames(frame_i);
        dim = 2; subcycled = true;
        
        % already loaded this
        output = MushyLayerOutput(dim, frame, output_dir, plot_prefix, subcycled);
        
        time = output.t;
        
        if timescale > -1
            % Convert to hours
            time = time*timescale/3600;
        end
        
        porosity = output.dataForComp(output.components.Porosity);
        %streamfunction = output.getStreamfunction(100, 1);
        
        xvel = output.dataForComp(output.components.xAdvectionvelocity);
        yvel = output.dataForComp(output.components.yAdvectionvelocity);
        
        Sl = output.dataForComp(output.components.Liquidconcentration);
        
        % Make dimensional
        Sl = S_e + Sl*(S_e - S_inf);
        
        probDomain = output.problemDomain;
        dx = probDomain.dxCoarse;
        
        numx = (probDomain.domainExtent.hi_i - probDomain.domainExtent.lo_i) + 1;
        numy = (probDomain.domainExtent.hi_j - probDomain.domainExtent.lo_j) + 1;
        
        x = linspace(double(probDomain.domainExtent.lo_i), double(probDomain.domainExtent.hi_i), numx)*dx;
        y = linspace(double(probDomain.domainExtent.lo_j), double(probDomain.domainExtent.hi_j), numy)*dx;
        
        
        xlo = double(probDomain.domainExtent.lo_i)*dx;
        xhi = double(probDomain.domainExtent.hi_i)*dx;
        width = xhi-xlo;
        
        
        %zoom_lo_z = zoom_lo_j*dx;
        %zoom_hi_z = zoom_hi_j*dx;
        %streamfunction = streamfunction.';
        
        chi = porosity;
        U = xvel.';
        V = yvel.';
        %psi = streamfunction;
        xzoom = x;
        yzoom = y;
        
         % Convert to cm (L is in metres)
        if L > -1
            xzoom = xzoom*L*100;
            yzoom = yzoom*L*100;
        end
        
        % Stretch slightly as otherwise the plot extent is e.g. 0 to 0.9875
        xzoom = xzoom*1.01;
        yzoom = yzoom*1.01;
        
        [X, Y] = meshgrid(xzoom, yzoom);
        
        
        % Mush liquid boundary
        %[mush_x, mush_y] = find(chi<1.0);
        %idx_liquid = 1-idx_mush;
        %refinement = 16;
        %chiRef = resizem(chi, refinement, 'bicubic');
        %Xref = resizem(X, refinement, 'bilinear');
        %Yref = resizem(Y, refinement, 'bilinear');
        %chiRef = smoothn(chiRef, 'robust');
        
        % Velocity
        %offset = 5.0; numContours = 6; maxVal = 90;
        %contourValsNeg = linspace(-maxVal, -offset, numContours);
        %contourValsPos = linspace(offset, maxVal, numContours);
        
        contourValsSl = linspace(min(min(Sl))*1.06, S_e, 9);
        
        h = figure();
        if hideFigures
            set(gcf,'Visible', 'off');
        end
        h.Units = 'pixels';
        figPosition = [100 100 2600 1600];
        h.Position = figPosition;
        
        %Region to plot: 
      
        
        starti = 1;
        endi = length(xzoom);
        startj = 1; %round(length(yzoom)*0.2);
        endj = length(yzoom);
        
        
        
        m = 2; n = 2;
        axPorosity = subplot(m, n, [1,3]);
        
        
        topAxisExtent = [0.05 0.25 0.5 0.7];
        bottomAxisExtent = [0.6 0.3 0.3 0.3];
        
        leftColorBarPos = [0.05 0.1 0.18 0.03];
        rightColorBarPos = leftColorBarPos;
        rightColorBarPos(1) = rightColorBarPos(1) + 0.25;
        
        %axPorosity = axes;
        % This should be the 'blues' colormap from visit
        colormap(axPorosity, makeColorMap([8 48 107]/255, [104, 172, 212]/255, [1 1 1]));
        caxis(axPorosity, [0 1]);
        imagesc(xzoom(starti:endi), yzoom(startj:endj), chi(starti:endi, startj:endj).');
        set(axPorosity,'dataAspectRatio',[1 1 1])
        set(axPorosity,'Ydir','normal');
        c = colorbar('southoutside');
        c.Ticks = [0.01, 0.5, 1.0];
        c.TickLabels = {'0', '0.5', '1'};
        %c.Location = 'southoutside';
        c.Position = leftColorBarPos;
        
        c.Label.String = '\chi';
        %c.Label.Rotation = 0;
        %c.Label.Position = [-1.0 0.8];
        
        
        
        set(axPorosity, 'position', topAxisExtent);
        
        titleStr = '';
        
        if RaC > -1
            titleStr = [titleStr, 'RaC = ', num2str(RaC), ', '];
        end
        
        titleStr = [titleStr, 't=', num2str(time)];
        
        if timescale > -1
            titleStr = [titleStr, 'hrs'];
        end
        
        title(titleStr);
        
        if L > -1
            xlabel('x (cm)');
            ylab = ylabel('y (cm)');
            %ylab.Position = [-0.8  1.2];
        else
            xlabel('x');
            ylab = ylabel('y');
        end
        
         % Add panel with parameters
        descr = {'Solidification of sea water in a Hele-Shaw cell.';
    ['Cell wall separation: ', num2str(str2double(inputs.d)*1e3), 'mm.'];
    ['Top temperature: -22 celcius.'];
    ['Initial temperature: ', inputs.bottomTemperature, ' celcius, salinity: ', inputs.initialComposition, 'g/kg.'];
    'Insulating side and bottom walls.';
    'Permeability: Kozeny-Carman.';
    'Similar setup to Middleton et. al. (2016)'};

    text(25, 18, descr, 'FontSize', 16);
        
        %set(get(gca,'YLabel'),'Rotation',0);
        %set(get(gca,'YLabel'),'Position',[0.0 0.8]);
         
        box on;
        
        %axis([0 20 10 20]);
        
        axSl = axes;
        
        %colormap(axSl, makeColorMap([0 0 0] , [230, 175, 0]/255, [200 0 0]/255));
        colormap(axSl, makeColorMap([0 0 0] ,  [200 0 0]/255));
        caxis(axSl, [-1 0]);
        [CSl, hSl] = contour(X(startj:endj, starti:endi), ...
            Y(startj:endj, starti:endi), ...
            Sl(starti:endi, startj:endj).', contourValsSl);
        
        %axis([0 20 10 20]);
        
        c2 = colorbar('southoutside');
        %c2.Location = 'southoutside';
        c2.Position = rightColorBarPos;
        c2.Ticks = [55  228];
        c2.TickLabels = {'35', '230'};
        
        c2.Label.String = 'S_l (g/kg)';
        %c2.Label.Position = [-2.0 0.8];
        
        
        hSl.LineWidth = 2.0;
        
        set(axSl,'dataAspectRatio',[1 1 1])
        set(axSl, 'position', topAxisExtent);
        
        %set(axSl,'TickLength',[0 0])
        
        axis(axSl, 'off');
        axSl.Visible = 'off';
        
       
        
        %axis([xzoom(1)+dx xzoom(end)-dx yzoom(1)+dx yzoom(end)-dx]);
        %ax = gca;
        %set(gca,'linewidth',2); %make box border bold to hide chimney outline
        %set(ax,'clim',[0 1]);
        %set(gca, 'Layer', 'top'); % put box border above imagesc
        
        linkaxes([axPorosity axSl]);
        
        diagAx = subplot(m, n, 4);
        % Plot of diagnostics (e.g. salt flux) to this point in time
        step = frame + 1 - initialFrame;
        
        step = min(step, length(pout.times));
        
        times = pout.times(1:step);
        if timescale > -1
            times = times*timescale/3600;
        end
        
        % Make dimensional. Only take y=0,4,8
        dimensionalSalinity =  S_e + pout.averageSalinity(:, 1:3)*(S_e - S_inf);
        
        
        plot(times, dimensionalSalinity(1:step, :));
        if timescale > -1
            xlabel('time (hrs)');
        else
            xlabel('time');
        end
        ylabel('Average S_l (g/kg)');
        
        
        maxTime = pout.times(end)*timescale/3600;
        
        axis([times(1) maxTime ...
            min(min(dimensionalSalinity)) max(max(dimensionalSalinity))*0.98]);
        box on;
        
        leg = {};
        
        % Height in cm
        if isfield(inputs, 'domain_length')
            height = str2double(inputs.domain_length)*L*100;
        else
            height = yzoom(end)*L*100;
        end
        
        for i = 1:3
            leg{i} = ['y = ', num2str((i-1)*height/5), ' cm'];
        end
        
        legend(leg, 'Location', 'northwest');
        
        set(diagAx, 'position', bottomAxisExtent);
        
        % Do this last
       
        dpi = 150;
        
        h.PaperUnits = 'inches';
        h.PaperPosition = [0 0 figPosition(3) figPosition(4)]/dpi;
        
        %set(h,'Units','Inches');
        %pos = get(h,'Position');
        %set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        outputDir = [pwd, '/videos/',videoName,'RaC', num2str(RaC), '/'];
        
        if ~exist(outputDir, 'dir')
            mkdir(outputDir);
        end
        
        print(h,[outputDir, 'frame',num2str(frame),'.png'],'-dpng','-r150')
        
    end
end %end if make frames


% Now make the video
if makeVideo
outputVideo = VideoWriter([pwd,'/videos/',videoName,'RaC', num2str(RaC), '.avi']);
outputVideo.FrameRate = 8;
open(outputVideo);

for ii = 1:length(frames)
    frame = frames(ii);
    img = imread([pwd, '/videos/',videoName,'RaC', num2str(RaC), '/frame',num2str(frame),'.png']);
    writeVideo(outputVideo,img)
end

close(outputVideo)

end
