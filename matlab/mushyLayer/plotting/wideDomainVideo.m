%clear all;
close all;

set(0, 'defaultlinelinewidth',2);
set(0, 'defaultaxeslinewidth',2);
set(0, 'defaultpatchlinewidth',2);
set(0, 'defaultAxesFontSize',22);


baseDir = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/';
dataFolder = '/media/parkinsonjl/FREECOM HDD/';
%output_dir = [dataFolder, 'mushyLayer-periodic/CR1.25RaC250.0Le200ChiCubedPermeabilitypts512-0/'];
output_dir = '/media/parkinsonjl/FREECOM HDD/mushyLayer-periodic/CR1.25RaC500.0Le200ChiCubedPermeabilitypts512-0-wavey/'; % wavey
makeFrames = true;
plotDiagnostics = true;
hideFigures = true;
makeVideo = true;
videoName = 'waveyDirectionalSolidification';

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

%frames = frames(30);

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

inputs = readInputs([output_dir, 'inputs']);
L = -1;

if isfield(inputs, 'height')
    L = inputs.height;
end

RaC = -1;
if isfield(inputs, 'rayleighComp')
   RaC = inputs.rayleighComp;
end

% numCellsStr = strsplit(inputs.num_cells, ' ');
% numCells = [];
% for i = 1:length(numCellsStr)
%     numCells(i) = str2num(numCellsStr{i});
% end
% aspectRatio = numCells(2)/numCells(1);
% 

if makeFrames
    
    for frame_i = 1:length(frames)
        frame = frames(frame_i);
        dim = 2; subcycled = true;
        
        output = MushyLayerOutput(dim, frame, output_dir, plot_prefix, subcycled);
        
        time = output.t;
        
        porosity = output.dataForComp(output.components.Porosity);
        permeability = output.dataForComp(output.components.Permeability);
        %streamfunction = output.getStreamfunction(100, 1);
        
        xvel = output.dataForComp(output.components.xAdvectionvelocity);
        yvel = output.dataForComp(output.components.yAdvectionvelocity);
        
        Sl = output.dataForComp(output.components.Liquidconcentration);
        
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
        
        contourValsSl = linspace(-0.98, 0, 7);
        
        h = figure();
        if hideFigures
            set(gcf,'Visible', 'off');
        end
        h.Units = 'pixels';
        h.Position = [100 100 1600 800];
        
        
        
        topAxisExtent = [0.25 0.55 0.65 0.38];
        bottomAxisExtent = topAxisExtent;
        bottomAxisExtent(2) = bottomAxisExtent(2) - 0.45;
        
        m = 2; n = 1;
        
        axPorosity = subplot(m, n, 1);
        
        
        % This should be the 'blues' colormap from visit
        colormap(axPorosity, makeColorMap([8 48 107]/255, [104, 172, 212]/255, [1 1 1]));
        caxis(axPorosity, [0 1]);
        imagesc(xzoom, yzoom, chi.');
        set(axPorosity,'dataAspectRatio',[1 1 1])
        set(axPorosity,'Ydir','normal');
        c = colorbar('eastoutside');
        c.Label.String = '\chi';
        %c.Label.Rotation = 0;
        %c.Label.Position = [1.5 0.8];
        c.Ticks = [0.01, 0.5, 1.0];
        c.TickLabels = {'0', '0.5', '1'};
        
        %axPorosity.XTick = [];
        %axPorosity.YTick = [];
        set(axPorosity, 'position', topAxisExtent);
        %set(axPorosity,'TickLength',[0 0])
        
        title(['t=', num2str(time)]);
        xlabel('x');
        ylabel('y');
        set(get(gca,'YLabel'),'Rotation',0)
        
        box on;
        
        axSl = axes;
        
        %colormap(axSl, makeColorMap([0 0 0] , [230, 175, 0]/255, [200 0 0]/255));
        colormap(axSl, makeColorMap([0 0 0] ,  [200 0 0]/255));
        caxis(axSl, [-1 0]);
        [CSl, hSl] = contour(X, Y, Sl.', contourValsSl);
        
        c2 = colorbar;
        c2.Label.String = 'S_l contours';
        c2.Location = 'westoutside';
        originalPosition = c2.Position;
        c2.Position = [0.08 topAxisExtent(2) 0.02 topAxisExtent(4)-0.05];
        c2.Ticks = [-0.97 -0.5 0.0];
        c2.TickLabels = {'-1.0', '-0.5', '0.0'};
        
        
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
        
        
        axBottom = subplot(m,n,2);
        colormap(axBottom, makeColorMap([8 48 107]/255, [104, 172, 212]/255, [1 1 1]));
        caxis(axBottom, [0 1]);
        imagesc(xzoom, yzoom, permeability.');
         cPerm = colorbar('eastoutside');
        cPerm.Label.String = ('\Pi = \chi^3');
        cPerm.Ticks = [0.01 0.5 0.99];
        cPerm.TickLabels = {'0', '0.5',  '1'};
        
        xlabel('x'); ylabel('y');
        
        set(axBottom,'dataAspectRatio',[1 1 1])
        set(axBottom,'Ydir','normal');
        set(axBottom, 'position', bottomAxisExtent);
        
       
        
        axFlow = axes;
        
        %streamline(X, Y, U, V);
        %hStream = streamslice(X, Y, U, V, 2, 'method', 'cubic');
        skip = 10;
        hQuiver = quiver(X(1:skip:end, 1:skip:end), Y(1:skip:end, 1:skip:end),...
            U(1:skip:end, 1:skip:end), V(1:skip:end, 1:skip:end));
        hQuiver.AutoScaleFactor = 1.5;
        hQuiver.Color = 'red';
        hQuiver.LineWidth = 2.0;
        %colormap(axFlow, [0 0 0; 0 0 0]);
        axis(axFlow, 'off');
        axFlow.Visible = 'off';
        set(axFlow,'dataAspectRatio',[1 1 1])
        set(axFlow, 'position', bottomAxisExtent);
        
        %hStream.Color = 'k';
        
        linkaxes([axBottom axFlow]);
        
        whiteNoise = true;
        initialPert = 'small white noise';
        if isfield(inputs, 'restart_wavenumber') && inputs.restart_wavenumber ~= 0
            whiteNoise = false;
            initialPert = ['sinusoidal perturbation, k = ', num2str(inputs.restart_wavenumber)];
        end
        
        descr = {'Directional solidification with dimensionless units:';
         ['Ra_C: ' , num2str(RaC)];
         ['Stefan:', num2str(inputs.stefan)];
         ['Le:', num2str(inputs.lewis)];
         ['Concentration Ratio:', num2str(inputs.compositionRatio)];
         ['Frame translation speed:', num2str(inputs.nonDimVel)];
         ['Horizontally periodic'];
         ['Bottom inflow: H = ', num2str(inputs.bottomEnthalpy),', C = ', num2str(inputs.bottomBulkConc)];
         ['Bottom outflow: extrapolation'];
         ['Top: H = ', num2str(inputs.topEnthalpy), ', C = ', num2str(inputs.topBulkConc)];
         ['Initial condition: '];
         ['   Diffusive steady state'];
         ['   + ', initialPert];
        };
    
        text(-0.28, 0.1, descr);
        
        
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