% This scripts makes use of the Streamfunction Color function:
% https://uk.mathworks.com/matlabcentral/fileexchange/69268-streamfunction-color

function Fig7FixedChill(output_dir, frames, figure_output_dir)

if nargin < 1
    output_dir = '/home/parkinsonjl/mnt/sharedStorage/Test/FixedChill-t100.0-Ra1e+07-Da5.0e-05-C5.00/Uniform-FixedChill-64--0/';
    output_dir = '/home/parkinsonjl/mushy-layer/execSubcycle/amr/';
    output_dir = '/home/parkinsonjl/mnt/sharedStorage/fixedChillAMR/TBottom1.3AMR4-2-0/';
end

if nargin < 2
    frames = [-1]; % just process the last frame by default
    frames = [3000, 4800, 17000]; % Testing
end

if nargin < 3
    figure_output_dir = '/home/parkinsonjl/mnt/sharedStorage/TestDiffusiveTimescale/';
end

% folders = dir(output_dir);
% 
% 
% for i=1:length(folders)
%     if strcmp(folders(i).name, '.') || strcmp(folders(i).name, '..')
%         continue
%     end
%     fprintf('Processing %s \n', folders(i).name);
%    
%     processSpecificFolder(fullfile(output_dir, folders(i).name), frames);
%     
% end

 processSpecificFolder(output_dir, frames, figure_output_dir);
 
end


function processSpecificFolder(output_dir, frames, figure_output_dir)



savePNG = true;
pcolorShadingInterp = false;
close all; 

domWidth = 1.0;
domHeight = 1.0;
plotHeight = 0.5; % plot top 0.5 of domain

plotAspectRatio = domWidth/plotHeight;

plotScreenWidth = 500;
plotScreenHeight = plotScreenWidth/plotAspectRatio;

% Start plotting stuff
h = figure();
h.Position = [200 200 plotScreenWidth+70 plotScreenHeight*length(frames)];

% Common for all plots
topFraction = plotHeight;
fluidVelScale = 200;
maxPsi = 20;
minPsi = -maxPsi; 
dPsi = (maxPsi-minPsi)/40; % need even number of divisions
v = linspace(minPsi, maxPsi, 40);%minPsi+0.5*dPsi:dPsi:maxPsi-0.5*dPsi; % Avoid the zero contour
SlLims = [-1.0, -0.0];

% Setup axes based on number of plotfiles
axPositions = {};
allAxes = {};

pcolorField = 'S'; %porosity, Sl

numFrames = length(frames);
pltWidth = 0.75;
pltHeight = 0.8/numFrames;

for i=1:length(frames)
   axPositions{end+1} = [0.07, 0.06 + (numFrames-i)*(pltHeight+0.06), pltWidth, pltHeight];
end

% Get all frames in folder
[actual_plot_prefix, allFrames] = getFrames(output_dir, '/');

% Load plotfile
for frame_i=1:length(frames)
    axPos = axPositions{frame_i}; %[0.08 0.08 0.8 0.75];
    
    %ml = getFinalPlotFile(output_dir);
    thisFrame = frames(frame_i);
    if thisFrame == -1
        ml = getFinalPlotFile(output_dir);
    else
        ml =  MushyLayerOutput(2, thisFrame, output_dir, actual_plot_prefix, true);
        
    end
    
    
    mlSmooth =  MushyLayerOutput(2, thisFrame, output_dir, actual_plot_prefix, true, 'bicubic');

    if isempty(isprop(ml, 'levelArray'))
        fprintf('Folder does not contain any plot files \n');
        return;
    end


    % Get data to plot
    [X,Y] = ml.grid();
    porosity = ml.dataForComp(ml.components.Porosity).';
    Sl = ml.dataForComp(ml.components.Liquidconcentration).';
    S = ml.dataForComp(ml.components.Bulkconcentration).';
    U = ml.dataForComp(ml.components.xDarcyvelocity).';
    V = ml.dataForComp(ml.components.yDarcyvelocity).';
    Streamfunction = ml.dataForComp(ml.components.streamfunction).';
    
    % Use level 1 data so the contours are smooth
    lev1 = ml.levelArray(1);
    Streamfunction1 = lev1.dataForComp(ml.components.streamfunction);
    
    U1 = lev1.dataForComp(ml.components.xDarcyvelocity);
    V1 = lev1.dataForComp(ml.components.yDarcyvelocity);
    [X1,Y1] = lev1.grid();
    
    porositySmooth = mlSmooth.dataForComp(ml.components.Porosity).';
    [Xsmooth,Ysmooth] = mlSmooth.grid();
    
    
    
    % Now trim data
    numypts = size(Y, 1);
    ypts = round((1-topFraction)*numypts):1:numypts;
    porosity = porosity(ypts, :);
    Sl = Sl(ypts, :);
    S = S(ypts, :);
    U = U(ypts, :); V= V(ypts, :);
    Streamfunction = Streamfunction(ypts, :);
    X=X(ypts, :);
    Y=Y(ypts, :);
    
    numypts_lev1 = size(Y1, 1);
    yptslev1 = round((1-topFraction)*numypts_lev1):1:numypts_lev1;
    Streamfunction1 = Streamfunction1(yptslev1, :);
   % porosity1 =porosity1(yptslev1, :);
    U1 =U1(yptslev1, :); V1=V1(yptslev1, :); X1=X1(yptslev1, :); Y1=Y1(yptslev1, :);
    
    numypts_smooth = size(Ysmooth, 1);
    ypts_smooth = round((1-topFraction)*numypts_smooth):1:numypts_smooth;
    porositySmooth = porositySmooth(ypts_smooth, :);
    Xsmooth = Xsmooth(ypts_smooth, :);
    Ysmooth = Ysmooth(ypts_smooth, :);
   
    % Get limits
    xl = [X(1,1) X(end,end)];
    yl = [Y(1,1) Y(end,end)];
    dx = X(1,2) - X(1,1);
    
    hold on;

    %allAxes{frame_i} = [];

    if frame_i == 1
          axPorosity{frame_i} = gca;
    else
          axPorosity{frame_i} = axes;
    end
  

    % Porosity first
    if strcmp(pcolorField, 'Sl')
        pcolor(X,Y,Sl);
        colormapLims = [-1, 0];
        colormapLabel = 'S_l';
        
        colormap(axPorosity{frame_i}, bluewhitered(257));
        
        cmapTickLabels = {'-1.0', '0.0'};
        cmapTicks = [-1, 0];
        
    elseif strcmp(pcolorField, 'porosity')
        pcolor(X,Y,porosity);
        
        colormap(axPorosity{frame_i}, blues);
        
        colormapLims = [0, 1];
        colormapLabel = '\chi';
        cmapTickLabels = {'0 (solid)', '1 (liquid)'};
        cmapTicks = [0 1];
        
    elseif strcmp(pcolorField, 'S')
         %colormapLims = [-1.5 -0.5];
        colormapLims = [-1.0, 1.0];
        colormapLabel = 'Bulk Concentration, \Theta';
        
        cmapTickLabels = {sprintf('%1.1f', colormapLims(1)-1), sprintf('%1.1f',colormapLims(2)-1)};
        cmapTicks = colormapLims;
        
        %minS = -2.0;
        %deltaS = -2*(minS+1.0);
        %deltaS = max(max(S)) - minS;
        pcolor(X,Y,S+1);
        
        caxis(colormapLims);
        
        % or coolwarmcmap ?
        colormap(axPorosity{frame_i}, bluewhitered(257));
        
    end
    
    if pcolorShadingInterp
    shading interp;
    end
    
  %  caxis(colormapLims);
    
    
    if frame_i == 2
        cbar = colorbar('Location', 'eastoutside');
    end
    
    %caxis(colormapLims);
    

    daspect([ 1 1 1]);
    xlim(xl);
    ylim(yl);
    box on;


    allAxes{frame_i} = [axPorosity{frame_i}];

    plotFlow = false;
    if plotFlow
        axQuiver = axes;

        % Modify data a little to make velocities more similar
        Uplot = sign(U).*sqrt(abs(U));
        Vplot = sign(V).*sqrt(abs(V));

        sx = 1:6:size(X, 1);
        sy = 1:6:size(X, 2);
        quiverc(X(sx,sy),Y(sx,sy),Uplot(sx,sy),Vplot(sx,sy), 0.8); 

        axQuiver.Visible = 'off';

        daspect([ 1 1 1]);
        xlim(xl);
        ylim(yl);

        allAxes{frame_i}(end+1) = axQuiver;
    end
    
    hold off;
    
    plotPorosityContours = true;
    if plotPorosityContours
        axPorosityContours{frame_i} = axes;
        %colormap(axPorosityContours{frame_i}, autumn);
       
        v = 0:0.2:0.99;
        
       [M, c] = contour(Xsmooth,Ysmooth,porositySmooth, v);
       c.LineWidth = 1;
       c.LineColor = 'k';
       axPorosityContours{frame_i}.Visible = 'off';
         
        daspect([1 1 1]);
        xlim(xl);
        ylim(yl);
        
        allAxes{frame_i}(end+1) = axPorosityContours{frame_i};
    end
    
    plotStreamlines = false;
    if plotStreamlines
        axStream{frame_i} = axes;

        % Modify data a little to make velocities more similar

        %maxPsi = max(max(abs(Streamfunction)));
       % minPsi = -maxPsi; 

       % dPsi = (maxPsi-minPsi)/10; % need even number of divisions
       % v = minPsi+0.5*dPsi:dPsi:maxPsi-0.5*dPsi; % Avoid the zero contour

        colormap(axStream{frame_i}, autumn);
        if frame_i == 2
            cVel = streamfunctioncolor(X1,Y1,U1,V1,Streamfunction1,v, fluidVelScale);
        %cbar.Label.String = 'Scaled fluid velocity';
        %cbar.Ticks = [0 1];
        else
             streamfunctioncolor(X1,Y1,U1,V1,Streamfunction1,v, fluidVelScale);
        end

        axStream{frame_i}.Visible = 'off';

        daspect([1 1 1]);
        xlim(xl);
        ylim(yl);

        allAxes{frame_i}(end+1) = axStream{frame_i};

    end

    % Now liquid salinity contours
    plotSl = false;
    if plotSl
        axSl = axes;
        contour(Xbase,Ybase,Sl_base);
        %pcolor(Xbase,Ybase,Sl_base);
        %xlim(xlbase);
        %ylim(ylbase);

        v = -1.5:0.2:-0.05;

        contour(X,Y,Sl_smooth, v);
        xlim(xl);
        ylim(yl);

        axSl.Visible='off';
        daspect([1 1 1]);

        %cSl = colorbar();
        allAxes{frame_i}(end+1) = axSl;
    end


    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Level outlines
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
   
    if length(ml.levelArray) > 1
        
        axPolyshapes = allAxes{frame_i}(end);
        
        hold on;

         if exist('polyshape')

        domainBox = polyshape([xl(1) xl(1) xl(end) xl(end)], [yl(1) yl(end) yl(end) yl(1)]);

        meshes = ml.getMeshes(domainBox);

        % Plot meshes 
        pshape = meshes(2);
        pshape = trimShape(pshape);
        
        mesh2 = plot(axPolyshapes, pshape); % meshes on chombo level 1

        %plot(domainBox);
        mesh2.FaceAlpha = 0;
        mesh2.EdgeColor = [1 0 1]; % magenta
        mesh2.LineWidth = 2;

        if length(meshes) > 2
            pshape = meshes(3);
            pshape = trimShape(pshape);
            
            mesh3 = plot(axPolyshapes, pshape); % meshes on chombo level 2
            mesh3.FaceAlpha = 0;
            mesh3.EdgeColor = [0 1 1]; %cyan
            mesh3.LineWidth = 2;
        end
        
        else
        fprintf('Not plotting level outlines as polyshape cannot be found \n');
         end
    
        %axPolyshapes.Visible='off';
        %daspect([1 1 1]);
        %allAxes{frame_i}(end+1) = axPolyshapes;
        
        hold off;
        
    end
   

    %axSl.Position = axPos;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Colorbars
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if frame_i == 2
        cbar.Label.String = colormapLabel;
        cbar.Label.Position(1) = cbar.Label.Position(1) -1.0;
        
        cbar.Position(2) = cbar.Position(2)-0.1;
        cbar.Position(4) = cbar.Position(4)+0.2;
        %cPorosity.Label.Position = [0.5 1.2];
        %cPorosity.Label.Rotation = 0;
        cbar.Ticks = cmapTicks;
        cbar.Limits = colormapLims;
        
        cbar.TickLabels = cmapTickLabels;
    end
%     elseif frame_i == 2
%     
%         cVel.Label.String = 'Fluid velocity';
%         cVel.Ticks = [0 1.0];
%        % maxU = sqrt(max(max(abs(U.^2 + V.^2))));
%         cVel.TickLabels = {'0', sprintf('%.1e', fluidVelScale)};
%         cVel.Label.Position = [1.0 0.5];
%     
%     end

    %cPorosity.Label.Rotation = 0;
        %cSl.Label.String = 'S_l';
    %cSl.Position = [0.85 cSl.Position(2) cSl.Position(3) cSl.Position(4)];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Link axes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 2:length(allAxes{frame_i})
        linkaxes([allAxes{frame_i}(1), allAxes{frame_i}(i)]);
    end

    set(allAxes{frame_i}(1), 'Layer', 'top');


    for i = 1:length(allAxes{frame_i})
        allAxes{frame_i}(i).Position = axPos;
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Axis labels
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %  box on;

    allAxes{frame_i}(1).Visible = 'on';
    axLabels{frame_i} = allAxes{frame_i}(1);

    %daspect([ 1 1 1]);

    %xlim(xl);
    %ylim(yl);
     format = '%1.1f';
     
 
    
    if frame_i == numFrames

     xlab = xlabel(allAxes{frame_i}(1), '$x$');
     xlab.Position(2) = xlab.Position(2) + 0.05;
     axLabels{frame_i}.XTick = xl;
     axLabels{frame_i}.XTickLabels = {sprintf(format, xl(1)-dx/4), sprintf(format, xl(2)+dx/4)};
    
    else
        axLabels{frame_i}.XTick = [];
    end
    
    ylab{frame_i} = ylabel(allAxes{frame_i}(1), '$z$');
    ylab{frame_i}.Position(1) = ylab{frame_i}.Position(1) + 0.05;
    axLabels{frame_i}.YTick = yl;

    axLabels{frame_i}.YTickLabels = {sprintf(format, yl(1)-dx/4), sprintf(format, yl(2)+dx/4)};
    
    
    
    % Add a/b/c and time labels
    
    textLabels = {'(a)', '(b)', '(c)'};
    thisLabel = [textLabels{frame_i}, ' $t = ',sprintf('%1.4f', ml.t) , '$'];
    text(0.02, 0.55, thisLabel, 'FontSize', 16);
    
end    

    % h = gcf;
    h.InvertHardcopy = 'off';
    h.Color = 'white';
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    filename = fullfile(figure_output_dir, 'Fig7FixedChillSimulation');
    
    if pcolorShadingInterp
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

    

    
end

% Trim polyshapes to stay within domain
function pshape_trimmed = trimShape(pshape)

dx = 1/128;
dz = 1/128;
domainPolyshape = polyshape([dx dx 1-dx 1-dx], [dz  1-dz 1-dz dz]);

%pshape_trimmed = intersect(pshape, domainPolyshape);

for i = 1:length(pshape)
     pshape_trimmed(i) = intersect(pshape(i), domainPolyshape);
end
stop =0;

end