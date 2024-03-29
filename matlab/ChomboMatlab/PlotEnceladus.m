%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A script to plot data from a 2D chombo simulation, or a series of
% simulations.
% Either make 2D plots, or plot the temperature and temperature flux at the
% top boundary of the domain.
% Make sure you update the options in the first c. 50 lines to do what you
% want!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PlotEnceladus
close all;

% Specify location of plot files you want to use
% Loads the file at
%     {output_dir}{actual_plot_prefix}{frame}.2d.hdf5
% where {thisFrame} is padded with 0's as needed
% so in this case,
% /home/parkinsonjl/mushy-layer/execSubcycle/enceladus/plt000500.2d.hdf5
output_dir = '/home/parkinsonjl/mushy-layer/execSubcycle/enceladus/'; % directory containing plot files
actual_plot_prefix = 'plt'; % prefix used to make plot files (all files follow the naming pattern {prefix}{frame}.2d.hdf5)
frames = [100, 500]; % frame(s) to plot

% Specify where to put the axes for the different plots
% For just 1 frame:
if length(frames) == 1
    axPositions = [[0.1, 0.1, 0.6, 0.75]];
elseif length(frames) == 2
    % For two frames:
    axPositions = [[0.1 0.55 0.6 0.3],
        [0.1 0.1 0.6 0.3]];
else
    % If using more than two frames, specify the positions yourself!
    % (and uncomment the two lines below)
    fprintf('Error - axis positions not specified \n');
    return
    
end

% Decide which plots to make
% Can either make 2D plots of the porosity and flow fields, or
% Just plot the temperature and temperature flux at the top boundary
% Only set one of these to 'true'
make_2d_plot = false;
make_top_boundary_plot = true;

if make_2d_plot && make_top_boundary_plot
    print('Error - cannot do 2D plot and top boundary plot simultaneously');
    return
end


% Options for the 2D plot
options.pcolorShadingInterp = true; % set shading interp on for pcolor plots
options.topFraction = 1.0; % How much of the domain to plot. Set to 1.0 to plot the whole domain, 0.5 to just plot the top half, etc.
options.fluidVelScale = 200; % For scaling fluid velocity arrows
options.pcolorField = 'porosity'; % Which field to plot on the color scale (other options: Sl (liquid salinity), S (bulk salinity))
options.timeTitle = true; % if true, add a title giving the time of the frame


% Make the figure with the size you want
h = figure();
h.Position = [200 200 1000 600];


for frame_i = 1:length(frames)
    
    thisFrame = frames(frame_i);
    options.axPos = axPositions(frame_i, :);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %  Load data
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    if thisFrame == -1
        ml = getFinalPlotFile(output_dir);
    else
        ml = MushyLayerOutput(2, thisFrame, output_dir, actual_plot_prefix, true);
        
    end
    
    
    if isempty(isprop(ml, 'levelArray')) || length(ml.levelArray) == 0
        
        options.maxSearchFrame = 100000; 
        if isfield(options, 'maxSearchFrame')
            while isempty(isprop(ml, 'levelArray')) || length(ml.levelArray) == 0 ...
                    && thisFrame < options.maxSearchFrame
                thisFrame = thisFrame + 1;
                fprintf('  Trying frame %d \n', thisFrame);
                ml = MushyLayerOutput(2, thisFrame, output_dir, actual_plot_prefix, true);
                
            end
        else
            fprintf('Could not find plot file \n');
            return;
        end
    end
    
    mlSmooth =  MushyLayerOutput(2, thisFrame, output_dir, actual_plot_prefix, true, 'bicubic');
    
    % Get matlab version
    [v, d] = version;
    year = str2num(d(end-4:end));
    if year > 2016
        recentMatlab = true;
    else
        recentMatlab = false;
    end
    
    % Get data to plot
    [X,Y] = ml.grid();
    porosity = ml.dataForComp(ml.components.Porosity).';
    Sl = ml.dataForComp(ml.components.Liquidconcentration).';
    S  = ml.dataForComp(ml.components.Bulkconcentration).';
    U  = ml.dataForComp(ml.components.xDarcyvelocity).';
    V  = ml.dataForComp(ml.components.yDarcyvelocity).';
    Streamfunction = ml.dataForComp(ml.components.streamfunction).';
    
    % Use level 1 data so the contours are smooth
    lev1 = ml.levelArray(1);
    Streamfunction1 = lev1.dataForComp(ml.components.streamfunction);
    
    U1 = lev1.dataForComp(ml.components.xDarcyvelocity);
    V1 = lev1.dataForComp(ml.components.yDarcyvelocity);
    [X1,Y1] = lev1.grid();
    
    porositySmooth = mlSmooth.dataForComp(ml.components.Porosity).';
    SlSmooth = mlSmooth.dataForComp(ml.components.Liquidconcentration).';
    TSmooth = mlSmooth.dataForComp(ml.components.Temperature).';
    [Xsmooth,Ysmooth] = mlSmooth.grid();
    
    % Now trim data
    numypts = size(Y, 1);
    min_y_pts = max(1, round((1-options.topFraction)*numypts));
    ypts = min_y_pts:1:numypts;
    porosity = porosity(ypts, :);
    Sl = Sl(ypts, :);
    S = S(ypts, :);
    U = U(ypts, :); V= V(ypts, :);
    Streamfunction = Streamfunction(ypts, :);
    SlSmooth = SlSmooth(ypts, :);
    TSmooth = TSmooth(ypts, :);
    X=X(ypts, :);
    Y=Y(ypts, :);
    
    numypts_lev1 = size(Y1, 1);
    min_point = max(1, round((1-options.topFraction)*numypts_lev1));
    yptslev1 = min_point:1:numypts_lev1;
    Streamfunction1 = Streamfunction1(yptslev1, :);
    
    U1 =U1(yptslev1, :); V1=V1(yptslev1, :); X1=X1(yptslev1, :); Y1=Y1(yptslev1, :);
    
    numypts_smooth = size(Ysmooth, 1);
    min_y_val = max(1, round((1-options.topFraction)*numypts_smooth));
    ypts_smooth = min_y_val:1:numypts_smooth;
    porositySmooth = porositySmooth(ypts_smooth, :);
    Xsmooth = Xsmooth(ypts_smooth, :);
    Ysmooth = Ysmooth(ypts_smooth, :);
    
    % Get limits
    xl = [X(1,1) X(end,end)];
    yl = [Y(1,1) Y(end,end)];
    dx = X(1,2) - X(1,1);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %  Make 2D plot
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    if make_2d_plot
        
        
        
        hold on;
        
        %allAxes{frame_i} = [];
        
        if frame_i == 1
            axPorosity{frame_i} = gca;
        else
            axPorosity{frame_i} = axes;
        end
        
        
        % Porosity first
        if strcmp(options.pcolorField, 'Sl')
            pcolor(X,Y,Sl);
            colormapLims = [-1, 0];
            colormapLabel = 'S_l';
            
            colormap(axPorosity{frame_i}, bluewhitered(257));
            
            cmapTickLabels = {'-1.0', '0.0'};
            cmapTicks = [-1, 0];
            
        elseif strcmp(options.pcolorField, 'porosity')
            pcolor(X,Y,porosity);
            
            colormap(axPorosity{frame_i}, blues);
            
            colormapLims = [0, 1];
            colormapLabel = '\chi';
            cmapTickLabels = {'0 (solid)', '1 (liquid)'};
            cmapTicks = [0 1];
            
        elseif strcmp(options.pcolorField, 'S')
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
        
        if options.pcolorShadingInterp
            shading interp;
        end
        
        %  caxis(colormapLims);
        
        
        
        cbar = colorbar('Location', 'eastoutside');
        
        
        
        %caxis(colormapLims);
        daspect([ 1 1 1]);
        xlim(xl);
        ylim(yl);
        box on;
        
        
        allAxes{frame_i} = [axPorosity{frame_i}];
        
        plotFlow = true;
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
        
        plotPorosityContours = false;
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
            
            v = -1.5:0.2:-0.05;
            
            contour(X,Y,SlSmooth, v);
            xlim(xl);
            ylim(yl);
            
            axSl.Visible='off';
            daspect([1 1 1]);
            
            %cSl = colorbar();
            allAxes{frame_i}(end+1) = axSl;
        end
        
        
        plotTemperature = false;
        if plotTemperature
            axT = axes;
            
            v = -1.5:0.2:1.5;
            
            contour(X,Y,TSmooth, v);
            xlim(xl);
            ylim(yl);
            
            colormap(axT, jet);
            
            axT.Visible='off';
            daspect([1 1 1]);
            
            %cSl = colorbar();
            allAxes{frame_i}(end+1) = axT;
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
                % [1 0 1]; % magenta
                mesh2.EdgeColor = [0 1 0]; % green
                mesh2.LineWidth = 2;
                
                if length(meshes) > 2
                    pshape = meshes(3);
                    pshape = trimShape(pshape);
                    
                    mesh3 = plot(axPolyshapes, pshape); % meshes on chombo level 2
                    mesh3.FaceAlpha = 0;
                    % [0 1 1]; %cyan
                    mesh3.EdgeColor = [1 0 1]; % magenta
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
        
        cbar.Label.String = colormapLabel;
        
        
        %         if isfield(options, 'cbarPos')
        %
        %         else
        %             cbar.Position(1) = cbar.Position(1)+0.05;
        %
        %             cbar.Position(2) = cbar.Position(2)-0.1;
        %             cbar.Position(4) = cbar.Position(4)*1.0;
        %             cbar.Position(3) = cbar.Position(3)*0.6;
        %
        %             if recentMatlab
        %
        %             cbar.Label.Position(1) = cbar.Label.Position(1) + 0.0;
        %
        %
        %             else
        %                 cbar.Label.Position(1) = cbar.Label.Position(1) + 1.0;
        %             end
        %
        %         end
        
        %cPorosity.Label.Position = [0.5 1.2];
        %cPorosity.Label.Rotation = 0;
        cbar.Ticks = cmapTicks;
        cbar.Limits = colormapLims;
        
        cbar.TickLabels = cmapTickLabels;
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Link axes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i = 2:length(allAxes{frame_i})
            linkaxes([allAxes{frame_i}(1), allAxes{frame_i}(i)]);
        end
        
        set(allAxes{frame_i}(1), 'Layer', 'top');
        
        
        for i = 1:length(allAxes{frame_i})
            allAxes{frame_i}(i).Position = options.axPos;
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Axis labels
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        allAxes{frame_i}(1).Visible = 'on';
        axLabels{frame_i} = allAxes{frame_i}(1);
        
        format = '%1.1f';
        
        %if frame_i == options.numFrames
            
            xlab = xlabel(allAxes{frame_i}(1), '$x$');
            xlab.Position(2) = xlab.Position(2) + 0.05;
            axLabels{frame_i}.XTick = xl;
            axLabels{frame_i}.XTickLabels = {sprintf(format, xl(1)-dx/4), sprintf(format, xl(2)+dx/4)};
            
       % else
       %     axLabels{frame_i}.XTick = [];
       % end
        
        ylab{frame_i} = ylabel(allAxes{frame_i}(1), '$z$');
        ylab{frame_i}.Position(1) = ylab{frame_i}.Position(1) + 0.05;
        axLabels{frame_i}.YTick = yl;
        
        axLabels{frame_i}.YTickLabels = {sprintf(format, yl(1)-dx/4), sprintf(format, yl(2)+dx/4)};
        
        % Add a/b/c and time labels
            thisLabel = [' $t = ',sprintf('%1.4f', ml.t) , '$'];
            if options.timeTitle
                title(allAxes{frame_i}(1), thisLabel);
            %else
            %    text(0.02, 0.55, thisLabel, 'FontSize', textFontSize);
            end
        
        % allAxes{frame_i}(1).FontName = 'Times'
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Get temperature and temperature flux at top boundary
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if make_top_boundary_plot
        
        
        y_top = Y(end, end);
        
        x_vals = X(1, :);
        y_vals = y_top + 0*x_vals;
        TemperatureOnTopBoundary = interp2(X,Y,TSmooth,x_vals, y_vals);
        
        [dTdx, dTdy] = gradient(TSmooth, dx);
        TFluxOnTopBoundary = interp2(X, Y, dTdx, x_vals, y_vals);
        
        if frame_i == 1
            ax = gca;
        else
            ax = axes;
        end
        ax.Position = options.axPos;
        
        %figure();
        hold on;
        plot(x_vals, TemperatureOnTopBoundary);
        ylabel('$T$');
        yyaxis right;
        plot(x_vals, TFluxOnTopBoundary);
        hold off;
        xlabel('$x$');
        ylabel('$dT/dy$');
        
        title(sprintf('$t = %f$', ml.t ));
        
        %legend({'$T$', '$dT/dy$'});
        box on;
        
        
    end
    
    
end % end loop over frames

% return success=true if we've made it to the end of the function

success = true;

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
