function [xzoom,yzoom] = plotChannel(output, axisExtent, drawColorbar, drawLabels, toffset, ...
    drawTime, drawSl, maxPsi, Npsi, mushyZoom, axisLabelPrecision, setOriginToZero, lims)

if nargin < 13
    lims = [];
end

if nargin < 12
    axisLabelPrecision = 1;
    setOriginToZero = [false false];
end

if nargin < 10
    mushyZoom = false;
end

if nargin < 9
    drawTime = false;
    drawSl = false;
    maxPsi = -1; % less than 0 means don't add psi
    Npsi = 5;
end

if nargin < 3
    drawColorbar = false;
    drawLabels = false;
    toffset = 0.0;
end

if length(output.levelArray) == 0
    xzoom = []; yzoom =[];
    return;
end

TESTING = false;

if TESTING
    xzoom = linspace(0,4,100);
    yzoom = linspace(0,1.0,100);
    
    [Y, X] = meshgrid(xzoom, yzoom);
    
    chi = sin(3.14*X.*X);
    streamfunction = 0.2*Y.*Y.*sin(X*3.14);
    Sl = -Y;
    
    t = 0.123556;
    
else
    
    porosity = output.dataForComp(output.components.Porosity);
    streamfunction = output.getStreamfunction(100, 1);
    
    xvel = output.dataForComp(output.components.xAdvectionvelocity);
    yvel = output.dataForComp(output.components.yAdvectionvelocity);
    
    Sl = output.dataForComp(output.components.Liquidconcentration);
    T = output.dataForComp(output.components.Temperature);
    
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
    streamfunction = streamfunction.';
    
    chi = porosity;
    U = xvel.';
    V = yvel.';
    %psi = streamfunction;
    
    xzoom = x;
    
    % Now zoom in on channel
    
    extraPadding = 3;
    % Find y limits to just include mushy layer
    [~, y_mush_i] = find(chi<1.0);
    [~, y_sol_i] = find(chi<1e-4);
    [~, y_liq_i] = find(chi==1.0);
    
    y_min = round(min(y_mush_i)) - extraPadding;
    %y_max = round(min(y_sol_i))-1;
    y_max = round(max(y_liq_i)) + extraPadding;
    
    
    % Make sure we don't go outside bounds
    y_min = max(y_min, 1);
    y_max = min(y_max, length(y));
    
    yzoom = y(y_min:y_max);
    
    % Find x limits to just include channel
    channelSide = +1;
    x_chanWall = find(chi(:,min(y_mush_i))<1.0);
    
    
    
    if channelSide > 0
        xzoom = x(x_chanWall-extraPadding:end);
        
        x_max = length(x);
        x_min = x_chanWall-extraPadding;
    else
        xzoom = x(1:x_chanWall+extraPadding);
        
        x_min = 1;
        x_max = x_chanWall+extraPadding;
    end
    
    
    
   % streamfunction = streamfunction(x_min:x_max, y_min:y_max);
   % chi = chi(x_min:x_max, y_min:y_max);
   % Sl = Sl(x_min:x_max, y_min:y_max);
   % T = T(x_min:x_max, y_min:y_max);
    
    
    if length(lims) > 0
        yl = lims(2,: );
        if yl(1) ~= yl(2)
            
        y_min = find(y>yl(1));
        y_min = min(y_min);
        
        y_max = find(y<yl(2));
        y_max = max(y_max);
        
        yzoom = y(y_min:y_max);
        
        end
        
        xl = lims(1, :);
        if xl(1) ~= xl(2)
            
            x_min = find(x>xl(1));
        x_min = min(x_min);
        
        x_max = find(x<xl(2));
        x_max = max(x_max);
        
        xzoom = x(x_min:x_max);
        
        end
        
        
        
       
    end
    
    
     streamfunction = streamfunction(x_min:x_max, y_min:y_max);
        chi = chi(x_min:x_max, y_min:y_max);
        Sl = Sl(x_min:x_max, y_min:y_max);
        T = T(x_min:x_max, y_min:y_max);
%     
    
    if setOriginToZero(1)
        xzoom = xzoom - xzoom(1);
    end
    
    if setOriginToZero(2)
        yzoom = yzoom - yzoom(1);
    end
    
    [X, Y] = meshgrid(xzoom, yzoom);
    
    
    t = output.t - toffset;
end


%%%%%%%%%%%%%%%%%%%%%%%%%
% Porosity
%%%%%%%%%%%%%%%%%%%%%%%%%

axPorosity = axes;

%make smooth chi:
dx = xzoom(2)-xzoom(1);
dy = yzoom(2)-yzoom(1);

zoomFactor = 10;
[xnew, ynew] = meshgrid(xzoom(1):dx/zoomFactor:xzoom(end), yzoom(1):dy/zoomFactor:yzoom(end));
chiZoom = interp2(xzoom,yzoom,chi.',xnew,ynew, 'spline');


contourValsPorosity = [0.999 0.999];

hold on
[CChi, hChi] = contour(X, Y, chi.', contourValsPorosity);
%[CChiSmooth, hChiSmooth] = contour(xnew, ynew, chiZoom, contourValsPorosity);
hold off

hChi.LineWidth = 4.0;
hChi.LineColor = 'r';

%hChiSmooth.LineWidth = 4.0;
%hChiSmooth.LineColor = 'r';



set(axPorosity,'dataAspectRatio',[1 1 1])
set(axPorosity, 'Layer', 'top'); % put box border on top
box on;

% done later
%set(axPorosity, 'position', axisExtent);

if drawColorbar
    cbar = colorbar('Location', 'eastoutside');
    cbar.Label.String = '\chi';
    cbar.Label.Rotation = 0.0;
    cbar.Label.Position = [1.6 0.6 0]; %cbar.Label.Position + [0.02 0.8 0.0];
    cbar.Ticks = [0.01 0.99];
    cbar.TickLabels = {'0', '1'};
end

if drawLabels
    xlab = xlabel('$x$');
    pos = get(xlab, 'Position');
    
    axEx = get(axPorosity, 'Position');
    
    %set(xlab, 'Position', pos+[0 0.08 0]);
    set(xlab, 'Units', 'normalized');
    set(xlab, 'Position', [0.5 axEx(2)-0.2 0.0]);
    
    ylab = ylabel('$z$');
    set(ylab, 'Rotation', 0);
    %pos = get(ylab, 'Position');
    %set(ylab, 'Position', [-0.08 0.15 0]);
    set(ylab, 'Units', 'normalized');
    set(ylab, 'Position', [axEx(1)-0.25 0.5 0.0]);
    
    axPorosity.XTick = [min(xzoom) max(xzoom)];
    axPorosity.YTick = [min(yzoom) max(yzoom)];
    
    formatStr = ['%1.',num2str(axisLabelPrecision),'f'];
    
    xminStr = sprintf(formatStr, min(xzoom));
    xmaxStr = sprintf(formatStr, max(xzoom));
    axPorosity.XTickLabels = {xminStr, xmaxStr};
    
    yminStr = sprintf(formatStr, min(yzoom));
    ymaxStr = sprintf(formatStr, max(yzoom));
    axPorosity.YTickLabels = {yminStr, ymaxStr};
    
    set(axPorosity, 'TickLength',[0 0])
    
else
    axPorosity.XTick = [];
    axPorosity.YTick = [];
    
    set(axPorosity, 'TickLength',[0 0])
    
end


% HACK
axPorosity.YTick = [min(yzoom) max(yzoom)];

formatStr = ['%1.',num2str(axisLabelPrecision),'f'];

yminStr = sprintf(formatStr, min(yzoom));
ymaxStr = sprintf(formatStr, max(yzoom));
axPorosity.YTickLabels = {yminStr, ymaxStr};

set(axPorosity, 'TickLength',[0 0])
% END HACK


set(axPorosity, 'position', axisExtent);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sl
%%%%%%%%%%%%%%%%%%%%%%%%%%%

axSl = axes;

if drawSl
    
    minSl = min(min(Sl));
    
    minRound = round(minSl, 1);
    if minRound < minSl
        minRound = minRound + 0.1;
    end
    
    numPoints = (-minRound/0.1) + 1;
    
    contourValsSl = linspace(minRound, 0, numPoints);
    contourValsSl = linspace(-1, 0, 11);
    
    %colormap(axSl, makeColorMap([0 0 0] , [230, 175, 0]/255, [200 0 0]/255));
   % salinityColormap = makeColorMap([0 0 0] ,  [200 0 0]/255);
   % salinityColormap = makeColorMap([1 1 1] ,  [0.5 0.5 0.5], [200 0 0]/255);
   % salinityColormap = makeColorMap([200 0 0]/255, [200 0 0]/255);
   % colormap(axSl, salinityColormap);
    
    colormap(axSl, parula);
    %colormap(axSl, flipud(autumn));
    
    
    %caxis(axSl, [-1 0]);
    [CSl, hSl] = contour(X, Y, Sl.', contourValsSl);
    
    caxis(axSl, [-1 0]); % Needs to come after contour()
    hSl.LineWidth = 2.0;
    hSl.LineStyle = '--';
    
    
end

set(axSl,'dataAspectRatio',[1 1 1])
set(axSl, 'position', axisExtent);
set(axSl,'TickLength',[0 0])

axis(axSl, 'off');
axSl.Visible = 'off';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Temperature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axT = axes;

colormap(axT, coolwarmcmap);

contourValsT = linspace(0, 1.0, 11);
[CT, hT] = contour(X, Y, T.', contourValsT);
caxis([0 1]);
hT.LineWidth = 2.0;
%hT.LineStyle = '--';
set(axT,'dataAspectRatio',[1 1 1])
set(axT, 'position', axisExtent);
set(axT,'TickLength',[0 0])

axis(axT, 'off');
axT.Visible = 'off';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Psi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



axPsi = axes;
if maxPsi > 0
    %maxPsi = 0.2; %max(max(abs(streamfunction)));
    psiVals = linspace(-maxPsi,maxPsi,2*Npsi);
    psiColormap = makeColorMap([1 1 1], [1 1 1]);
    psiColormap = makeColorMap([0 0 0], [0 0 0]);
    
    colormap(axPsi, psiColormap);
    %[cPsi,hPsi] = contour(X,Y,streamfunction.', psiVals);
    %hPsi.LineWidth = 2.0;
    % Want to avoid the psi=0 contour which is always a bit wierd
    minPsi = maxPsi/(10*Npsi);
    contourValsNeg = linspace(-maxPsi,-minPsi,Npsi);
    contourValsPos = linspace(minPsi,maxPsi,Npsi);
    hold on;
    [CconNeg, hconNeg] = contour(X, Y, streamfunction.', contourValsNeg);
    [CconPos, hconPos] = contour(X, Y, streamfunction.', contourValsPos);
    
    hconNeg.LineWidth = 1.0;
    hconPos.LineWidth = 1.0;
    hconPos.LineStyle = '-';
    
    %Make arrows!
    %plotArrows(CconNeg, hconNeg, 'w', 0.6);
    %plotArrows(CconPos, hconPos, 'w', 0.6);
    hold off;
    
    %hPsi.LineStyle = '--';
    set(axPsi,'dataAspectRatio',[1 1 1])
    set(axPsi, 'position', axisExtent);
    set(axPsi,'TickLength',[0 0])
    
end



axis(axPsi, 'off');
axPsi.Visible = 'off';

% Text annotation
if drawTime
    textPos = [axisExtent(1)+0.0045 axisExtent(2)+axisExtent(4)-0.0905 0.1 0.09];
    an = annotation('textbox', textPos, 'String',sprintf('t=%1.3f',t),'LineStyle','none');
else
    fprintf('t=%1.3f\n',t);
end



linkaxes([axPorosity axSl]);

linkaxes([axT axSl]);






end