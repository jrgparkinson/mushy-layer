function [xzoom,yzoom] = plotSingleChannelReflected(output, axisExtent, drawColorbar, drawLabels, toffset, ...
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
    
    if mushyZoom
        % Find y limits to just include mushy layer
        %[~, y_i_min] = find(chi
        temp = 0;
        [~, y_mush_i] = find(chi<1.0);
        [~, y_sol_i] = find(chi<1e-5);
        
        y_min = round(min(y_mush_i)*0.9);
        y_max = round(min(y_sol_i)*1.1);
        
        % Make sure we don't go outside bounds
        y_min = max(y_min, 1);
        y_max = min(y_max, length(y));
        
        yzoom = y(y_min:y_max);
        streamfunction = streamfunction(:, y_min:y_max);
        chi = chi(:, y_min:y_max);
        Sl = Sl(:, y_min:y_max);
        
    else
        yzoom = y;
    end
    
    
    if length(lims) > 0
        yl = lims(2,: );
        y_min = find(y>yl(1));
        y_min = min(y_min);
        
        y_max = find(y<yl(2));
        y_max = max(y_max);
        
        yzoom = y(y_min:y_max);
        streamfunction = streamfunction(:, y_min:y_max);
        chi = chi(:, y_min:y_max);
        Sl = Sl(:, y_min:y_max);
    end
    
    
    if setOriginToZero(1)
       xzoom = xzoom - xzoom(1); 
    end
    
    if setOriginToZero(2)
       yzoom = yzoom - yzoom(1); 
    end
    
    [X, Y] = meshgrid(xzoom, yzoom);
    
   
    t = output.t - toffset;
end


standardLinewidth = 1.0;


% Now reflect data
X2 = [-fliplr(X) X];
Y2 = [Y Y];
Sl2 = [Sl; NaN*flipud(Sl)];
streamfunction2 = [NaN*streamfunction; flipud(streamfunction)];
chi2 = [chi; flipud(chi)];
chiRight = [NaN*chi; flipud(chi)].';

X = X2;
Y = Y2;
Sl = Sl2.';
streamfunction = streamfunction2.';
chi = chi2.';

xzoom = [-fliplr(xzoom) xzoom];

fprintf('Max abs(streamfunction): %1.3f \n', max(max(abs(streamfunction))));

%pcolor(X,Y,chi); colorbar();

axPorosity = axes;
% This should be the 'blues' colormap from visit
blues =  makeColorMap([8 48 107]/255, [104, 172, 212]/255, [1 1 1]);

colormap(axPorosity,flipud(blues));

%colormap(axPorosity, flipud(parula));
caxis(axPorosity, [0 1]);
%imagesc(xzoom, yzoom, chi.');
pcolor(X,Y,chi);
%set(axPorosity,'dataAspectRatio',[1 1 1])
%set(axPorosity,'Ydir','normal');

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
    lScale = 10; % 10 dm in a metre
    xlab = xlabel('$x (cm)$');
    pos = get(xlab, 'Position');
    
    axEx = get(axPorosity, 'Position');
    
    %set(xlab, 'Position', pos+[0 0.08 0]);
    set(xlab, 'Units', 'normalized');
    set(xlab, 'Position', [0.5 axEx(2)-0.18 0.0]);
    
    ylab = ylabel('$z (cm)$'); oldPos = ylab.Position;
    %ylab.Position = [oldPos(1) + 0.1 oldPos(2)];
    ylab.Position = [oldPos(1) + 0.0 oldPos(2)];
    %set(ylab, 'Rotation', 0);
    %set(ylab, 'Units', 'normalized');
    %set(ylab, 'Position', [axEx(1)-0.18 0.5 0.0]);
    
    %%pos = get(ylab, 'Position');
    %%set(ylab, 'Position', [-0.08 0.15 0]);
  
    yTop = 0.0;
    yBottom = min(yzoom)-max(yzoom);
    
    
    axPorosity.XTick = [min(min(X)) max(max(X))];
    axPorosity.YTick = [yBottom yTop];
    
    formatStr = ['%1.',num2str(axisLabelPrecision),'f'];
    
    xminStr = sprintf(formatStr, min(min(X))*lScale);
    xmaxStr = sprintf(formatStr, max(max(X))*lScale);
    axPorosity.XTickLabels = {xminStr, xmaxStr};
    
   
    
    yminStr = sprintf(formatStr, yBottom*lScale);
    ymaxStr = sprintf(formatStr, yTop*lScale);
    axPorosity.YTickLabels = {yminStr, ymaxStr};
    
    set(axPorosity, 'TickLength',[0 0])
    
else
    axPorosity.XTick = [];
    axPorosity.YTick = [];

    set(axPorosity, 'TickLength',[0 0])

end

set(axPorosity,'dataAspectRatio',[1 1 1])

% HACK
  yTop = 0.0;
    yBottom = min(yzoom)-max(yzoom);
    
axPorosity.YTick = [min(yzoom) max(yzoom)];
    
    formatStr = ['%1.',num2str(axisLabelPrecision),'f'];

    yminStr = sprintf(formatStr, yBottom);
    ymaxStr = sprintf(formatStr,yTop);
    axPorosity.YTickLabels = {yminStr, ymaxStr};
    
    set(axPorosity, 'TickLength',[0 0])
% END HACK


set(axPorosity, 'position', axisExtent);

axChiContour = axes;



    contourValsChi = linspace(0, 1, 10); %num porosity contours
  
   chiColorMap = makeColorMap([0 0 0], [0 0 0]);
    colormap(axChiContour, chiColorMap);
    caxis(axChiContour, [0 1]);
    [Cchi, hchi] = contour(X, Y, chiRight, contourValsChi);
    hchi.LineWidth = 1.0;

set(axChiContour,'dataAspectRatio',[1 1 1])
set(axChiContour, 'position', axisExtent);
set(axChiContour,'TickLength',[0 0])

axis(axChiContour, 'off');
axChiContour.Visible = 'off';


axSl = axes;

if drawSl
    contourValsSl = linspace(-0.98, 0, 10);

    %colormap(axSl, makeColorMap([0 0 0] , [230, 175, 0]/255, [200 0 0]/255));
    salinityColormap = makeColorMap([0 0 0] ,  [200 0 0]/255);
    salinityColormap = makeColorMap([1 1 1] ,  [0.5 0.5 0.5], [200 0 0]/255);
    %salinityColormap = makeColorMap([200 0 0]/255, [200 0 0]/255);
    salinityColormap = makeColorMap([1 0 1], [1 0 1]);
    colormap(axSl, salinityColormap);
    caxis(axSl, [-1 0]);
    [CSl, hSl] = contour(X, Y, Sl, contourValsSl);
    hSl.LineWidth = standardLinewidth;

end

set(axSl,'dataAspectRatio',[1 1 1])
set(axSl, 'position', axisExtent);
set(axSl,'TickLength',[0 0])

axis(axSl, 'off');
axSl.Visible = 'off';


axPsi = axes;
if maxPsi > 0
    %maxPsi = 0.2; %max(max(abs(streamfunction)));
    psiVals = linspace(-maxPsi,maxPsi,2*Npsi);
    %psiColormap = makeColorMap([1 1 1], [1 1 1]);
    psiColormap = makeColorMap([1 0 1], [1 0 1]);
    gray = 0.5;
    psiColormap = makeColorMap([gray gray gray], [gray gray gray]);
    
    colormap(axPsi, psiColormap);
    %[cPsi,hPsi] = contour(X,Y,streamfunction.', psiVals);
    %hPsi.LineWidth = standardLinewidth;
    % Want to avoid the psi=0 contour which is always a bit wierd
    minPsi = maxPsi/(10*Npsi);
    contourValsNeg = linspace(-maxPsi,-minPsi,Npsi);
    contourValsPos = linspace(minPsi,maxPsi,Npsi);
    
    %spacing = contourValsPos(end)-contourValsPos(end-1);
    %fineSpacing = minPsi/3;
    %Nsmall = round(contourValsPos(2)/fineSpacing);
    
    %contourValsNegSmall = linspace(contourValsNeg(end-1),fineSpacing,Nsmall);
    %contourValsPosSmall = linspace(fineSpacing,contourValsPos(2),Nsmall);
    
    hold on;
    [CconNeg, hconNeg] = contour(X, Y, streamfunction, contourValsNeg);
    [CconPos, hconPos] = contour(X, Y, streamfunction, contourValsPos);

    hconNeg.LineWidth = standardLinewidth;
    hconPos.LineWidth = standardLinewidth;
    hconPos.LineStyle = '-';
    
    %[CconNegSmall, hconNegSmall] = contour(X, Y, streamfunction, contourValsNegSmall);
    %[CconPosSmall, hconPosSmall] = contour(X, Y, streamfunction, contourValsPosSmall);
    % hconNegSmall.LineWidth = standardLinewidth;
    %hconPosSmall.LineWidth = standardLinewidth;
    %hconPosSmall.LineStyle = '--';
    %hconNegSmall.LineStyle = '--';
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





end
