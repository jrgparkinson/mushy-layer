function [xzoom,yzoom] = plotWideDomain(output, axisExtent, drawColorbar, drawLabels, toffset, ...
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
    drawLabels = [false false];
    toffset = 0.0;
end

lscale = 0.01;
drawPsi = false;

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
    streamfunction = NaN*porosity; %output.dataForComp(output.components.streamfunction);
    
    xvel =  NaN*porosity; %output.dataForComp(output.components.xAdvectionvelocity);
    yvel =  NaN*porosity; %output.dataForComp(output.components.yAdvectionvelocity);
    
    Sl = output.dataForComp(output.components.Liquidconcentration);
    
    [dsldx, dsldz] = gradient(Sl);
    schlieren = dsldx;
    
    schlieren = Sl;
    maxSlSchlieren = -0.9;
    schlieren(schlieren > maxSlSchlieren) = NaN;
    
    [~,mush_j] = find(porosity<1);
    mushBottom = min(mush_j);
    avChiHoriz = nanmean(porosity, 1);
   
    mushBottom = min(find(avChiHoriz<0.8));
    
    schlieren(porosity < 1) = NaN;
    
    
    
    schlieren(:, mushBottom:end) = NaN;
    
    %get rid of outliers
    %mn = nanmean(nanmean(schlieren)); std = nanstd(nanstd(schlieren));
    %schlieren(schlieren < mn-8*std) = NaN;
    %schlieren(schlieren > mn+8*std) = NaN;
    
    %schlieren(schlieren > -0.98) = NaN;
    
    
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
   % streamfunction = streamfunction.';
    
    chi = porosity;
    %U = xvel.';
   % V = yvel.';
    %psi = streamfunction;
    
    xzoom = x;
    
    if mushyZoom
        % Find y limits to just include mushy layer
        %[~, y_i_min] = find(chi
        temp = 0;
        [~, y_mush_i] = find(chi<1.0);
        [~, y_sol_i] = find(chi<1e-5);
        
        y_max = round(min(y_sol_i)*1.1);
        if length(y_max) == 0
            y_max = length(y);
        end
        
        mushDepth = y_max - min(y_mush_i);
        
        y_min = round(min(y_mush_i)-mushDepth);
        
        
        % Make sure we don't go outside bounds
        y_min = max(y_min, 1);
        y_max = min(y_max, length(y));
        
        yzoom = y(y_min:y_max);
        streamfunction = streamfunction(:, y_min:y_max);
        chi = chi(:, y_min:y_max);
        Sl = Sl(:, y_min:y_max);
        schlieren = schlieren(:, y_min:y_max);
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
        schlieren = schlieren(:, y_min:y_max);
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




axPorosity = axes;
% This should be the 'blues' colormap from visit
blues = makeColorMap([8 48 107]/255, [104, 172, 212]/255, [1 1 1]);
parulaBlueEnd = [0.0165    0.4266    0.8786]; %[0.2081    0.1663    0.5292]
bluesParula = makeColorMap(parulaBlueEnd, [1 1 1]); % Trying to fade into parula

colormap(axPorosity, flipud(blues));
%colormap(axPorosity, flipud(parula));
caxis(axPorosity, [0 1]);
imagesc(xzoom, yzoom, chi.');
set(axPorosity,'dataAspectRatio',[1 1 1])
set(axPorosity,'Ydir','normal');

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

formatStr = ['%1.',num2str(axisLabelPrecision),'f'];
 axEx = get(axPorosity, 'Position');
 
if drawLabels(1)
    xlab = xlabel('$x$ (cm)');
    pos = get(xlab, 'Position');
    
   
    set(xlab, 'Units', 'normalized');
    set(xlab, 'Position', [0.5 axEx(2)-0.15 0.0]);

    axPorosity.XTick = [min(xzoom) max(xzoom)];
       
   
    xminStr = sprintf(formatStr, min(xzoom));
    xmaxStr = sprintf(formatStr, max(xzoom));
    axPorosity.XTickLabels = {xminStr, xmaxStr};
    
    set(axPorosity, 'TickLength',[0 0])
else
     axPorosity.XTick = [];
end

% ylabels
if drawLabels(2)
    
     ylab = ylabel('$z$ (cm)');
     
     oldPos = ylab.Position;
     ylab.Position = [oldPos(1)+0.1 oldPos(2)];
     
     yTop = 0.0;
    yBottom = min(yzoom)-max(yzoom);
    
     axPorosity.YTick = [yBottom yTop];
     
    yminStr = sprintf(formatStr, yBottom);
    ymaxStr = sprintf(formatStr, yTop);
    axPorosity.YTickLabels = {yminStr, ymaxStr};
    
set(axPorosity, 'TickLength',[0 0])


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





else
     axPorosity.YTick = [];
end


if sum(drawLabels) == 0
  %  axPorosity.XTick = [];
  %  axPorosity.YTick = [];

    set(axPorosity, 'TickLength',[0 0])

end





set(axPorosity, 'position', axisExtent);

axSl = axes;

if drawSl
    contourValsSl = linspace(-0.98, 0, 7);

    %colormap(axSl, makeColorMap([0 0 0] , [230, 175, 0]/255, [200 0 0]/255));
    salinityColormap = makeColorMap([0 0 0] ,  [200 0 0]/255);
    salinityColormap = makeColorMap([1 1 1] ,  [0.5 0.5 0.5], [200 0 0]/255);
    salinityColormap = makeColorMap([200 0 0]/255, [200 0 0]/255);
    colormap(axSl, salinityColormap);
    caxis(axSl, [-1 0]);
    [CSl, hSl] = contour(X, Y, Sl.', contourValsSl);
    hSl.LineWidth = 2.0;

end

set(axSl,'dataAspectRatio',[1 1 1])
set(axSl, 'position', axisExtent);
set(axSl,'TickLength',[0 0])

axis(axSl, 'off');
axSl.Visible = 'off';

if drawPsi 
axPsi = axes;
if maxPsi > 0
    %maxPsi = 0.2; %max(max(abs(streamfunction)));
    psiVals = linspace(-maxPsi,maxPsi,2*Npsi);
    psiColormap = makeColorMap([1 1 1], [1 1 1]);
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

end


%Schlieren
axSchlieren = axes;
colormap(axSchlieren, parula);
pcolor(X,Y,schlieren.');
caxis([min(min(schlieren)) maxSlSchlieren]);



set(axSchlieren,'dataAspectRatio',[1 1 1])
set(axSchlieren, 'position', axisExtent);
set(axSchlieren,'TickLength',[0 0])

axis(axSchlieren, 'off');
axSchlieren.Visible = 'off';

% Text annotation
if drawTime
    textPos = [axisExtent(1)+0.0045 axisExtent(2)+axisExtent(4)-0.0905 0.1 0.09];
    an = annotation('textbox', textPos, 'String',sprintf('t=%1.3f',t),'LineStyle','none');
else
    fprintf('t=%1.3f\n',t);    
end



linkaxes([axPorosity axSl]);





end