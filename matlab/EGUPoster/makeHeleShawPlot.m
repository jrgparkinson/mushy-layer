
function [axPorosity, axSl] = makeHeleShawPlot(dataFile, depth, doLabel, axPos, SlMax)

blues = makeColorMap([8 48 107]/255, [104, 172, 212]/255, [1 1 1]); blues = flipud(blues);

%slMap = makeColorMap([  0.9763    0.9831    0.0538] , [ 0.0488    0.5772    0.8228], [8 48 107]/255); 
slMap = makeColorMap([  0.9763    0.9831    0.0538] , [ 0.6488    0.7772    0.2228], [8 48 107]/255); 

axPorosity = axes;

load(dataFile);

zmin = 900;
zmax=1024;

zmin = round(1024*(1-depth/8));


Xplot = X(zmin:zmax, :);
Zplot = Z(zmin:zmax, :);
chiPlot = porosity(:, zmin:zmax).';
SlPlot = Sl(:, zmin:zmax).';

pcolor(Xplot,Zplot,chiPlot);

colormap(axPorosity, (blues));

axPorosity.XTick = [];
if doLabel
    axPorosity.YTick = [min(min(Z(zmin:zmax, :))) max(max(Z)) ];
    bottomLabel = sprintf('%1.1f', -depth);
    axPorosity.YTickLabels = { bottomLabel, '0.0'};
    ylab = ylabel('$z$ (cm)');
    oldPos = ylab.Position;
    %ylab.Position = [oldPos(1)+0.08 oldPos(2)];
else
    axPorosity.YTick = [];
end

daspect([1 1 1]);

axPorosity.Position = axPos;

caxis([0 1]);

%Sl axes
axSl = axes;

smoothTransition = false;

%SlMax = 0.16;
SlDiff  =  min(SlPlot+(1-SlMax), 0);

if smoothTransition
SlDiff(chiPlot<0.8) = NaN;
else
    SlDiff(chiPlot<1.0) = NaN;
end

SlDiff = -min(min(SlDiff))+SlDiff;



colormap(axSl, flipud(slMap));
%minSl = -1;
avPorosity = mean(chiPlot.^6, 1);

chiFilter = repmat(avPorosity,256,1); %porosity.^2;

if smoothTransition
SlDiff =SlDiff.*chiFilter;
end

%min(min(SlDiff(:, zmin:zmax)))

pcolor(Xplot, Zplot, SlDiff);
%caxis([min(min(SlDiff)) max(max(SlDiff))]);
caxis([0 SlMax]);
%colorbar('Location', 'eastoutside');

daspect([1 1 1]);

axis(axSl, 'off');
axSl.Visible = 'off';

axSl.Position = axPos;  

linkaxes([axPorosity axSl])

dimensionalTime = t*800/60;
if dimensionalTime > 48 && dimensionalTime < 51
    dimensionalTime = 50;
end
timeStr = sprintf('$t=%1.0f$ mins', dimensionalTime);
text(0.08, 8-depth+0.15,timeStr, 'Color', [1 1 1], 'FontSize', 14);


end
