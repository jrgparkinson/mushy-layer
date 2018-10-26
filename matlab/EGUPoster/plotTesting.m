close all;
load('/home/parkinsonjl/convection-in-sea-ice/MushyLayer/matlab/toPlot23600.mat')
blues = makeColorMap([8 48 107]/255, [104, 172, 212]/255, [1 1 1]); blues = flipud(blues);

axPos = [0.1 0.2 0.7 0.7];

slMap = makeColorMap([  0.9763    0.9831    0.0538] , [ 0.0488    0.5772    0.8228], [8 48 107]/255); 
h = figure();
h.Position = [200 200 600 400];
z = Z(:, 1); zmin = 800; zmax = 1024;

axPorosity = axes;
colormap(axPorosity, blues);
pcolor(X(zmin:zmax, :),Z(zmin:zmax, :),porosity(:, zmin:zmax).');

%daspect([1 1 1]);

axPorosity.Position = axPos;

axSl = axes;

smoothTransition = false;

SlDiff  =  min(Sl+0.9, 0);

if smoothTransition
SlDiff(porosity<0.8) = NaN;
else
    SlDiff(porosity<0.99) = NaN;
end

SlDiff = -min(min(SlDiff))+SlDiff;




%Sl(porosity<0.8) = NaN;
colormap(axSl, flipud(slMap));
%minSl = -1;
avPorosity = mean(porosity.^6, 1);

chiFilter = repmat(avPorosity,256,1); %porosity.^2;

if smoothTransition
SlDiff =SlDiff.*chiFilter;
end

%min(min(SlDiff(:, zmin:zmax)))

pcolor(X(zmin:zmax, :),Z(zmin:zmax, :),SlDiff(:, zmin:zmax).');
caxis([min(min(SlDiff(:, zmin:zmax))) max(max(SlDiff))]);
colorbar('Location', 'eastoutside');

%daspect([1 1 1]);

axis(axSl, 'off');
axSl.Visible = 'off';

axSl.Position = axPos;

linkaxes([axPorosity axSl])