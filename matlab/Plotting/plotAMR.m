close all; 

% Load plotfile
output_dir = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/output/';
plot_prefix = 'fixedChill-periodic-CR20.000RaC10Le100KozenyPermeabilityDa1.0e-02R1.0e-03pts256-St5-domWidth32.0-';
frame = 864 ; %001024;
dim = 2;
subcycled=true;
amrFile = MushyLayerOutput(dim, frame, output_dir, plot_prefix, subcycled);
amrFile_smooth = MushyLayerOutput(dim, frame, output_dir, plot_prefix, subcycled, 'bicubic');

% Get data to plot
[X,Y] = amrFile.grid();
porosity = amrFile.dataForComp(amrFile.components.Porosity).';
Sl = amrFile.dataForComp(amrFile.components.Liquidconcentration).';
Sl_smooth = amrFile_smooth.dataForComp(amrFile.components.Liquidconcentration).';

[Xbase,Ybase] = amrFile.levelArray(1).grid();
 % Think there's a bug, need to shift this:
 coarseDx = amrFile.levelArray(1).dx;
%Xbase= Xbase-amrFile.levelArray(1).dx;
Sl_base = amrFile.levelArray(1).dataForComp(amrFile.components.Liquidconcentration);

xl = [X(1,1) X(end,end)];
yl = [Y(1,1) Y(end,end)];

xlbase = [Xbase(1,1) Xbase(end,end)];
ylbase = [Ybase(1,1) Ybase(end,end)];

xl=xlbase;
yl=ylbase;

% Start plotting stuff
h = figure();
h.Position = [200 200 7;00 800];
hold on;

axPorosity = gca;

% Porosity first
pcolor(X,Y,porosity);
cPorosity = colorbar();

daspect([ 1 1 1]);
xlim(xl);
ylim(yl);
box on;

% Now liquid salinity contours
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

cSl = colorbar();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Level outlines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%domainBox = polyshape([X(1,1) X(1,1) X(end,end) X(end,end)], [Y(1,1) Y(end,end) Y(end,end) Y(1,1)]);
domainBox = polyshape([xl(1) xl(1) xl(end) xl(end)], [yl(1) yl(end) yl(end) yl(1)]);

meshes = amrFile.getMeshes(domainBox);

% Plot meshes 
mesh2 = plot(axPorosity, meshes(2)); % meshes on chombo level 1
mesh3 = plot(axPorosity, meshes(3)); % meshes on chombo level 2
%plot(domainBox);

mesh2.FaceAlpha = 0;
mesh3.FaceAlpha = 0;

mesh2.EdgeColor = [1 0 0];
mesh3.EdgeColor = [1 0 1];

mesh2.LineWidth = 2;
mesh3.LineWidth = 2;


% Stop plotting stuff now
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Link axes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
linkaxes([axPorosity, axSl]);

set(axPorosity, 'Layer', 'top');

axPos = [0.07 0.15 0.6 0.75];
axPorosity.Position = axPos;
axSl.Position = axPos;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Colorbars
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cPorosity.Label.String = '\chi';

cSl.Label.String = 'S_l';
cSl.Position = [0.85 cSl.Position(2) cSl.Position(3) cSl.Position(4)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Axis labels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

box on;


axLabels = axPorosity;

%daspect([ 1 1 1]);

%xlim(xl);
%ylim(yl);

xlabel(axPorosity, '$x$');
ylabel(axPorosity, '$y$');

axLabels.XTick = xl;
axLabels.YTick = yl;

format = '%1.1f';
axLabels.XTickLabels = {sprintf(format, xl(1)), sprintf(format, xl(2))};
axLabels.YTickLabels = {sprintf(format, yl(1)), sprintf(format, yl(2))};
