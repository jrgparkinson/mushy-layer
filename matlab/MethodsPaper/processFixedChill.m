% This scripts makes use of the Streamfunction Color function:
% https://uk.mathworks.com/matlabcentral/fileexchange/69268-streamfunction-color

function processFixedChill(output_dir)

if nargin < 1
    output_dir = '/home/parkinsonjl/mnt/sharedStorage/Test/FixedChill-t100.0-Ra1e+07-Da5.0e-05-C5.00/Uniform-FixedChill-64--0/';
    output_dir = '/home/parkinsonjl/mushy-layer/execSubcycle/';
end

folders = dir(output_dir);

for i=1:length(folders)
    fprintf('Processing %s \n', folders(i).name);
    processSpecificFolder(fullfile(output_dir, folders(i).name));
end

end


function processSpecificFolder(output_dir)

close all; 

savePNG = true;

% Load plotfile
ml = getFinalPlotFile(output_dir);

if isempty(isprop(ml, 'levelArray'))
    fprintf('Folder does not contain any plot files \n');
    return;
end

% Get data to plot
[X,Y] = ml.grid();
porosity = ml.dataForComp(ml.components.Porosity).';
Sl = ml.dataForComp(ml.components.Liquidconcentration).';
U = ml.dataForComp(ml.components.xDarcyvelocity).';
V = ml.dataForComp(ml.components.yDarcyvelocity).';
Streamfunction = ml.dataForComp(ml.components.streamfunction).';

xl = [X(1,1) X(end,end)];
yl = [Y(1,1) Y(end,end)];

dx = X(1,2) - X(1,1);


% Start plotting stuff
h = figure();
h.Position = [200 200 900 700];
axPos = [0.08 0.08 0.8 0.75];
hold on;

%allAxes = [];

axPorosity = gca;

% Porosity first
pcolor(X,Y,porosity);
colormap(axPorosity, blues);
cPorosity = colorbar('Location', 'northoutside');

daspect([ 1 1 1]);
xlim(xl);
ylim(yl);
box on;


allAxes = [axPorosity];

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

    allAxes(end+1) = axQuiver;
end

plotStreamlines = true;
if plotStreamlines
    axStream = axes;
    
    % Modify data a little to make velocities more similar

    %streamslice(X,Y,U,V)
    pow = 1.0;
    Uplot = sign(U).*(abs(U).^pow);
    Vplot = sign(V).*(abs(V).^pow);
    
    vmag = sqrt(Uplot.^2+Vplot.^2);  
    dx = X(1, 2) - X(1,1);
    startx = xl(1):3.8*dx:xl(end);
    %starty = Y(1,1)*startx;
    
    [sy,sx] = meshgrid(Y(1,1),startx);
    
    %colormap(axStream, spring);
    %colormap(axStream, autumn);
    %cmap = colormap(axStream);
    %Streamcolor(X,Y,Uplot,Vplot,sx,sy, vmag);
    %colormap(axStream, hot);
    %minPsi = min(min(Streamfunction));
    maxPsi = max(max(abs(Streamfunction)));
    minPsi = -maxPsi; 
    
    dPsi = (maxPsi-minPsi)/10; % need even number of divisions
    v = minPsi+0.5*dPsi:dPsi:maxPsi-0.5*dPsi; % Avoid the zero contour
    
    colormap(axStream, autumn);
    cbar = streamfunctioncolor(X,Y,U,V,Streamfunction,v);
    cbar.Label.String = 'Scaled fluid velocity';
    cbar.Ticks = [0 1];

    axStream.Visible = 'off';
    
    
    daspect([ 1 1 1]);
    xlim(xl);
    ylim(yl);

    allAxes(end+1) = axStream;
    
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
    allAxes(end+1) = axSl;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Level outlines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(ml.levelArray) > 1
    
    domainBox = polyshape([xl(1) xl(1) xl(end) xl(end)], [yl(1) yl(end) yl(end) yl(1)]);

    meshes = amrFile.getMeshes(domainBox);

    % Plot meshes 
    pshape = shiftMesh(meshes(2), dx, xl, k, d);
    mesh2 = plot(allAxes(1), pshape); % meshes on chombo level 1

    %plot(domainBox);
    mesh2.FaceAlpha = 0;
    mesh2.EdgeColor = [1 0 1];
    mesh2.LineWidth = 2;

    if length(meshes) > 2
        pshape = shiftMesh(meshes(3),dx, xl,  k, d);
    mesh3 = plot(allAxes(1), pshape); % meshes on chombo level 2
    mesh3.FaceAlpha = 0;
    mesh3.EdgeColor = [0 1 1];
    mesh3.LineWidth = 2;
    end
end

% Stop plotting stuff now
hold off;


%axSl.Position = axPos;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Colorbars
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cPorosity.Label.String = 'Porosity';
cPorosity.Label.Position = [0.5 1.2];
cPorosity.Label.Rotation = 0;
cPorosity.Ticks = [0 1];
cPorosity.Limits = [0 1];
cPorosity.TickLabels = {'0 (solid)', '1 (liquid)'};

% cVel.Label.String = 'Fluid velocity';
% cVel.Ticks = [0 1];
% maxU = sqrt(max(max(abs(U.^2 + V.^2))));
% cVel.TickLabels = {'0', sprintf('%.0f', maxU)};
% cVel.Label.Position = [1.0 0.5];
%cPorosity.Label.Rotation = 0;
    %cSl.Label.String = 'S_l';
%cSl.Position = [0.85 cSl.Position(2) cSl.Position(3) cSl.Position(4)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Link axes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 2:length(allAxes)
linkaxes([allAxes(1), allAxes(i)]);
end

set(allAxes(1), 'Layer', 'top');


for i = 1:length(allAxes)
allAxes(i).Position = axPos;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Axis labels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

box on;


axLabels = allAxes(1);

%daspect([ 1 1 1]);

%xlim(xl);
%ylim(yl);

xlab = xlabel(allAxes(1), '$x$');
ylab = ylabel(allAxes(1), '$y$');

xlab.Position(2) = xlab.Position(2) + 0.02;
ylab.Position(1) = ylab.Position(1) + 0.02;

axLabels.XTick = xl;
axLabels.YTick = yl;


%format = '%1.1f';
axLabels.XTickLabels = {sprintf(format, xl(1)-dx/2), sprintf(format, xl(2)+dx/2)};
axLabels.YTickLabels = {sprintf(format, yl(1)-dx/2), sprintf(format, yl(2)+dx/2)};


% h = gcf;
h.InvertHardcopy = 'off';
h.Color = 'white';
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
filename = fullfile(output_dir, 'finalState.png');

%%if savePDF
%   print(h,[filename, '.pdf'],'-dpdf','-r0')
% end

if savePNG
    print(h,filename,'-dpng','-r400')
end
    
    
end


