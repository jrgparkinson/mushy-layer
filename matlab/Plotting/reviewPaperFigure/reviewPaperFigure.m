%clear all;
function reviewPaperFigure

% Boring setup stuff
close all;
set(0, 'defaultlinelinewidth',2);
set(0, 'defaultaxeslinewidth',2);
set(0, 'defaultpatchlinewidth',2);
fontsize = 18;
set(0, 'defaultAxesFontSize',fontsize);
SlMax = 0.085; %55g/kg


% Where the data files are stored
dataFolder = '/home/parkinsonjl/mushy-layer/matlab/EGUPoster/';
data_dir = [dataFolder, 'T-15/'];

% Main options to change

%savePlots = false;
savePNG = false;

saveEPS = false;
savePDF = false;

%Possible frames: 2000, 13000, 19000, 31500
resStr = 'lowRes';
frameLeft = 2000;     frameRight = 13000;
depth = 1.2;
duplicate = false;

resStr = 'hiRes';
data_dir = '/home/parkinsonjl/mushy-layer/matlab/Plotting/reviewPaperFigure/';
frameLeft = 9000;     
frameRight = 48400; % Ideally make this 100000, once processed
frameRight = 76100; 
depth = 4.8;
duplicate = true;


% Style options: 
% 1 - just porosity, 
% 2 - porosity/sl mix, 
% 3 - porosity + sl contours
styles = {'Porosity', 'PorositySl', ...
    'PorosityEqualContours', 'PorosityUnequalContours', ...
    'PorositySlUnequalTContours'};
styleLeft = 5;          styleRight = 5;

%depth = 1.2;
doLabel = false;

% Now do the plot. A few options.
% Either a) plot the 'fake hele shaw'
% Or b) plot just the porosity.
saveName = ['reviewPaper-', resStr,  num2str(frameLeft), styles{styleLeft}, '_', num2str(frameRight), styles{styleRight}];

fprintf('Save name: %s \n', saveName);

h = figure();
plotHeight = 330;
plotWidth = 800;
set(h, 'Position', [200 300 plotWidth plotHeight]);

width = 0.46;
height = 0.8;
ax_bottom = 0.17;

cbar_bottom = 0.13;
cbar_height = 0.06;

m = 1; n=2;

subplot(m, n, 1);
axPosLeft = [0.03 ax_bottom width height];
dataFile = [data_dir, 'toPlot', num2str(frameLeft),'.mat'];
[axesLeft, cbarPorosity] = makePlot(dataFile, depth, doLabel, axPosLeft, SlMax, styleLeft, 'Porosity', duplicate);
%title('(a)');
cbarPorosity.Position = [0.03 cbar_bottom width cbar_height];
cbarPorosity.Ticks = [0 1];
cbarPorosity.Label.String = '\chi';

cbarPorosity.Label.Position = [0.5 0];

subplot(m, n, 2);
axPosRight = [0.52 ax_bottom width height];
dataFile = [data_dir, 'toPlot', num2str(frameRight),'.mat'];
[axesRight, cbarSl] = makePlot(dataFile, depth, doLabel, axPosRight, SlMax, styleRight, 'Sl', duplicate);
%title('(b)');
cbarSl.Position = [0.52 cbar_bottom width cbar_height];
cbarSl.Label.String = 'S_l';
cbarSl.Ticks = [0 0.08];
cbarSl.Label.Position = [0.04 0];



%title(['$T_a=-15^\circ$ C'], 'FontSize', fontsize);

%xlabel('$x$');
%ylabel('$z$');


    
    h = gcf;
    h.InvertHardcopy = 'off';
    h.Color = 'white';
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    filename = ['/home/parkinsonjl/convection-in-sea-ice/MushyLayer/matlab/Plotting/', saveName];
    
    if savePDF
        print(h,[filename, '.pdf'],'-dpdf','-r0')
    end
    
    if savePNG
        print(h,[filename, '.png'],'-dpng','-r400')
    end
    
    if saveEPS
        print(h,[filename, '.eps'],'-depsc','-r0')
    end
    


end

function [plotAxes, cbar] = makePlot(dataFile, depth, doLabel, axPos, SlMax, style, cbarField, duplicate)

plotAxes = [];
cbar = '';

%blues = makeColorMap([8 48 107]/255, [104, 172, 212]/255, [1 1 1]); blues = flipud(blues);
slMap = makeColorMap([  0.9763    0.9831    0.0538] , [ 0.6488    0.7772    0.2228], [8 48 107]/255); 

axPorosity = gca;

load(dataFile);

zmin = 900;
zmax = size(X, 1);

if duplicate
    zmax = round(zmax*0.95);
end

zmin = round(1024*(1-depth/8));

Xplot = X(zmin:zmax, :);
Zplot = Z(zmin:zmax, :);
chiPlot = porosity(:, zmin:zmax).';
SlPlot = Sl(:, zmin:zmax).';
TPlot = T(:, zmin:zmax).';

if duplicate
    dx = Xplot(end, end) - Xplot(end, end-1);
    Xplot = [Xplot; Xplot(2:end, :)+Xplot(end, end)-dx];
    %dx = Xplot(end, end) - Xplot(end-1, end);
    %for j=1:size(Xplot, 1)
    %   Xplot(end+1, :) = Xplot(end, :) + dx; 
    %end
    
    %Xplot = [Xplot; Xplot];
    
    Zplot = [Zplot; Zplot(2:end, :)];
    
   % Zplot = [Zplot; Zplot+Zplot(end, end)];
    
    %chiPlot = [chiPlot; chiPlot(2:end, :)];
    %SlPlot = [SlPlot; SlPlot(2:end, :)];
    %TPlot = [TPlot; TPlot(2:end, :)];
    
    chiPlot = repmat(chiPlot, 1, 2);
    SlPlot = repmat(SlPlot, 1, 2);
    TPlot = repmat(TPlot, 1, 2);
    
    x = linspace(0, size(chiPlot, 2)*dx, size(chiPlot, 2));
    z = linspace(0,  size(chiPlot, 1)*dx, size(chiPlot, 1));
    [Xplot, Zplot] = meshgrid(x,z);
    
end

pcolor(Xplot, Zplot, chiPlot);
 
if style == 2 || style == 5
   colormap(axPorosity, blues);
else
   colormap(axPorosity, viridis); 
end

axPorosity.XTick = [];
axPorosity.YTick = [];

daspect([1 1 1]);

if strcmp(cbarField, 'Porosity')
   cbar = colorbar(axPorosity, 'Location', 'southoutside'); 
end

axPorosity.Position = axPos;

caxis(axPorosity, [0 1]);

plotAxes(end+1) = axPorosity;



%Sl axes
if style > 1
    
   

    if style == 2 || style == 5
        axSlColor = axes;
        
        smoothTransition = false;
        SlDiff  =  min(SlPlot+(1-SlMax), 0);

        if smoothTransition
            SlDiff(chiPlot<0.8) = NaN;
        else
            SlDiff(chiPlot<1.0) = NaN;
        end

        SlDiff = -min(min(SlDiff))+SlDiff;

        colormap(axSlColor, flipud(slMap));
        avPorosity = mean(chiPlot.^6, 1);
        chiFilter = repmat(avPorosity,256,1); %porosity.^2;

        if smoothTransition
            SlDiff = SlDiff.*chiFilter;
        end

        pcolor(Xplot, Zplot, SlDiff);

        caxis([0 SlMax]);
        
        daspect([1 1 1]);

    axis(axSlColor, 'off');
    axSlColor.Visible = 'off';
    
    if strcmp(cbarField, 'Sl')
        cbar = colorbar(axSlColor, 'Location', 'southoutside'); 
    end

    axSlColor.Position = axPos;  

    linkaxes([axPorosity axSlColor])

    caxis(axPorosity, [0 1]);
    
    plotAxes(end+1) = axSlColor;
    
    end
    
    if style == 3 || style == 4 || style == 5
        
        axSl = axes;
        
        if style == 3 || style == 4
            contourPlotField = SlPlot;
        elseif style == 5
            contourPlotField = TPlot; % TODO: should be temperature
        end
      
        minSl = min(min(contourPlotField));
        maxSl = max(max(contourPlotField));
        
        contourValsFine = linspace(minSl*1.05, maxSl, 40); %num porosity contours
        contourValsCoarse = linspace(minSl*1.05, maxSl, 10); %num porosity contours
        
        equalSpacing = linspace(0, 1, 10);
        unevenSpacing = equalSpacing.^3; %(exp(equalSpacing).^2-1)/(exp(1)^2-1);
        
        if style == 3 || style==5
           spacing = equalSpacing;
        else
           spacing = unevenSpacing;
        end
        
        if style == 5
            linearSpace = linspace(0, 1, 10);
            scaling = 10^2.5;
            contourVals = log10(1+scaling*linearSpace)/log10(1+scaling); 
        else
            contourVals = minSl+0.005 + (maxSl-minSl)*spacing;
        end
        % Find contour vals concentrated around plume
        %fineFrac = 1/4;
        % contourVals = [contourValsFine(1:round(end*fineFrac)) ...
        %     contourValsCoarse(round(end*fineFrac):end)];
        
        cMap = makeColorMap([0 0 0], [0 0 0]);
        cMap = makeColorMap([1 0 0], [1 0 0]);
        
        colormap(axSl, cMap);
        %   caxis(axSl, [minSl maxSl]);
        [Cchi, hchi] = contour(Xplot, Zplot, contourPlotField, contourVals);
        
        caxis(axSl, [minSl maxSl]);
        
        hchi.LineWidth = 1.0;
        
        daspect([1 1 1]);

        axis(axSl, 'off');
        axSl.Visible = 'off';
        axSl.Position = axPos;  

        linkaxes([axPorosity axSl])

        caxis(axPorosity, [0 1]);
        
        plotAxes(end+1) = axSl;
              
    end

    
end



dimensionalTime = t*800/60;
if dimensionalTime > 48 && dimensionalTime < 51
    dimensionalTime = 50;
end
%timeStr = sprintf('$t=%1.0f$ mins', dimensionalTime);
%text(0.08, 8-depth+0.15,timeStr, 'Color', [1 1 1], 'FontSize', 14);



end

