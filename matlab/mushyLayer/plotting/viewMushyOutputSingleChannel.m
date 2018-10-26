clear all;
close all;

set(0, 'defaultlinelinewidth',2);
set(0, 'defaultaxeslinewidth',2);
set(0, 'defaultpatchlinewidth',2);
set(0, 'defaultAxesFontSize',22);
% Plot upper and lower branches

% Length scale (metres)
L = 2;

baseDir = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/';

dataFolder = '/media/parkinsonjl/FREECOM HDD/';

RaCs = [250, 500, 1000];
frames = [28000, 31200, 28600];

%RaCs = [250]; frames = [28000];

for RaCi = 1:length(RaCs)
    RaC = RaCs(RaCi);
    output_dir = [ dataFolder, 'mushyLayer-periodic/CR1.25RaC',num2str(RaC),'.0Le200ChiCubedPermeabilitypts512-0/'];
    plot_prefix = ['mushyLayer-periodic-CR1.25RaC',num2str(RaC),'.0Le200ChiCubedPermeabilitypts512-'];
    
    frame = frames(RaCi);
    dim = 2; subcycled = true;
    
    output = MushyLayerOutput(dim, frame, output_dir, plot_prefix, subcycled);
    
    porosity = output.dataForComp(output.components.Porosity);
    T = output.dataForComp(output.components.Temperature);
    Sl = output.dataForComp(output.components.Liquidconcentration);
    streamfunction = output.getStreamfunction(3000, 1);
    
    xvel = output.dataForComp(output.components.xAdvectionvelocity);
    yvel = output.dataForComp(output.components.yAdvectionvelocity);
    
    zPorosity = output.horizontallyAverage(output.components.Porosity, false);
    V_averaged = output.horizontallyAverage(output.components.yAdvectionvelocity, true);
    
    probDomain = output.problemDomain;
    dx = probDomain.dxCoarse;
    
    numx = (probDomain.domainExtent.hi_i - probDomain.domainExtent.lo_i) + 1;
    numy = (probDomain.domainExtent.hi_j - probDomain.domainExtent.lo_j) + 1;
    
    x = linspace(double(probDomain.domainExtent.lo_i), double(probDomain.domainExtent.hi_i), numx)*dx;
    y = linspace(double(probDomain.domainExtent.lo_j), double(probDomain.domainExtent.hi_j), numy)*dx;
    
    xlo = double(probDomain.domainExtent.lo_i)*dx;
    xhi = double(probDomain.domainExtent.hi_i)*dx;
    width = xhi-xlo;
    
    zoom_lo_i = 5+64; % There's a bit of an offset to the first channel
    zoom_hi_i = zoom_lo_i+66;
    zoom_lo_j = 40;
    zoom_hi_j = 96;
    
    figureWidth = 600;
    figureHeight = round(figureWidth*(zoom_hi_j - zoom_lo_j)/64);
    
    %zoom_lo_z = zoom_lo_j*dx;
    %zoom_hi_z = zoom_hi_j*dx;
    streamfunction = streamfunction.';
    
    chi = porosity(zoom_lo_i:zoom_hi_i, zoom_lo_j:zoom_hi_j);
    U = xvel(zoom_lo_i:zoom_hi_i, zoom_lo_j:zoom_hi_j).';
    V = yvel(zoom_lo_i:zoom_hi_i, zoom_lo_j:zoom_hi_j).';
    psi = streamfunction(zoom_lo_i:zoom_hi_i, zoom_lo_j:zoom_hi_j);
    xzoom = x(zoom_lo_i:zoom_hi_i);
    yzoom = y(zoom_lo_j: zoom_hi_j);
    Tzoom = T(zoom_lo_i:zoom_hi_i, zoom_lo_j:zoom_hi_j).';
    Slzoom = Sl(zoom_lo_i:zoom_hi_i, zoom_lo_j:zoom_hi_j).';
    
    zPorosityzoom = zPorosity(zoom_lo_j: zoom_hi_j);
    V_averagedzoom = V_averaged(zoom_lo_j: zoom_hi_j);
    [X, Y] = meshgrid(xzoom, yzoom);
    
    
    
    %Work out edge of mush-liquid boundary
    chi;
    
    %[mush_x, mush_y] = find(chi<1.0);
    %idx_liquid = 1-idx_mush;
    refinement = 4*round(RaC/250);
    refinement = 16;
    chiRef = resizem(chi, refinement, 'bicubic');
    Xref = resizem(X, refinement, 'bilinear');
    Yref = resizem(Y, refinement, 'bilinear');
    chiRef = smoothn(chiRef, 'robust');
    
    offset = 1.0; numContours = 8; maxVal = 80;
    contourValsNeg = linspace(-maxVal, -offset, numContours);
    contourValsPos = linspace(offset, maxVal, numContours);
    
    
    TcontourVals = -0.5:0.1:1.5;
    SlcontourVals = -0.95:0.1:0;  
    
    h = figure();
    set(h, 'Position', [200 200 figureWidth figureHeight]);
    %colormap bone;
    % This should be the 'blues' colormap from visit
    colormap(makeColorMap([8 48 107]/255, [104, 172, 212]/255, [1 1 1]));
    box on;
    hold on;
    imagesc(xzoom, yzoom, chi.');
    
    %hStream=streamslice(X, Y, U, V, 0.1, 'method', 'cubic');
    iStart = round(length(X(1, :))/2)+1; % Channel isn't exactly central

    psi = psi.';
    
    %[CconNeg, hconNeg] = contour(X(:, 1:iStart), Y(:, 1:iStart), psi(:, 1:iStart), contourValsNeg);
    
    iStart = iStart - 4;
    %[CconPos, hconPos] = contour(X(:, 1:iStart), Y(:, 1:iStart), psi(:, 1:iStart), contourValsPos);
    
    iStart = iStart + 4 ;
    
    [CconNeg, hconNeg] = contour(X, Y, psi, contourValsNeg);
    [CconPos, hconPos] = contour(X, Y, psi, contourValsPos);
    
    
    %[Tcon, Th] = contour(X(:, iStart:end), Y(:, iStart:end), Tzoom(:, iStart:end), TcontourVals);
    %[Slcon, Slh] = contour(X(:, iStart:end), Y(:, iStart:end), Slzoom(:, iStart:end), SlcontourVals);
    
    %Make arrows!
    plotArrows(CconNeg, hconNeg, 'k', 0.8);
    plotArrows(CconPos, hconPos, 'k', 0.8);
    
    
    
    %h = quiver(hconPos.XData(ny, nx), hconPos.YData(ny, nx), U(ny, nx), V(ny, nx));
    
    [cml, hml] = contour(Xref, Yref, chiRef.', [-0.6 0.99]);
    %[cml2, hml2] = contour(X, Y, chi.', [-0.1 0.99]);
    
    hold off;
    
    axis([xzoom(1)+dx xzoom(end)-dx yzoom(1)+dx yzoom(end)-dx]);
    
    ax = gca;
       
    
    set(gca,'linewidth',2); %make box border bold to hide chimney outline
    set(ax,'clim',[0 1]);
    set(gca, 'Layer', 'top'); % put box border above imagesc
    
    hconNeg.LineColor = 'k';
    hconNeg.LineWidth = 2;
    hconPos.LineColor = 'k';
    hconPos.LineWidth = 2;
   
    
    hml.LineWidth = 2;
    hml.LineColor = 'r';
    hml.LineStyle = '-';
    
%     Th.LineColor = 'k';
%     Th.LineStyle = '--';
%     Th.LineWidth = 2;
     
%     Slh.LineColor = 'k';
%     Slh.LineStyle = '--';
%     Slh.LineWidth = 2;
    
    % hml2.LineWidth = 1;
    % hml2.LineColor = 'g';
    % hml2.LineStyle = ':';
    
    %ax.YTick =
    n = -2;
    
    %disp_lo_z = round(zoom_lo_z*10^(-n))/(10^(-n));
    %disp_hi_z = round(zoom_hi_z*10^(-n))/(10^(-n));
    %disp_hi_x = round(zoom_hi_x*10^(-n))/(10^(-n));
    
    %ax.YTickLabel ={num2str(disp_lo_z),num2str(disp_hi_z)};
    ax.XTick = [];
    ax.YTick = [];
    %ax.XTickLabel ={num2str(disp_hi_x), '0'};
    
    set(ax,'dataAspectRatio',[1 1 1])
    set(ax,'Ydir','normal');
    
    set(gca, 'position', [0.01 0.01 0.98 0.98]);
    %set(gca, 'position', [-0.05 0.5 0.9 0.73]);
    %ylabel('$z$', 'Interpreter','latex');
    %xlabel('$x$', 'Interpreter','latex');
    
    %set(ax,'TickLength',[0 0])
    %set(get(gca,'YLabel'),'Rotation',0, ...
    %   'Units', 'Normalized', 'Position', [-0.8, 0.5, 0])
    
    % c = colorbar();
    % c.Label.String = '\Pi';
    % c.Label.Rotation = 0;
    % c.Label.Position = [2.0 0.6];
    % c.Position = [0.8 0.22 0.05 0.7];
    % c.Ticks = [0.01 0.99];
    % c.TickLabels = {'0', '1'};
    % c.TickLength = 0;
    
    set(gca,'TickLength',[0 0])
    
    h =gcf;
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(h,['RaC', num2str(RaC), 'singleChannel.pdf'],'-dpdf','-r0')
    
    %
    %     figure();
    %     hold on;
    %     yyaxis left;
    %     plot(yzoom, zPorosityzoom);
    %     ylabel('porosity');
    %
    %     yyaxis right;
    %     plot(yzoom, V_averagedzoom);
    %     ylabel('w');
    %     hold off;
end
