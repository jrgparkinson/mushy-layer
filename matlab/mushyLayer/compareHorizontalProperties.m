clear all;
close all;

set(0, 'defaultlinelinewidth',2);
set(0, 'defaultaxeslinewidth',2);
set(0, 'defaultpatchlinewidth',2);
set(0, 'defaultAxesFontSize',22);

% Length scale (metres)
L = 2;

axisLabels = true;

baseDir = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/';

dataFolder = '/media/parkinsonjl/FREECOM HDD/mushyLayerLowC-upperBranch-insulating/';

files = []; legendStr = {};
y_wcrit = [];
y_chicrit = [];
y_psicrit = [];
chicrit = 0.4; wcrit = 0.8; psicrit = 0.5;
RaCs = [];

i = 0;

i = i+1;
files(i).output_dir = 'CR1.25RaC75Le200ChiCubedPermeabilitypts48-0/';
files(i).plot_prefix = 'mushyLayerLowC-upperBranch-insulating-CR1.25RaC75Le200ChiCubedPermeabilitypts48-';
files(i).frame = 59592;
files(i).legend = 'Ra = 75';
files(i).Ra = 75;

%
i = i+1;
files(i).output_dir = 'CR1.25RaC100Le200ChiCubedPermeabilitypts40-0/';
files(i).plot_prefix = 'mushyLayerLowC-upperBranch-insulating-CR1.25RaC100Le200ChiCubedPermeabilitypts40-';
files(i).frame = 50771;
files(i).legend = 'Ra = 100';
files(i).Ra = 100;

% Optimal state
% i=i+1;
%  files(i).output_dir = 'CR1.25RaC150Le200ChiCubedPermeabilitypts40-0/';
%  files(i).plot_prefix = 'mushyLayerLowC-upperBranch-insulating-CR1.25RaC150Le200ChiCubedPermeabilitypts40-';
%  files(i).frame = 45631;
%  files(i).legend = 'Ra = 150';
%  files(i).Ra = 150;

% Optimal state
i = i+1;
files(i).output_dir = 'CR1.25RaC200Le200ChiCubedPermeabilitypts32-0/';
files(i).plot_prefix = 'mushyLayerLowC-upperBranch-insulating-CR1.25RaC200Le200ChiCubedPermeabilitypts32-';
files(i).frame = 43334;
files(i).legend = 'Ra = 200';
files(i).Ra = 200;

% Optimal state
% i = i+1;
% files(i).output_dir = 'CR1.25RaC300Le200ChiCubedPermeabilitypts28-0/';
% files(i).plot_prefix = 'mushyLayerLowC-upperBranch-insulating-CR1.25RaC300Le200ChiCubedPermeabilitypts28-';
% files(i).frame = 44705;
% files(i).legend = 'Ra = 300';
% files(i).Ra = 300;

% Optimal state
i = i+1;
files(i).output_dir = 'CR1.25RaC400Le200ChiCubedPermeabilitypts28-0/';
files(i).plot_prefix = 'mushyLayerLowC-upperBranch-insulating-CR1.25RaC400Le200ChiCubedPermeabilitypts28-';
files(i).frame = 141902;
files(i).legend = 'Ra = 400';
files(i).Ra = 400;

% files(2).output_dir = 'CR1.25RaC250Le200ChiCubedPermeabilitypts32-0/';
% files(2).plot_prefix = 'mushyLayerLowC-insulating-CR1.25RaC250Le200ChiCubedPermeabilitypts32-';
% files(2).frame = 37322;
% files(2).legend = 'Ra = 250';

% Optimal state
% i = i+1;
%   files(i).output_dir = 'CR1.25RaC500Le200ChiCubedPermeabilitypts28-0/';
%   files(i).plot_prefix = 'mushyLayerLowC-upperBranch-insulating-CR1.25RaC500Le200ChiCubedPermeabilitypts28-';
%   files(i).frame = 142094;
%   files(i).legend = 'Ra = 500';
%   files(i).Ra = 500;



% Optimal state
i = i+1;
files(i).output_dir = 'CR1.25RaC600Le200ChiCubedPermeabilitypts56-0/';
files(i).plot_prefix = 'mushyLayerLowC-upperBranch-insulating-CR1.25RaC600Le200ChiCubedPermeabilitypts56-';
files(i).frame = 082852;
files(i).legend = 'Ra = 600';
files(i).Ra = 600;


%Optimal state
%  i = i+1;
%    files(i).output_dir = 'CR1.25RaC700Le200ChiCubedPermeabilitypts48-0/';
%    files(i).plot_prefix = 'mushyLayerLowC-upperBranch-insulating-CR1.25RaC700Le200ChiCubedPermeabilitypts48-';
%    files(i).frame = 117900;
%    files(i).legend = 'Ra = 700';
%    files(i).Ra = 700;

%optimal state
% i = i+1;
% files(i).output_dir = 'CR1.25RaC800Le200ChiCubedPermeabilitypts48-0/';
% files(i).plot_prefix = 'mushyLayerLowC-upperBranch-insulating-CR1.25RaC800Le200ChiCubedPermeabilitypts48-';
% files(i).frame = 122778;
% files(i).legend = 'Ra = 800';
% files(i).Ra = 800;

i = i+1;
files(i).output_dir = 'CR1.25RaC800Le200ChiCubedPermeabilitypts48-0/';
files(i).plot_prefix = 'mushyLayerLowC-upperBranch-insulating-CR1.25RaC800Le200ChiCubedPermeabilitypts48-';
files(i).frame = 160877;
files(i).legend = 'Ra = 800';
files(i).Ra = 800;

% i = i+1;
% files(i).output_dir = 'CR1.25RaC800Le200ChiCubedPermeabilitypts40-0/';
% files(i).plot_prefix = 'mushyLayerLowC-upperBranch-insulating-CR1.25RaC800Le200ChiCubedPermeabilitypts40-';
% files(i).frame = 147915;
% files(i).legend = 'Ra = 800';
% files(i).Ra = 800;


% files(3).output_dir = 'CR1.25RaC600Le200ChiCubedPermeabilitypts32-0/';
% files(3).plot_prefix = 'mushyLayerLowC-upperBranch-insulating-CR1.25RaC600Le200ChiCubedPermeabilitypts32-';
% files(3).frame = 122778;
% files(3).legend = 'Ra = 600';


% files(3).output_dir = 'CR1.25RaC800Le200ChiCubedPermeabilitypts64-0/';
% files(3).plot_prefix = 'mushyLayerLowC-upperBranch-insulating-CR1.25RaC800Le200ChiCubedPermeabilitypts64-';
% files(3).frame = 081882;
% files(3).legend = 'Ra = 800';

%  files(1).output_dir = 'CR1.25RaC600Le200ChiCubedPermeabilitypts32-0/';
%  files(1).plot_prefix = 'mushyLayerLowC-upperBranch-insulating-CR1.25RaC600Le200ChiCubedPermeabilitypts32-';
%  files(1).frame = 124863;
%  files(1).legend = 'Ra600-32';



%  files(1).output_dir = 'CR1.25RaC600Le200ChiCubedPermeabilitypts48-0/';
%  files(1).plot_prefix = 'mushyLayerLowC-upperBranch-insulating-CR1.25RaC600Le200ChiCubedPermeabilitypts48-';
%  files(1).frame = 94637;
%  files(1).legend = 'Ra600-48';
%
%   files(2).output_dir = 'CR1.25RaC600Le200ChiCubedPermeabilitypts56-0/';
%  files(2).plot_prefix = 'mushyLayerLowC-upperBranch-insulating-CR1.25RaC600Le200ChiCubedPermeabilitypts56-';
%  files(2).frame = 082852;
%  files(2).legend = 'Ra600-56';
%
%
%   files(3).output_dir = 'CR1.25RaC600Le200ChiCubedPermeabilitypts64-0/';
%  files(3).plot_prefix = 'mushyLayerLowC-upperBranch-insulating-CR1.25RaC600Le200ChiCubedPermeabilitypts64-';
%  files(3).frame = 63666;
%  files(3).legend = 'Ra600-64';

h = figure();
set(h, 'Position', [200 200 800 400]);
set(h,'defaultAxesColorOrder',[[0 0 0 ]; [0 0 0]]);

plots = [];
logPsiMaxVec = NaN*ones(length(files), 128);

%colors = {'b', 'r', 'k', 'm', 'g', 'c', 'y', 'b', 'r', 'k'};
%colors =  linspecer(length(files)+5);
colors = [   0    0.4470    0.7410
    0.8500    0.3250    0.0980
    
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840
    0.1     0.1     0.1
    0.9290    0.6940    0.1250];
%axes('NextPlot','replacechildren', 'ColorOrder',C);

logPsiMax = NaN*ones(length(files), 100);
yLogPsiMax = [];

box on;
hold on;

for file_i = 1:length(files)
    RaCs(file_i) = files(file_i).Ra;
    
    dim = 2; subcycled = true;
    
    output = MushyLayerOutput(dim, files(file_i).frame, [dataFolder, files(file_i).output_dir],...
        files(file_i).plot_prefix, subcycled);
    
    
    
    zPorosity = output.horizontallyAverage(output.components.Porosity, false);
    V_averaged = output.horizontallyAverage(output.components.yAdvectionvelocity, true);
    
    Vfield = output.dataForComp(output.components.yAdvectionvelocity);
    V_max = output.horizontallyMax(Vfield, true);
    
    psi = output.getStreamfunction(200, 1).';
    psi_max = output.horizontallyMax(psi, true);
    
    probDomain = output.problemDomain;
    dx = probDomain.dxCoarse;
    
    numx = (probDomain.domainExtent.hi_i - probDomain.domainExtent.lo_i) + 1;
    numy = (probDomain.domainExtent.hi_j - probDomain.domainExtent.lo_j) + 1;
    
    x = linspace(double(probDomain.domainExtent.lo_i), double(probDomain.domainExtent.hi_i), numx)*dx;
    y = linspace(double(probDomain.domainExtent.lo_j), double(probDomain.domainExtent.hi_j), numy)*dx;
    
    
    xlo = double(probDomain.domainExtent.lo_i)*dx;
    xhi = double(probDomain.domainExtent.hi_i)*dx;
    width = xhi-xlo;
    
    %[X, Y] = meshgrid(x, y);
    
    y_wcrit(file_i) = NaN;
    y_chicrit(file_i)  = NaN;
    y_psicrit(file_i)= NaN;
    for j = 1:numy
        
        if V_max(j) < wcrit && isnan(y_wcrit(file_i))
            y_wcrit(file_i) = y(j);
        end
        
        if zPorosity(j) < chicrit && isnan(y_chicrit(file_i))
            y_chicrit(file_i) =  y(j);
        end
        
        if psi_max(j) < psicrit && isnan(y_psicrit(file_i)) && y(j) > 0.05
            y_psicrit(file_i) = y(j);
        end
    end
    
    
    %color = colors{file_i};
    color = squeeze(colors(file_i, :));
    % Plot
    yyaxis right;
    %plots(end+1) = plot(y, zPorosity, '-', 'color',colors(file_i,:));
    %plots(end+1) = plot(y, exp(psi_max)/exp(max(psi_max)), '-', 'color',colors(file_i,:));
    plots(end+1) = plot(y, psi_max, '-', 'color',colors(file_i,:));
    
    yyaxis left;
    %plot(y, V_averaged*8, [color, ':']);
    
    %plot(y, log10(psi_max./max(psi_max)), [color, '--']);
    
    plot(y, log10(psi_max), '--', 'color',colors(file_i,:));
    
    logPsiMaxVec(file_i, :) = log10(psi_max);
    
    temp = log(psi_max);
    
    
    start =  round(length(temp)*0.35); finish =  round(length(temp)*0.6);
    num = finish- start;
    logPsiMax(file_i, 1:num+1) = temp(start:finish);
    yLogPsiMax = y(start:finish);
    
    legendStr{file_i} = files(file_i).legend;
    
end

hold off;

set(gca,'linewidth',2); %make box border bold to hide chimney outline

ax = gca;

lgd =  legend(plots, legendStr);
lgd.FontSize = 18;

if axisLabels
    yyaxis right;
    ylabel('$\bar{\chi}$', 'Interpreter','latex');
    set(get(gca,'YLabel'),'Rotation',0, ...
        'Units', 'Normalized', 'Position', [1.06, 0.78, 0])
    %ax.YTick=[0, 0.5, 1];
    %ax.YTickLabels={'0', '0.5', '1'};
    
    yyaxis left;
    %ylabel({'$\bar{|w|}$,','max$|w|$'}, 'Interpreter','latex');
    %set(get(gca,'YLabel'),'Rotation',0, ...
    %    'Units', 'Normalized', 'Position', [-0.16, 0.55, 0])
    ylabel('log(max$|\psi|)$', 'Interpreter','latex');
    set(get(gca,'YLabel'),'Rotation',90, ...
        'Units', 'Normalized', 'Position', [-0.12, 0.65, 0])
    ax.YTick=[0.5, 1, 1.5, 2];
    ax.YTickLabels={'0.5', '1', '1.5', '2', ''};
    
    %axis([0 0.2 -2 2]);
    y_lo = 0.5; y_hi =2;
    axis([0 0.2 y_lo y_hi]);
    
    xlabel('$y$', 'Interpreter','latex');
    set(gca, 'Position', [0.12 0.22 0.57 0.73]);
    lgd.Position = [0.78 0.22 0.2 0.73];
    
    patch_lo = [0.08, 0.5];
    patch_hi = [0.095, 1.0];
    
    p=patch([patch_lo(1) patch_lo(1) patch_hi(1) patch_hi(1)],...
        [patch_lo(2) patch_hi(2) patch_hi(2) patch_lo(2)],'k');
    set(p,'FaceAlpha',0.1);
    p.LineStyle = '-';
    p.LineWidth = 1.5;
    
    
    inset_lo = [0.48, 0.58];
    inset_size = [0.19 0.34];
    
   % Convert to actual units
    axis_inset_lo = [0.127, 1.245];
    axis_inset_size = [0.066 0.69];
    
    %Dotted lines to connect patches
    hold on;
    plot([patch_lo(1) axis_inset_lo(1)], ...
        [patch_hi(2) axis_inset_lo(2) + axis_inset_size(2)], ...
        ':k');
    plot([patch_hi(1) axis_inset_lo(1) + axis_inset_size(1)], ...
        [patch_lo(2) axis_inset_lo(2)], ...
        ':k');
    hold off; 
    % Inset axes
    axes('Position',[inset_lo(1) inset_lo(2)  inset_size(1)  inset_size(2)])
    box on;
    hold on;
    for i = 1:length(files)
        %  color = colors{i};
        plot(y,squeeze(logPsiMaxVec(i, :)), '--', 'color',colors(i,:));
    end
    hold off;
    axis([0.08 0.095 0.5 1.0]);
    ax = gca;
    ax.XTick = [0.0805 0.0885];
    ax.YTick = [0.52 0.98];
    ax.XTickLabels = {'0.08', '0.09'};
    ax.YTickLabels = {'0.5', '1.0'};
    
    p=patch([patch_lo(1) patch_lo(1) patch_hi(1) patch_hi(1)],...
        [patch_lo(2) patch_hi(2) patch_hi(2) patch_lo(2)],'k');
    set(p,'FaceAlpha',0.1);
    p.LineStyle = '-';
    p.LineWidth = 0.5;
    
    
    
    
    %     xlabel('y');
    %     ylabel('log(max$|\psi|)$', 'Interpreter','latex');
    %     set(get(gca,'YLabel'),'Rotation',90, ...
    %         'Units', 'Normalized', 'Position', [-0.12, 0.65, 0])
    
    
    set(ax,'fontsize',16)
    
    %h =gcf;
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(h,['verticalCrossSectionsWithLabelsLog.pdf'],'-dpdf','-r0')
    
else
    ax.XTick = [];
    yyaxis left; ax.YTick = [];
    yyaxis right;ax.YTick = [];
    set(ax,'TickLength',[0 0]);
    
    set(gca, 'Position', [0.03 0.03 0.92 0.92]);
    
    
    h =gcf;
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(h,['verticalCrossSections.pdf'],'-dpdf','-r0')
    
end



%figure();
% plot(RaCs, y_psicrit,RaCs, y_chicrit);
% legend(['psi =',num2str(psicrit)], ['chi = ', num2str(chicrit)]);
 

 figure(); plot(y(end-(length(logPsiMax)-1):end), logPsiMax);
 
 %Fit for each RaC:
 p = [];
 for i = 1:length(files)
     temp = logPsiMax(i, :);
     numNan = nnz(isnan(temp));
     thisLogPsi = temp(1:end-numNan);
    fit = polyfit(yLogPsiMax, thisLogPsi, 1)
    p(i, :) = fit;
     
 end
 
 
 figure(); 
 
 
%
% ylabel('y');
% xlabel('Ra');