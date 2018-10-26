clear all;
close all;

% Plot upper and lower branches

% Length scale (metres)
L = 2;
RaC = 400;
ConcRatio = 1.25;

baseDir = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/';

dataFolder = '/media/parkinsonjl/FREECOM HDD/';
%upperBranch = 'mushyLayerLowC-upperBranch-insulating';
upperBranch = 'optimalStates-highRes';

widths = 20:1:50;

fluxes = NaN*widths;
fluxesLower = NaN*widths;

if exist('poutUpper', 'var') == 0
    poutUpper = Pout();
end

if exist('poutLower', 'var') == 0
    poutLower = Pout();
end

% Variable width by keeping constant aspect ratio
if exist('poutAlternative', 'var') == 0
    poutAlternative = Pout();
end

for i=1:length(widths)
    if (length(poutUpper) < i || ~poutUpper(i).steadyState)
        fname = [dataFolder , upperBranch, '/CR',num2str(ConcRatio), 'RaC',num2str(RaC), ...
        'Le200ChiCubedPermeabilitypts', ...
        num2str(widths(i)),'-0/pout.0'];
    poutUpper(i) = Pout(fname);
    end
    
    if (length(poutLower) <i || ~poutLower(i).steadyState)
    poutLower(i) = Pout([dataFolder , 'mushyLayerLowC-lowerBranch-insulating/CR',num2str(ConcRatio), 'RaC',num2str(RaC), ...
        'Le200ChiCubedPermeabilitypts', ...
        num2str(widths(i)),'-0/pout.0']);
    end
end

alternativeWidths = [0.047 0.0475 0.048 0.0485 0.049 0.0495];
alternativeWidths = [];
alternativeFluxes = alternativeWidths*NaN;
for i=1:length(alternativeWidths)
    poutAlternative(i) = Pout([dataFolder , upperBranch, '/CR',num2str(ConcRatio), 'RaC',num2str(RaC), ...
        'Le200ChiCubedPermeabilitypts32width', ...
        num2str(alternativeWidths(i)),'-0/pout.0']);
    
    if length(poutAlternative(i).times) ~= 0
       alternativeFluxes(i) = poutAlternative(i).fluxBottom(end);
    end
end

for i=1:length(widths)
    %/media/parkinsonjl/FREECOM HDD/mushyLayerLowC-periodic/CR1.25RaC250Le200ChiCubedPermeabilitypts96-0
    poutNames{i} = num2str(0.2*widths(i)/128);
    if length(poutUpper(i).times) ~= 0

        fluxes(i) = poutUpper(i).fluxBottom(end);
        
    end
    
    if length(poutLower(i).times) ~= 0
       fluxesLower(i) = poutLower(i).fluxBottom(end);
    end
    
    
end


%Fill lower branch where we haven't done simulations
for i=(length(widths)-1):-1:1
   if isnan(fluxesLower(i)) &&  fluxesLower(i+1) < 1.0001
       fluxesLower(i) = 1.0;
   end
end




% poutNames = {'Ra200', 'Ra250', 'Ra300'};
lineStyles = {'-', '--',':', '-.'};
lineColours = {'r', 'b', 'g', 'm', 'black'};
lineColours = {'m', 'c', 'r', 'g', 'b', 'k'};
h = figure();
set(h, 'Position', [200 200 600 400]);


legendStr = {};


actualWidths = widths*(0.2/128);
%actualWidths = actualWidths*L*2; %convert to metres
%actualWidths =actualWidths*100; %(cm)
hold on;
plot(actualWidths(~isnan(fluxes)), fluxes(~isnan(fluxes))-1, 'x-');
plot(actualWidths(~isnan(fluxesLower)), fluxesLower(~isnan(fluxesLower))-1, 'o-');
plot(alternativeWidths, alternativeFluxes-1, 'd-');
hold off;
axis([0 0.2 0 (max(max(fluxes), max(alternativeFluxes))-1)*1.1]);

%legend({'Upper branch', 'Lower branch'});

%title('Steady state vertical solute flux');

xlabel('$\lambda = l/(2H)$', 'Interpreter','latex');
ylabel('$F_s - V$', 'Interpreter','latex');

set(get(gca,'YLabel'),'Rotation',0, ...
    'Units', 'Normalized', 'Position', [-0.18, 0.55, 0])

set(gca, 'position', [0.25 0.22 0.7 0.73]);

%grid on; 
box on;

% For saving figure as PDF
 h =gcf;
 set(h,'Units','Inches');
  pos = get(h,'Position');
 set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
 print(h,['soluteFluxRaC',num2str(RaC),'CR',num2str(ConcRatio), '.pdf'],'-dpdf','-r0')

