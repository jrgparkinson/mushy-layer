clear all;
close all;

% Plot upper and lower branches

% Length scale (metres)
L = 2;

RaCs = [200, 400];

h = figure();
set(h, 'Position', [200 200 600 400]);


hold on;

for RaC_i = 1:length(RaCs)
    
RaC = RaCs(RaC_i);

baseDir = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/';

dataFolder = '/media/parkinsonjl/FREECOM HDD/';

widths = 0:1:160;

fluxes = NaN*widths;
fluxesLower = NaN*widths;

if exist('poutUpper', 'var') == 0
    poutUpper = Pout();
end

if exist('poutLower', 'var') == 0
poutLower = Pout();
end

for i=1:length(widths)
    if (length(poutUpper) < i || ~poutUpper(i).steadyState)
    poutUpper(i) = Pout([dataFolder , 'mushyLayerLowC-upperBranch-insulating/CR1.25RaC',num2str(RaC), ...
        'Le200ChiCubedPermeabilitypts', ...
        num2str(widths(i)),'-0/pout.0']);
    end
    
    if (length(poutLower) <i || ~poutLower(i).steadyState)
    poutLower(i) = Pout([dataFolder , 'mushyLayerLowC-lowerBranch-insulating/CR1.25RaC',num2str(RaC), ...
        'Le200ChiCubedPermeabilitypts', ...
        num2str(widths(i)),'-0/pout.0']);
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



actualWidths = widths*(0.2/128);
%actualWidths = actualWidths*L*2; %convert to metres
%actualWidths =actualWidths*100; %(cm)

plot(actualWidths(~isnan(fluxes)), fluxes(~isnan(fluxes))-1, 'x-');

%plot(actualWidths(~isnan(fluxesLower)), fluxesLower(~isnan(fluxesLower))-1, 'o-');


end

hold off;
axis([0 0.15 0 0.3]);

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
% h =gcf;
% set(h,'Units','Inches');
%  pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(h,['soluteFluxRaC',num2str(RaC),'CR1.25.pdf'],'-dpdf','-r0')

