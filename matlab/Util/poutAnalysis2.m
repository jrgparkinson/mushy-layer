%clear all;
close all;

baseDir = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/';

% pout32Sponge = Pout([baseDir , 'output-sponge/pout.0']);
% pout32NoSponge = Pout([baseDir , 'output-noSponge/pout.0']);
% poutAMR = Pout(['/home/parkinsonjl/gyre/mushyLayer-periodic-adaptive/', ...
%     'CR5.0RaC80.0Le200.0ChiCubedPermeabilityDa1e-07-3/pout.0']);
% poutAMR = Pout([baseDir , 'output-AMR/pout.0']);
% 
% poutAMR2 = Pout([baseDir, 'pout.0']);

baseDir2 = '/home/parkinsonjl/gyre/mushyLayer-insulating/';
%pout32 = Pout([baseDir2, 'CR5.0RaC80.0Le200.0ChiCubedPermeabilitypts32-0/pout.0']);



% pout(1) = Pout('/home/parkinsonjl/gyre/mushyLayer-insulating/CR5.0RaC80.0Le20.0ChiCubedPermeabilitypts16-2/pout.0');
% pout(2) = Pout('/home/parkinsonjl/gyre/mushyLayer-insulating/CR5.0RaC80.0Le20.0ChiCubedPermeabilitypts32-1/pout.0');
% pout(3) = Pout('/home/parkinsonjl/gyre/mushyLayer-insulating/CR5.0RaC80.0Le20.0ChiCubedPermeabilitypts64-1/pout.0');

poutNames = {'16-1proc', '32', '64', '128', '256'};
gyre4base = '/home/parkinsonjl/gyre4/convection-in-sea-ice/test/';
testDir  = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/run/';
%dataFolder = '/home/parkinsonjl/data/';
dataFolder = '/media/parkinsonjl/FREECOM HDD/';

pout(1) = Pout([dataFolder , 'mushyLayerLowC-periodic/CR1.25RaC200Le200ChiCubedPermeabilitypts128-0/pout.0']);
pout(2) =  Pout([dataFolder , 'mushyLayerLowC-periodic/CR1.25RaC250Le200ChiCubedPermeabilitypts128-1/pout.0']); 
 %pout(3) = Pout([dataFolder , 'mushyLayerLowC-periodic/CR1.25RaC275Le200ChiCubedPermeabilitypts128-0/pout.0']);
 pout(3) = Pout([dataFolder , 'mushyLayerLowC-periodic/CR1.25RaC300Le200ChiCubedPermeabilitypts128-0/pout.0']);
 %pout(5) =  Pout([dataFolder , 'mushyLayerLowC-periodic/CR1.25RaC500Le200ChiCubedPermeabilitypts128-1/pout.0']);
%pout(6) = Pout([dataFolder , 'mushyLayerLowC-periodic/CR1.25RaC750Le200ChiCubedPermeabilitypts128-0/pout.0']);
%pout(7) = Pout([dataFolder , 'mushyLayerLowC-periodic/CR1.25RaC1000Le200ChiCubedPermeabilitypts128-0/pout.0']);
 
 fluxScale = 0.5;
 
 % poutNames = {'Ra200', 'Ra250', 'Ra275', 'Ra300', 'Ra500',  'Ra750', 'Ra1000'};
 Ra_arr = [200, 250, 300];
  poutNames = {'Ra200', 'Ra250', 'Ra300'};
  lineStyles = {'-', '--',':', '-.'};
lineColours = {'r', 'b', 'g', 'm', 'black'};
lineColours = {'m', 'c', 'r', 'g', 'b', 'k'};
h = figure();
set(h, 'Position', [200 200 1400 800]);
hold on;


legendStr = {};

for pouti = 1:length(pout)
    
    % Don't plot if we haven't got data
    if length(pout(pouti).times) == 0
        continue;
    end

    smoothedSponge = NaN*pout(pouti).fluxSponge;

    smoothSize = 30;
    for i = 1 + smoothSize/2 : length(smoothedSponge) - smoothSize/2

        sample = pout(pouti).fluxSponge(i-smoothSize/2:i+smoothSize/2);
        smoothedSponge(i) = sum(sample)/length(sample);
    end

    col_i = mod(pouti, length(lineColours));
    
    
    %colour = lineColours{col_i+1};
    colour = lineColours{col_i+1};
    linestyle = lineStyles{floor(pouti/length(lineColours))+1};
    
    
              %  plot(pout(pouti).times, fluxScale*pout(pouti).fluxTop, [lineStyles{2}, colour]);  
              % legendStr{end+1} = [poutNames{pouti}, ' top'];
              
                plot(pout(pouti).times, (fluxScale*pout(pouti).fluxBottom)-1.0,  [linestyle, colour]);
               legendStr{end+1} = [poutNames{pouti}];
            
            
end
            legend(legendStr, ...
               'Location', 'bestoutside' );
           %legend({'Top', 'Sponge top (averaged)'}, ...
           %    'Location', 'northeast' );
        
           title('Vertical solute flux at bottom of domain');
           
           xlabel('time');
           ylabel('F_s - V');
           
           grid on; box on;
           
           
           % Compare steady state fluxes at top of domain
           topSteadyFluxes = [];
           %res = [];
           Ra = [];
           for i=1:length(pout)
               if pout(i).steadyState
                    topSteadyFluxes(end+1) = fluxScale*pout(i).fluxTop(end) - 1; 
                    %res(i) = 1/(2^(i+2));
                    Ra(end+1) = Ra_arr(i);
                   
               end
           end
           
           figure();
           plot(Ra, topSteadyFluxes, 'x-');
           xlabel('Ra');
           ylabel('steady state flux');
           box on; grid on;

% For saving figure as PDF
% h =gcf;
% set(h,'Units','Inches');
% pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(h,'soluteFlux.pdf','-dpdf','-r0')

         