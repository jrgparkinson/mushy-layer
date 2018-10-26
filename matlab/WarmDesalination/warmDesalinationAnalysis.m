% Analyse warm desalination fluxes
clear all; close all;

folder = getDataDir('springtime/');
%HBottom = {'6.02', '6.3', '6.05', '6.005', '6.0001', '6.5', '6.15'};
HBottom = [6.005, 6.05, 6.1];

leg = {};



%accumulatedSaltLoss = 

figure();

hold on;

for i = 1:length(HBottom)
    
    leg{end+1} = ['$T_{ocean} = ',num2str(HBottom(i)-6),'$'];
    
    %plot_prefix = ['CR1.15RaC10000Le200ChiCubedPermeabilityDa0.0R1.0pts128-domWidth5.0HBottom',HBottom{i},'-0'];
    plot_prefix = ['CR1.15RaC10000Le200ChiCubedPermeabilityDa0.0R1.0pts128-domWidth2.5HBottom',num2str(HBottom(i)),'profile1-0'];
       
    diagFile = [folder, plot_prefix, '/diagnostics.out'];
    diags = getDiagnostics(diagFile);
   
    temp = 0;
    
    plot(diags.time,max(-diags.Fs_bottom,1e-10) , '-'); %
    
    
end

hold off

set(gca, 'yscale', 'log');

xlabel('t (approximately days)');
ylabel('$F_s$ bottom (log scale)');
%ylabel('Average liquid salinity');
title('Positive flux corresponds to sinking salt');
box on;
legend(leg, 'Location', 'eastoutside');

%ylim([-0.1 0]);
%xlim([0 10]);