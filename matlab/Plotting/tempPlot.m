figure();
hold on;
for CR_i = 1:length(CR_arr)
    
    fluxForCR = optimalFluxMat(CR_i, :);
    %Ra = Ra_arr(Ra_i);
    Pi = averagePermMat(CR_i, :);
    Pi = optimalPermMat(CR_i, :);
    chi = optimalPorosMat(CR_i, :);
    
    CR = CR_arr(CR_i);
    
    Racrit = Ra_crit(CR_i);
    
    Ra_star = (Ra_arr-Racrit)/Racrit;
    Ra_star = Ra_arr;
    
    actualFluxes= ~isnan(fluxForCR);
    if sum(actualFluxes) > 0
        %plot1 = plot((CR./Ra_star).^(1/2), Pi(actualFluxes));
        %plot1 = plot(log10(Ra_star(actualFluxes)), log10(Pi(actualFluxes)));
        plot1 = plot(log10(Ra_star(actualFluxes)/CR), log10(chi(actualFluxes)));
    end
end
hold off;
box on;
%xlabel('$(\mathcal{C}/\mathcal{R})^{1/2}$'); ylabel('$\Pi$');
xlabel('$\mathcal{R}^*/CR$'); ylabel('$\chi$');
title('Median permeability');
%ylim([0 100]);
