function confinementScaling

clear all; close all;

set(groot, 'defaultAxesFontName', 'times');
set(groot, 'defaultAxesLineStyleOrder', {'-','--',':'});
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLineLineWidth', 2.0);

CR = 1e-3:1e-3:1e1;
CR = logspace(-3,3,10); %[0.01 0.1 1.0 10];
St = 5;
thetaInf = 0.8;
thetai = 0.0;

CR_leg = {};
chi_leg = {};

confinementChi = logspace(-1,0,5); %[0.01 0.05 0.1 0.2 0.5 0.9];
confinementChi = [0.01 0.05 0.1 0.2 0.5 0.9];

confinementSize = NaN*ones(length(confinementChi), length(CR));

h=figure;
set(h, 'Position', [300 300 700 400]);
m=1;n=2;
subplot(m,n,1);
hold on;

for CR_i = 1:length(CR)
    [chi, theta, z, h] = solution(St, CR(CR_i), thetaInf, thetai);
    
    plot(chi, z);
    
    CR_leg{CR_i} =['$\mathcal{C}=',sprintf('%.2e',CR(CR_i)),'$'];
    
    
    % Determine where chi=0.1;
   
        
    for chi_i = 1:length(confinementChi)
        
        chi_leg{chi_i} = ['$\chi=',num2str( confinementChi(chi_i)),'$'];
        
    confinement_i = find(chi < confinementChi(chi_i));
    confinement_i = max(confinement_i);
    
    confinement_depth = z(confinement_i);
    
     try
    confinementSize(chi_i, CR_i) = h - confinement_depth;
    
    catch e
    
    end
  
    end
    
    
    
    
    
end

xlabel('$\chi$'); ylabel('$z$');
box on;

legend(CR_leg, 'Location', 'northwest');

hold off;

subplot(m,n,2);

plot(log10(CR), log10(confinementSize), 'x-');
xlabel('log$_{10}(\mathcal{C})$'); ylabel('log$_{10}(h_{confinement})$');
box on;

legend(chi_leg, 'Location', 'southeast');



end


function [chi, theta, z, h] = solution(St, CR, thetaInf, thetai)


h = mushZ(St, CR, thetaInf, thetai, 0);

spacing = min(CR/10, 0.01);
theta = -1:spacing:0;

z = mushZ(St, CR, thetaInf, thetai, theta);

chi = (CR-thetai)./(CR - theta);

end

function z = mushZ(St, CR, thetaInf, thetai, theta)

A = 0.5*(CR + thetaInf + St);
B = sqrt(A^2 - CR*thetaInf - St*thetai);

alpha = A + B;
beta = A-B;

z = ((alpha-CR)/(alpha-beta) )*log((alpha+1)./(alpha-theta)) + ((CR-beta)/(alpha-beta))*log((beta+1)./(beta-theta));



end