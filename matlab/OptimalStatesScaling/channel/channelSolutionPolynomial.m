function [T,S,psi] = channelSolutionPolynomial(Le,Ra,T0,S0,dSdz,doPlots,a)
close all;


alpha = 2*(T0-S0)/(a*dSdz);
beta = -Le*a^3/12 - a^3/6;
gamma = (1+0.25*Le*dSdz*Ra*a^3)/(Ra*Le*dSdz);

A3 = (a-alpha)/(beta+gamma);
A1 = (alpha+beta*A3)/a;

B3 = -0.5*Le*dSdz*a*A3;
B2 = (A1 + (1/6)*a^2*A3)*dSdz-0.5*a*B3;

A1 = 1;
A3 = -1/a^2;
B2 = 1/a;
B3 = -1/a^2;

x = linspace(0,a,100);

T=T0+0*x;
S = S0 + 0.5*B2*x.^2 + (1/6)*B3*x.^3;
psi = A1*x + (1/6)*A3*x.^3;

if doPlots
figure();
hold on;

plot(x,  T, '-');
plot(x, psi, '-');
plot(x, S, '-');

hold off;

legend('T', 'psi', 'S');
end

end