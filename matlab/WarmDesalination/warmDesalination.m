% Jardon figure 1


function warmDesalination
close all;

T = linspace(-4,-1,100);

h = 0.4; % metres
h = 1;
Si = 5;


figure();

leg = {};

hold on;
for Si = 5:1:10
[Ra,~] = computeRa(T, Si, h);
plot(T, Ra, '-');
%semilogy(T, Ra, '-');
leg{end+1} = ['S=',num2str(Si)];

end

plot(T, 10+0*T, 'k--');

hold off
set(gca, 'yscale', 'log');
ylabel('Ra');

yyaxis right;

hold on;
for Si = 5:1:10
[~,chi] = computeRa(T, Si, h);
plot(T, chi, '-');
%semilogy(T, Ra, '-');
%leg{end+1} = ['S=',num2str(Si)];

end

hold off

ylabel('$\chi$');

xlabel('T');

legend(leg, 'Location', 'northwest');

box on;

end

function [Ra,chi] = computeRa(T, Si, h)

p= getPhysicalConstants();
Sbr = S(T);
mu = 0.054; %liquidus slope

%Sbr = -T/mu;
Sw = 30;

chi = -mu*Si./T;
Pi = computePermeability(chi);

p.kappa_l = 1.2e-7; %m^2 s^-1
p.eta = 2.55e-3;

Ra = p.g*h*p.rho_l*p.beta.*(Sbr-Sw).*Pi./(p.kappa_l*p.eta);

%Ra = g*(h-z);
end

function Pi = computePermeability(chi)

Pi = 1e-17*(1e3*chi).^(3.1);

end

function Sbr = S(Ti)
Sbr = -1.2-21.8*Ti-0.919*Ti.^2-0.0178*Ti.^3;
end