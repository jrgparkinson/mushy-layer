function [ p ] = computeNonDimVals( h, K0, S0 )

if nargin < 1
    h = 0.1;
end

if nargin < 2
    K0 = 1e-8;
end

if nargin < 3
    S0 = 30;
end

%COMPUTERM Summary of this function goes here
%   Detailed explanation goes here

p = getPhysicalConstants();

p.RaC = p.beta*p.rho_l*p.g*p.deltaC*h^3/(p.kappa_l * p.eta);
p.Da = K0/h^2;
p.RmC = p.RaC*p.Da;
p.Pr = p.eta/(p.rho_l*p.kappa_l);
p.timescale = h^2/p.kappa_l;
p.velScale = p.kappa_l/h;
p.Fsscale = p.deltaC*p.velScale;
p.CR = S0/(p.Se-S0);


fprintf('RaC = %1.5e \n', p.RaC);
fprintf(['Da = ', num2str(p.Da), '\n']);
fprintf(['RmC = ', num2str(p.RmC), '\n']);
fprintf(['Pr = ', num2str(p.Pr), '\n']);
fprintf(['Timescale = ', num2str(p.timescale), ' seconds = ',num2str(p.timescale/(3600*24)),' days \n']);
fprintf(['Velocity scale = ', num2str(p.velScale), ' m/s = ',num2str(p.velScale*100),' cm/s \n']);
fprintf(['Fs scale = ', num2str(p.Fsscale), ' m/s = ',num2str(p.Fsscale*1e4),' x 10^(-4) m/s \n']);

end

