function [ p ] = getPhysicalConstants( )
%GETPHYSICALCONSTANTS Summary of this function goes here
%   Detailed explanation goes here

p.kappa_l = 1.25e-7; %m^2 s^-1
p.eta = 1.55e-3; %kg m^-1 s^-1
p.deltaC = 230; %g/kg
p.rho_l = 1028; %kg m^-3
p.beta = 7.86e-4; % (g/kg)^-1
p.g = 9.8; %m/s
p.alpha = 3.87e-5; % K^-1
p.Se = 230; % Eutectic salinity

end

