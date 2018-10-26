function [ Hs,He,Hl ] = boundingEnergies( H, S, params )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

porosity_eutectic = 1-S/params.concRatio;

Hs = params.cp * max(0, (1/params.pc)*(S-params.concRatio));
He = porosity_eutectic*params.stefan;
Hl = params.stefan + S;

end

