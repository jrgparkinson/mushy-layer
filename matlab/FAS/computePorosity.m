function [ porosity ] = computePorosity( H, S, params )
%POROSITY Summary of this function goes here
%   Detailed explanation goes here

[ Hs,He,Hl ] = boundingEnergies( H, S, params );

porosity = 0*H;

for i=1:length(porosity)
    
if H(i) < Hs(i)
    porosity(i) = 0;
elseif H(i) <= He(i)
    porosity(i) = H(i)/params.stefan;
elseif H(i) < Hl(i)
    A = params.concRatio*(params.cp - 1) + params.stefan*(params.pc - 1);
    B = params.concRatio*(1-2*params.cp) + H(i)*(1-params.pc) + S(i)*(params.cp - 1) - params.pc*params.stefan;
    C = (params.concRatio - S(i))*params.cp + params.pc*H(i);
    
    porosity(i) = (-B - sqrt(B^2 - 4*A*C))/(2*A);
else
    porosity(i) = 1;
end

end

end

