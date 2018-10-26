function [ T ] = temperature( H, S, params )
%TEMPERATURE Summary of this function goes here
%   Detailed explanation goes here

[ Hs,He,Hl ] = boundingEnergies( H, S, params );

chi = porosity(H ,S, params);

% initialize T;
T = H*0;
 
for i=1:length(H)
    if H(i) < Hs(i)
        T(i) = H(i)/params.cp;
    elseif H(i) <= He(i)
        T(i) = 0;
    elseif H(i) < Hl(i)
        T(i) = (S(i) - params.concRatio*(1-chi(i)))/(chi(i) + params.pc*(1-chi(i)));
    else
        T(i) = H(i) - params.stefan;
    end
end



end

