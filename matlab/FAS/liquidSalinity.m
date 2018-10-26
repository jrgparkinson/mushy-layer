function [ Sl ] = liquidSalinity( H, S, params )
%TEMPERATURE Summary of this function goes here
%   Detailed explanation goes here

[ Hs,He,Hl ] = boundingEnergies( H, S, params );

chi = porosity(H ,S, params);

% initialize T;
Sl = H*0;
 
for i=1:length(H)
    if H(i) < Hs(i)
        Sl(i) = 0;
    elseif H(i) <= He(i)
        Sl(i) = 0;
    elseif H(i) < Hl(i)
        Sl(i) = (S(i) - params.concRatio*(1-chi(i)))/(chi(i) + params.pc*(1-chi(i)));
    else
        Sl(i) = S(i);
    end
end



end

