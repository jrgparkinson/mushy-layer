function [ f, df ] = heatEqNonLinear( H, Hn, S, params)
%HEATEQNONLINEAR Summary of this function goes here
%   Detailed explanation goes here

T = temperature(H, S, params);
diffusion = 0*T;
frameAdv = 0*T;

for i=2:length(T)-1
   diffusion(i) = ( T(i+1) + T(i-1) - 2*T(i) )/params.dx^2;
   frameAdv(i) = params.V*(H(i-1)-H(i))/params.dx;
end

f = H - Hn - params.dt*(diffusion + frameAdv);


[ Hs,He,Hl ] = boundingEnergies( H, S, params );
chi = porosity(H,S,params);

dTdH = 0*T;

for i=1:length(H)
    if H(i) <= Hs(i)
        dTdH(i) = 1/params.cp;
    elseif H(i) <= He(i)
        dTdH(i) = 0;
    elseif H(i) < Hl(i)
        %Assuming pc=0 for now
        A = params.concRatio*(params.cp-1) + params.stefan*(-1);
        B = params.concRatio*(1-2*params.cp) + H(i) + S(i)*(params.cp-1) - params.pc*params.stefan;
        C = (params.concRatio-S(i))*params.cp;
        
        dChidt = (-1/2*A)*(1 - B/sqrt(B^2 - 4*A*C));
        
        dTdH(i) =(S(i)-params.concRatio)*(-1/chi(i)^2)*dChidt ;%(S(i) - params.concRatio*(1-chi(i)))/(chi(i) + params.pc*(1-chi(i)));
    else
        dTdH(i) = 1;
    end
end





laplace_g = 0*dTdH;
for i=2:length(dTdH)-1
   laplace_g(i) = ( dTdH(i+1) + dTdH(i-1) - 2*dTdH(i) )/params.dx^2;
end

% d(frame adv)/dH = 0
df = 1 - params.dt*laplace_g;


end






