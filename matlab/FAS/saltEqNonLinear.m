function [ f, df ] = saltEqNonLinear( S, Sn, H, params)
%SALTEQNONLINEAR Summary of this function goes here
%   Detailed explanation goes here

chi = porosity(H, S, params);
Sl = liquidSalinity(H,S,params);

diffusion = 0*S;
frameAdv = 0*S;

for i=2:length(S)-1
    chiPlusHalf = (chi(i+1)+chi(i))/2;
    chiMinusHalf = (chi(i-1)+chi(i))/2;
      
   diffusion(i) = (chiPlusHalf*(Sl(i+1)-Sl(i))  ...
       - chiMinusHalf*(Sl(i)-Sl(i-1) ))/(params.Le*params.dx^2);
   frameAdv(i) = params.V*(S(i-1)-S(i))/params.dx;
end

f = S - Sn - params.dt*(diffusion + frameAdv);


[ Hs,He,Hl ] = boundingEnergies( H, S, params );

dchidS = 0*S;
dSldS = 0*S;

for i=1:length(H)
    if H(i) < Hs(i)
        dchidS(i) = 0;
        dSldS(i) = 0;
    elseif H(i) <= He(i)
        dchidS(i) = 0;
        dSldS(i) = 0;
    elseif H(i) < Hl(i)
        %Assuming pc=0 for simplicity
        A = params.concRatio*(params.cp-1) + params.stefan*(-1);
        B = params.concRatio*(1-2*params.cp) + H(i) + S(i)*(params.cp-1) - params.pc*params.stefan;
        C = (params.concRatio-S(i))*params.cp;
        
        Bprime = params.cp;
        Cprime = -params.cp;
        
        dchidS(i) = (1/2*A) * (-Bprime + 0.5*(2*B^2*Bprime - 4*A*Cprime)/sqrt(B^2-4*A*C));
        
        dSldS(i) = ((1+params.concRatio*dchidS(i))*chi(i) - dchidS(i)*(S(i)-params.concRatio*(1-chi(i))) )/chi(i)^2 ;
    else
        dchidS(i) = 0;
        
        dSldS(i) = 1;
    end
end


saltFlux = 0*S;

for i=2:length(dchidS)-1
    %saltFlux(i) = (1/(2*params.dx))* dchidS(i) * ( (Sl(i+1)-Sl(i-1)) + chi(i)*(dSldS(i+1)-dSldS(i-1))  );
    saltFlux(i) = (1/(2*params.dx))* ( dchidS(i) * (Sl(i+1)-Sl(i-1)) + chi(i)*(dSldS(i+1)-dSldS(i-1))  );
end

%Now do salt fluxes at boundary using one sided diffs
saltFlux(1) = (1/params.dx)* dchidS(1) * ( (Sl(2)-Sl(1)) + chi(1)*(dSldS(2)-dSldS(1))  );
saltFlux(end) = (1/params.dx)* dchidS(end) * ( (Sl(end)-Sl(end-1)) + chi(end)*(dSldS(end)-dSldS(end-1))  );


saltDiffusion = 0*dchidS;
for i=2:length(dchidS)-1
   saltDiffusion(i) = ( saltFlux(i+1)-saltFlux(i-1) )/(2*params.dx);
end

% d(frame adv)/dH = 0
df = 1 - (params.dt/params.Le)*saltDiffusion;


end






