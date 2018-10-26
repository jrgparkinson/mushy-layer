function [ p ] = computeWettlauferNonDimVals( S0, Ttop )

% note that S0 is a wt%, where 1wt% = 10g/kg
% Ttop in degrees C

if nargin < 1
    S0 = 3.5;
end

if nargin < 2
    Ttop = -20;
end
p = getPhysicalConstants();
h = 0.01; % length scale 10cm
K0 = 1e-8; %1e-8
d = 3e-4; %3mm middleton hele-shaw cell

Ce = 230;
Te = -23;
Ci = S0*10;
deltaC = Ce - Ci;
deltaT = Tl(Ci) - Te;

Tbottom = -2;

%Params to compute
p.thetaTop = (Ttop - Te)/deltaT;
p.thetaBottom = (Tbottom-Te)/deltaT;
p.Thetai = (Ci - Ce)/deltaC;
p.CR = Ce/deltaC;
p.RaS = p.beta*p.rho_l*p.g*deltaC*h^3/(p.kappa_l * p.eta);

p.Da = K0/h^2;
p.RmS = p.Da*p.RaS;
p.Le = 200;
p.Pr = p.eta/(p.rho_l*p.kappa_l);
p.timescale = h^2/p.kappa_l;
p.St = 5;
p.chiTop = computePorosityMush(p, p.thetaTop, p.Thetai);
p.HTop = p.St*p.chiTop + p.thetaTop;
p.HBottom = p.St + p.thetaBottom;

p.RaT = p.alpha*p.rho_l*p.g*deltaC*h^3/(p.kappa_l * p.eta);
p.RmT = p.RaT*p.Da;

%temp1 = p.Rm/deltaC;
%temp2 = p.Ra/(deltaC);


p.heleShawPerm = d^2/(12*K0);
p.reluctance = 1/p.heleShawPerm;


end


function chi = computePorosityMush(p, T, S)
Sl = -T;
chi = (S+p.CR)/(Sl+p.CR);
end

function T = Tl(Ci)
T = - 0.1*Ci;
end
