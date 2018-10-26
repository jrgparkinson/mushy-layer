% Assuming a purely diffusive solution with no salt diffusion, determine
% the far field temperature which would give the enforced enthalpy at the
% bottom of the domain, HBottom

% Based on the diffusive solution given by Worster (1991)

% It is assumed that the eutectic is at z=1, with the ocean extending down
% to z = -infinity
% Note that in most of my simulations, zbottom is 0.6 as the domain
% height is 0.4
% e.g. 
% E.g. computeFarFieldT(6.005,0.6,1.25,5,1)


function thetaInf = computeFarFieldT(HBottom, zbottom, CR, St, V)


Tinterface = 1.0;
TBottom = HBottom - St;

% Try different theta inf until one works

thetaInf = TBottom + 0.001;

while thetaInf < 10.0

A = 0.5*(CR + St + thetaInf);
B = sqrt(A^2 - CR*thetaInf - St);
alpha = A +B;
beta = A-B;

hV = V - ((alpha-CR)/(alpha-beta))*log((alpha)/(alpha-Tinterface)) - ((CR-beta)/(alpha-beta))*log(beta/(beta-Tinterface));
h = hV/V;

TBottomPredicted = thetaInf + (Tinterface - thetaInf)*exp(V*(zbottom-h));

fprintf('theta inf = %1.5f, predicted h = %1.5f, Predicted T bottom = %1.5f \n', thetaInf, h, TBottomPredicted);

if abs(TBottomPredicted - TBottom) < 0.01
    
    break;
end

thetaInf = thetaInf + 0.001;


end


% Plot the solution (in the liquid region) to check it makes sense

z = linspace(-5, h, 100);
T = thetaInf + (Tinterface - thetaInf)*exp(V*(z-h));

figure();
plot(z, T);
xlabel('z');
ylabel('T');






end