function [ h, thetaFinal, chi ] = directionSolidificationAnalyticSol(xGrid, params)
%DIRECTIONSOLIDIFICATIONANALYTICSOL Summary of this function goes here
%   Detailed explanation goes here

h = mushyZ(params.thetaInterface, params);

N = length(xGrid);
dtheta = (params.thetaInf-params.Te)/(10*N);

%Create a fine grid of thetas
thetaGrid = params.Te:dtheta:params.thetaInf;
%zMushCalc = NaN*thetaGrid;
zMushCalc = [];
thetaMushy = [];

for theta_i = 1:length(thetaGrid)
    theta = thetaGrid(theta_i);
    if theta < params.thetaInterface
        zMushCalc(end+1) = mushyZ(theta, params);
        thetaMushy(end+1) = theta;
    end
end


%Combine solutions in two different regions
for x_i =1:length(xGrid)
    x = xGrid(x_i);
    
    if x>=h
        %Interpolate
        %thetaTotal[z_i] =
        theta = interp1(zMushCalc,thetaMushy,x, 'spline');
        
        %thetKatz = (theta - Te)/DT;
        %solidFrac(z_i) = (TiL - T)/(DT*CompRatio + TiL - T);
        chi(x_i) = (params.concRatio - 1)/(params.concRatio - theta);
    else
        
        theta = params.thetaInf + (params.thetaInterface - params.thetaInf )*exp((x-h) * params.V);
        
        chi(x_i) = 1;
    end
    
    thetaFinal(x_i) = theta;
    
    %thetaKatz(x_i) = (theta - Te)/DT;
    %thetaWorster(x_i) = (theta-TiL)/DT;
    
    
end

%chi = 1-solidFrac;

end


function z = mushyZ(theta, params)

if params.stefan == 0
    z = 0;
    return;
end

A = 0.5*(params.concRatio + params.thetaInf + params.stefan);
B = sqrt(A*A - params.concRatio*params.thetaInf - params.stefan*params.thetaInterface);
alpha = A+B;
beta = A-B;

%z = (1/params.V)* (  ((alpha-params.concRatio) / (alpha-beta)) * log((alpha)/(alpha-theta)) + ...
%    ((params.concRatio - beta) / (alpha-beta)) * log((beta)/(beta-theta))   );

z = 1 - (1/params.V)* (  ((alpha-params.concRatio) / (alpha-beta)) * log((alpha)/(alpha-theta)) + ...
    ((params.concRatio - beta) / (alpha-beta)) * log((beta)/(beta-theta))   );
end

