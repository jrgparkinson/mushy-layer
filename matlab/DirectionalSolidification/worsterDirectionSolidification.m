clear all;
clf;

CompRatio = 0.1;
Stefan = 5;
V = 1e-6;
diffusivity = 1.429*10^(-7);
numPoints=8048;
height = 1;

nonDimVWorster = 1;
LWorster = diffusivity/V;

nonDimVKatz = V*height/diffusivity;
LKatz = height;

%Non dimensionalisation
domainHeight = height/LWorster;


thetaInf = 2.0;
thetaEutectic = -1;
thetaInterface = 0;

A = 0.5*(CompRatio + thetaInf + Stefan)
B = sqrt(A^2 - CompRatio*thetaInf)
alpha = A+B
beta = A-B

theta = thetaEutectic:0.01:thetaInterface;
zCalc = theta*NaN;
zApprox = theta*NaN;

h = zMushWorster(CompRatio, thetaInf, Stefan, 0);

dz = domainHeight/(numPoints-1);
zGrid = 0:dz:domainHeight;
thetaInterp = zGrid;
thetaApprox = zGrid;
phi = zGrid;
phiApprox = zGrid;

h_i = round(1.05*h/dz);

boundaryWidth = log((alpha+1)/(alpha-CompRatio));

for theta_i = 1:length(theta)
    t = theta(theta_i);
    zCalc(theta_i) = zMushWorster(CompRatio, thetaInf, Stefan, t);
    zApprox(theta_i) = log((alpha+1)/(alpha-t));
end

for z_i = 1:length(zGrid)
    if (zGrid(z_i) <= h)
       thetaInterp(z_i) = interp1(zCalc,theta,zGrid(z_i), 'spline');
       phi(z_i) =  -thetaInterp(z_i)./(CompRatio-thetaInterp(z_i));
       
       thetaApprox(z_i) =  min(interp1(zApprox,theta,zGrid(z_i), 'spline'), -0.05);
       
       phiApprox(z_i) = -thetaApprox(z_i)./(CompRatio-thetaApprox(z_i));
    else
       thetaInterp(z_i) = thetaInf + (thetaInterface - thetaInf)*exp(h-zGrid(z_i));
       phi(z_i) = 0;
    end
    
   
    
    
end


nonDimV = nonDimVWorster;



%Check interpolation is done right (lines should overlap)
% figure(1);
% hold on;
% plot(zCalc, theta);
% plot(zGrid, thetaInterp);


%Check governing equation is satisfied (residual should = 0)
H = thetaInterp + Stefan.*(1-phi);

figure(3);
hold on;
plot(zGrid(1:h_i), phi(1:h_i));
plot(zGrid(1:h_i), thetaInterp(1:h_i));
plot(zGrid(1:h_i), thetaApprox(1:h_i));
plot(zGrid(1:h_i), phiApprox(1:h_i));
plot([boundaryWidth boundaryWidth], [-1 1]);

legend({'\phi', '\theta', '\theta_{approx}', '\phi_{approx}'});

laplacian_theta = 4*del2(thetaInterp, dz);
dHdz = gradient(H, dz);

residual = laplacian_theta + nonDimV*dHdz;




