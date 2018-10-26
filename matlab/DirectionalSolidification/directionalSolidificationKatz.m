clear all;
clf;

CompRatio = 1.4;
Stefan = 5.7;
V = 1e-6;
diffusivity = 1.429*10^(-7);
numPoints=512;
height = 1;


nonDimV = V*height/diffusivity;
LKatz = height;

%Non dimensionalisation
domainHeight = height/LKatz;


thetaInf = 1.4;
thetaEutectic = 0;
thetaInterface = 1;


theta = thetaEutectic:0.01:thetaInterface;
zCalc = theta*NaN;

h = zMushKatz(CompRatio, thetaInf, Stefan, diffusivity, nonDimV, height, thetaInterface)

dz = domainHeight/(numPoints-1);
zGrid = 0:dz:domainHeight;
thetaInterp = zGrid;
phi = zGrid;


for theta_i = 1:length(theta)
    t = theta(theta_i);
    zCalc(theta_i) = zMushKatz(CompRatio, thetaInf, Stefan, diffusivity, nonDimV, height, t);
end

for z_i = 1:length(zGrid)
    if (zGrid(z_i) <= h)
       thetaInterp(z_i) = interp1(zCalc,theta,zGrid(z_i), 'spline');
       phi(z_i) =  (1-thetaInterp(z_i))/(CompRatio-thetaInterp(z_i));
    else
       thetaInterp(z_i) = thetaInf + (thetaInterface - thetaInf)*exp(nonDimV*(h-zGrid(z_i)));
       phi(z_i) = 0;
    end
end

% Check interpolation is done right (lines should overlap)
% figure(1);
% hold on;
% plot(zCalc, theta);
% plot(zGrid, thetaInterp);


%Check governing equation is satisfied (residual should = 0)
H = thetaInterp + Stefan.*(1-phi);

figure(3);
hold on;
plot(zGrid, phi);
plot(zGrid, thetaInterp);

laplacian_theta = 4*del2(thetaInterp, dz);
VdHdz = nonDimV * gradient(H, dz);

residual = laplacian_theta + VdHdz;

figure(5);
hold on;
plot(zGrid, laplacian_theta);
plot(zGrid, VdHdz);
plot(zGrid, residual);

legend('lap(theta)','v*dH_dz', 'residual')



