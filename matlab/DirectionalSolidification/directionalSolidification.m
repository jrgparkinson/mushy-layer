%clear all;
clf;

% Directional solidification analytic solution
% z is in real units here unless explicitly stated

%Physical constants
CompRatio = 1.4; %1.4
Tinf = 15;
Te = -20;
DT = 25;
TiL = Te + DT;
Stefan = 5.7; %5.7
V = 1 * 10^(-6);
diffusivity = 1.429*10^(-7);
TInterface = TiL;
height = 1; %metres

    %Scalings
    LWorster = diffusivity/V;
    LKatz = height; %*1.7
    
    numPoints = 256;
    
    dz = height/(numPoints);
    zGrid = (0:dz:height) + dz/2;
    zGrid = zGrid(1:end-1); 

    zKatz = zGrid/LKatz;        dzKatz = dz/LKatz;
    zWorster = zGrid/LWorster;  dzWorster = dz/LWorster;
    TGrid = (Te:(0.001):TInterface);
    
    V = 1*(10^(-6));
    
    [thetaKatz, thetaWorster, TTotal, solidFrac] = analyticSoln(zGrid, TGrid, ...
    TInterface, diffusivity, V, CompRatio, Stefan, Tinf, TiL, Te, DT, height);

    %figure(2);
    %plot(zKatz, thetaKatz);

    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Check solution
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dz = dzWorster;
    theta = thetaWorster;
    %nonDimV = V*(height/diffusivity); %Katz 2008/my model
    nonDimV = 1; %Worster 1991
    z = zWorster;
    
        
    %Lets load the data from Chombo model
    %data = csvread('stefanZeroAnalytic.data');
%     data = csvread('fullBm1Analytic.data');
%     Stefan = 5.7;
%     z = data(:,1);
%     theta = data(:, 2);
%     dz = z(2)-z(1);
%     solidFrac = data(:, 3);
    
    %figure(1);
    %plot(z, theta);
    
    H = theta + (1-solidFrac)*Stefan;
    %dthetadz = gradient(theta, dz);
    laplacian_theta = 4*del2(theta, dz);
    dHdz = gradient(H, dz);
    
    residual = laplacian_theta + nonDimV*dHdz;
    
    figure(5);
    hold on;
    plot(z, laplacian_theta);
    plot(z, nonDimV*dHdz);
    plot(z, residual);
    
    %figure(10);
    %plot(zKatz, H);
    
    
   

%{
figure(1);
hold on;
for i=1:4
    V = 2*(10^(-6))/(2^(i-1));
    [thetaKatz, solidFrac] = analyticSoln(zGrid, TGrid, ...
    TInterface, diffusivity, V, CompRatio, Stefan, Tinf, TiL, Te, DT, height);
    plot(zKatz, thetaKatz);
end
%}


% figure(2);
% hold on;
% for i=1:4
%     V = 2*(10^(-6))/(2^(i-1));
%     [thetaKatz, solidFrac] = analyticSoln(zGrid, TGrid, ...
%     TInterface, diffusivity, V, CompRatio, Stefan, Tinf, TiL, Te, DT, height);
%     plot(zKatz, solidFrac);
% end  


%Generate analytic soln for different grid resolutions
% for i=0:5
%     
%     % Grids - setup to same structure as chombo
%     numPoints = 32 * 2^i;
%     
%     dz = height/(numPoints);
%     zGrid = (0:dz:height) + dz/2;
%     zGrid = zGrid(1:end-1); 
% 
%     zKatz = zGrid/LKatz;
%     zWorster = zGrid/LWorster;
%     TGrid = (Te:(0.001):TInterface);
% 
%     fileID = fopen(strcat('analyticSoln', num2str(numPoints), '.data'),'w');
%     
%     V = 1*(10^(-6));
%     [thetaKatz, solidFrac] = analyticSoln(zGrid, TGrid, ...
%     TInterface, diffusivity, V, CompRatio, Stefan, Tinf, TiL, Te, DT, height);
% 
% %     printOut = zeros(length(zGrid), 3);
% %     for j=1:length(zGrid)
% %         printOut(j, :) = [zGrid(j), thetaKatz(j), solidFrac(j)]
% %     end
%     
%     printOut = [zGrid; thetaKatz; solidFrac];
%     
%     
%     fprintf(fileID,'%2.8f, %2.8f, %2.8f \n',printOut);
%  
%     fclose(fileID);
% end
% 
% plot(zKatz, solidFrac);
% 
