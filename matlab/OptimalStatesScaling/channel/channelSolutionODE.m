function [T,S,psi,areturn,x] = channelSolutionODE(Le,Ra,T0,S0,dSdz,dTdz,doPlots,ainit, stepSize)


if doPlots
    %close all
end

%a = 0.012;
if nargin < 9
ainit = 0.16/sqrt(Le); % Works for Ra=50, T0=0.07, S0=0.12
ainit = 0.139/sqrt(Le); % Works for Ra=50, T0=0.05,S0=0.12

Ra_exponent = 0.68; %0.7
if Ra > 100
    Ra_exponent = 0.6; %0.6
elseif Ra > 200
    Ra_exponent = 0.05;
end


ainit = ainit*((50/Ra)^(Ra_exponent)); % Works for Ra=50
%doPlots = false;

stepSize = (ainit/500)*(50/Ra); %5e-6; %stepSize = 2e-4;
end

areturn = NaN;

if doPlots
figure();
title(['Le = ', num2str(Le), 'Ra=',num2str(Ra)]);
hold on;
end


for a=ainit:stepSize:(ainit + 1e3*stepSize)
fprintf('%1.7f, ', a);
%global S0; global A; global Le; global dSdz; global dPsiDx;

A = Ra*Le*dSdz;
options = bvpset('NMax',100000);
x = linspace(0,a,200);

momsolinit = bvpinit(x, [0 1 1]);
saltsolinit = bvpinit(x, [S0, 0]);
heatsolinit = bvpinit(x, [T0, 0]);

momSol = bvp4c(@(x,y)momEq(x,y,A), @(x,y)momBC(x,y), momsolinit, options);

psiVals = deval(momSol,x);

dPsiDx = psiVals(2, :);
saltSol = bvp4c(@(odeX,y)saltEq(odeX,y,Le,dSdz,x,dPsiDx), @(x,y)saltBC(x,y,S0), saltsolinit);
Svals = deval(saltSol, x);

heatSol = bvp4c(@(odeX,y)heatEq(odeX,y,dTdz,x,dPsiDx), @(x,y)heatBC(x,y,T0), heatsolinit);
Tvals = deval(heatSol, x);

S = Svals(1, :); T=Tvals(1, :); psi = psiVals(1,:);



if doPlots
plot(x,psiVals(1,:), '--')
%plot(x,-psi(2,:))
plot(x,Svals(1, :), '-')
plot(x,Tvals(1, :), '-.')
ylim([0 S0]);
drawnow;
end

fprintf('diff = %1.7f, ', abs(S(end)-T(end)));
if abs(S(end)-T(end)) < 1e-3
    areturn = a;
    fprintf('\nFound a = %1.7f \n', a);
    break;
end

end

if doPlots
%legend('psi', 'v', 'S');
hold off;
end
% psi''' = A (1-psi' )
% set psi1 = psi, psi2 = psi', psi3 = psi''
% psi1' = psi2
% psi2' = psi3
% psi3' = A (1-psi2)
% 


end

function dpsidx = momEq(x,psi,A)
dpsidx = [psi(2); psi(3); A*(1-psi(2))];
end

function res = momBC(psi_0,psi_a)
res = [psi_0(1)-0.1; ... %psi(0) = 0
    psi_a(2)-0.0;... %psi'(a) = 0
    psi_0(3)]; %psi''(0) = 0
end

% Salt eq
%S'' =  Le dSdz * (1-psi')
% S1 = S, S2 = S'
% S1' = S2
% S2' = Le dSdz * (1-psi')
function dSdx = saltEq(x, S,Le,dSdz,xin, dpsidx)
dpsi = interp1(xin, dpsidx, x);
dSdx = [S(2); Le*dSdz*(1-dpsi)];
end

function res = saltBC(S_0, S_a,S0)
%S0 = 0.12;
res = [S_0(1) - S0; ... % fixed value at middle of channel
    S_0(2)-0.0]; %no flux at middle of channel
end


% T eq
%S'' =  dTdz * (1-psi')
% S1 = S, S2 = S'
% S1' = S2
% S2' = dTdz * (1-psi')
function dSdx = heatEq(x, T,dTdz,xin, dpsidx)
dpsi = interp1(xin, dpsidx, x);
dSdx = [T(2); dTdz*(1-dpsi)];
end

function res = heatBC(T_0, T_a,T0)
%T0 = 0.07;
res = [T_0(1) - T0; ... % fixed value at middle of channel
    T_0(2)-0.0]; %no flux at middle of channel
end