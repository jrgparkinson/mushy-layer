function [T,S,psi,x] = channelSolutionSin(Le,Ra,T0,S0,dSdz,dTdz,doPlots,a)

if nargin < 8
    a = channelWidth(Le,Ra,T0,S0,dSdz,dTdz);
end

gamma = Le*dSdz;
beta = sqrt(Ra*gamma);

alpha = Le/100000;

% For computing a, need
a_rhs = beta^2 * (S0-T0)/(dTdz-Le*dSdz);

%a_pred_res = sec(beta*a)*(1-cos(beta*a)) - a_rhs;
%a = acos(1/(a_rhs+1)) / beta

x = linspace(0,a,100);

%sec(beta*a)

T=T0 + (dTdz*sec(beta*a)/beta^2)*(1-cos(beta*x));
S = S0 + (gamma*sec(beta*a)/beta^2)*(1-cos(beta*x)); % + (alpha*x.^2)/a^2;
psi = x - (sec(beta*a)/beta)*sin(beta*x);
dpsidx = 1 - sec(beta*a)*cos(beta*x);
V = -dpsidx;


if doPlots
    
    figure();
    hold on;
    plot(x,S,'-');
    plot(x,psi,'-');
    plot(x,T,'-');
   %  plot(x,V,'-');
     
    legend('S', 'psi', 'T');
    %legend('psi', 'V');
    hold off;
   
    
end



end