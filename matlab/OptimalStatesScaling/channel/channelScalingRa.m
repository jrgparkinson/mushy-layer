close all;

getData = false;

if getData
Ra = 50;
dSdz = 4.0; dTdz = -4.0;
S0 = 0.12; T0 = 0.05; % previous results mainly for T0=0.07

%Le = [50 75 100 150 200 300 400 500 750 1000];
Le = 200;

% Looks a(Ra) looks like a~1/Ra for small Ra, and 1/(Ra^n) where n < 1 for
% large Ra

% For Ra too small, we don't find a channel
Ra = [50 45 40 39 38 37];
Ra = [50 100 200 300 400 500 600 800 1000];
a_guess = [0.00981 0.00613 0.00411 0.0033 0.002835 0.00253 0.0023 0.001986 0.001772];
%Ra = [800 1000];
%a_guess = [0.001986 0.001772];
%Ra = [400];
%a_guess = [0.002835];
delta_a = 1e-6;
    

a = Ra*NaN;
for i=1:length(Ra)
    fprintf('\nRa = %d \n', Ra(i));
    
    delta_a = 5e-6*(50/Ra(i));
    %delta_a = 1e-6;
    
    [T,S,psi,a(i)] = channelSolution(Le,Ra(i),T0,S0,dSdz,dTdz,true,a_guess(i), delta_a);
  
end
end

figure()
xVals = log10(1./Ra); yVals = log10(a);

%0.5 gradient line
ySqrt = yVals;
for i=length(yVals)-1:-1:1
   ySqrt(i) = ySqrt(i+1)*((xVals(i)/xVals(i+1))^(0.5));
end

xFit = xVals(~isnan(yVals)); yFit = yVals(~isnan(yVals));
p = polyfit(xFit, yFit, 1);
yLinear = polyval(p, xFit);
hold on;
plot(xVals, yVals, 'x');
%plot(xFit, yLinear, '-');
plot(xVals, ySqrt, '-');
hold off;
xlabel('log 10 (1/Ra)'); ylabel('log10(a)');
legend('Data', ['y=',num2str(p(1)), 'x + C'], 'Location', 'southeast');

