

getData = true;

if getData
Ra = 50;
dSdz = 4.0; dTdz = -4.0;
S0 = 0.12; T0 = 0.07;

Le = [50 75 100 150 200 300 400 500 750 1000];
a = Le*NaN;
for i=1:length(Le)
    fprintf('\nLe = %d \n', Le(i));
    
    [T,S,psi,a(i)] = channelSolution(Le(i),Ra,T0,S0,dSdz,dTdz,true);
  
end
end

figure()
xVals = log10(1./Le); yVals = log10(a);
p = polyfit(xVals, yVals, 1);
yFit = polyval(p, xVals);
hold on;
plot(log10(1./Le), log10(a), 'x');
plot(xVals, yFit, '-');
hold off;
xlabel('log 10 (1/Le)'); ylabel('log10(a)');
legend('Data', ['y=',num2str(p(1)), 'x + C'], 'Location', 'southeast');

