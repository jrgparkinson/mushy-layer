function chimneyAnalyticSol

close all;

Le = 200;
Ra = 50;
dSldz =  1.0; % 4 observed, 0.5 reproduces results

% Values at the midpoint of the chimney
T0 = 0.5; %0.7
S0 = T0+0.12; %-0.6

Le_arr= 50:50:1000;
Le_arr = [200];
a_arr = Le_arr*NaN;

for i = 1:length(Le_arr)
    
[T, Sl, psi, a_arr(i), Sa, x] = analyticSolChannel(Le_arr(i), Ra, dSldz, T0, S0);

end


% figure();
% xFit = log10(1./Le_arr); yFit = log10(a_arr);
% p = polyfit(xFit, yFit, 1)
% plot(xFit, yFit);
% hold on;
% 
% yVals = polyval(p, xFit);
% plot(xFit, yVals, '--');
% hold off;
% 
% xlabel('log10(1/Le)'); ylabel('log10($a$)');
% 



h=figure();
set(h, 'Position', [200 200 1000 700]);
hold on;


plot(x, psi, '-');
plot(x, T, '-');
plot(x, Sl, '-');
plot([a_arr(end) a_arr(end)], [min([min(Sl), min(psi), min(T)]) max([max(Sl), max(T), max(psi)])], '--');

hold on;

legend('psi','T','Sl', 'a', 'Location', 'eastoutside');


end

