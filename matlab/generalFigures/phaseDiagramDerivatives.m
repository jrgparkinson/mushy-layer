% Plot 3d phase diagram e.g. T(H, S)
close all;
clear all;



St_arr = 1:0.2:10;
CR = 0.25;
cp = 1;

averagedTdH = NaN*ones(1, length(St_arr));
averagedSldS = NaN*ones(1, length(St_arr));

for St_i = 1:length(St_arr)
    St = St_arr(St_i);
    [ averagedTdH(St_i), averagedSldS(St_i) ] = calcPhaseDiagDerivs( St, CR, cp );
end


figure()
hold on;
plot(St_arr, averagedTdH, 'x-');
plot(St_arr, averagedSldS, 'x-');
hold off;
box on;
title(['CR = ', num2str(CR), ', cp = ', num2str(cp)]);
legend({'dTdH', 'dSldS'});
xlabel('Stefan');




St = 5;
CR_arr = 0.01:0.02:1.0;

averagedTdH = NaN*ones(1, length(CR_arr));
averagedSldS = NaN*ones(1, length(CR_arr));

for cp_i = 1:length(CR_arr)
    CR = CR_arr(cp_i);
    [ averagedTdH(cp_i), averagedSldS(cp_i) ] = calcPhaseDiagDerivs( St, CR, cp );
end

figure()
hold on;
plot(CR_arr, averagedTdH, 'x-');
plot(CR_arr, averagedSldS, 'x-');
hold off;
box on;
title(['St = ', num2str(St), ', cp = ', num2str(cp) ]);
legend({'dTdH', 'dSldS'});
xlabel('CR');










St = 5;
CR = 0.25;
cp_arr = 0.1:0.1:5;

averagedTdH = NaN*ones(1, length(cp_arr));
averagedSldS = NaN*ones(1, length(cp_arr));

for cp_i = 1:length(cp_arr)
    cp = cp_arr(cp_i);
    [ averagedTdH(cp_i), averagedSldS(cp_i) ] = calcPhaseDiagDerivs( St, CR, cp );
end

figure()
hold on;
plot(cp_arr, averagedTdH, 'x-');
plot(cp_arr, averagedSldS, 'x-');
hold off;
box on;
title(['St = ', num2str(St), ', CR = ', num2str(CR)]);
legend({'dTdH', 'dSldS'});
xlabel('cp');