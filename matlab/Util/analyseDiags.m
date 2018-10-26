clear all; close all;
diagFile = ['/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/diagnostics.out'];

analyticNu = 3.07;

d = getDiagnostics(diagFile);

% Plot steady state 
% figure();
% plot(log10(d.dUdt));
% hold on;
% plot(log10(d.dSdt));
% plot(log10(d.dTdt));
% hold off;
% 
% legend('dU/dt', 'dS/dt', 'dT/dt', 'Location', 'eastoutside');

figure();
plot(d.time, d.Nusselt)
xlabel('$t$');
ylabel('Nu');

leg = {};
leg{1} = 'Calculated';

if ~isnan(analyticNu)
    hold on
    plot(d.time, d.time*0 + analyticNu);
    hold off
    leg{2} = 'Analytic';
end

legend(leg, 'Location', 'eastoutside');

