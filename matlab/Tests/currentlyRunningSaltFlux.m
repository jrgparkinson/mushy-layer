clear all;

data_dir = getDataDir('optimalStates-Brinkman/');

folders = {'CR1.1RaC125Le200ChiCubedPermeabilitypts64-unsteadyVerySmallFlux',
    'CR11.0RaC4000Le200KozenyPermeabilitypts80-steady'   };

folder = folders{2};

d = getDiagnostics([data_dir, folder, '/diagnostics.out']);

t = d.time;
F = -1 - d.Fs_vertical_av;

average_dt = mean(t(2:end)-t(1:end-1));

dFdt = NaN*F;
for i=2:length(F)
    dFdt(i) = abs((F(i)-F(i-1))/(t(i)-t(i-1)));
end

d2Fdt = NaN*dFdt;
for i=2:length(d2Fdt)
    d2Fdt(i) = average_dt*abs((dFdt(i)-dFdt(i-1))/(t(i)-t(i-1)));
end

h = figure();
set(h, 'Position', [300 300 900 600]);
plot(t, F, 'x');
title({folder, ''});
xlabel('time');
ylabel('F');
box on; grid on;

N = round(length(t)/2);

p = polyfit(t(end-N :end), F(end-N :end), 2);
predicted_t = t;
dt = t(end) - t(end-1);
for i=1:length(t)
    predicted_t(end+1) = predicted_t(end) + dt;
end

hold on;
predictedFluxPlot = plot(predicted_t, polyval(p, predicted_t), '-');

if max(F) < 0
    ylim([-max(abs(F)) max(abs(F))]);
else
    deltaF = max(F) - min(F);
    
    ylim([min(F)-0.1*deltaF max(F) + 0.1*deltaF]);
end
hold off;

ml = getFinalPlotFile([data_dir, folder]);
frameAdv = 1.0;
Le = 200;
katzUnits = true;

calculatedF = ml.computeVerticalSoluteFlux(frameAdv, Le, katzUnits);


hold on;
calculatedPlot = plot(predicted_t, predicted_t*0 + calculatedF, '-');
hold off;

%yyaxis right;
%hold on
%dFdtPlot = plot(t, log10(abs(dFdt./F)), 'b-');
%d2FdtPlot = plot(t, log10(abs(d2Fdt./F)), 'r-');
%hold off
%ylabel('$log10(dF/dt)$');

%legend([predictedFluxPlot, calculatedPlot, dFdtPlot, d2FdtPlot], ...
%{'Predicted flux', 'Computed Flux', '$(1/F)dF/dt$', '$(1/F)d2F/dt2$'});
