% Predict how long it will take simulation to reach steady state
% /home/parkinsonjl/mnt/sharedStorage/optimalStates-highRes-maxSteadyState/CR1.08RaC200Le200ChiCubedPermeabilitypts60-0/diagnostics.out

function timeLeft = timeToSteady(diagsFile)
diags = getDiagnostics(diagsFile);

figure();
plot(log10(diags.time), log10(diags.dSdt), 'x');
xlabel('log10(time)');
ylabel('log10(dS/dt)');

% fit linear trend to log-log plot
p = polyfit(log10(diags.time), log10(diags.dSdt), 1);

hold on;
plot(log10(diags.time), p(2) + log10(diags.time)*p(1), '--');

% Solve for d/dt = 10^-3
criteria = -5;
logTime = (criteria-p(2))/p(1);
predictedFinish = exp(logTime);

predictedTimeRemaining = predictedFinish - diags.time(end);

dt = mean(diags.time(2:end) - diags.time(1:end-1));
wallClockPerTimestep = 1; % should calculate this properly rather than guessing 

predictedPhysicalTime = (predictedTimeRemaining/dt)*wallClockPerTimestep;

fprintf('Predicted run time remaining: %2.2f seconds = %2.2f hours \n', predictedPhysicalTime, predictedPhysicalTime/3600);


end