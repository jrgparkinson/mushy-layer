close all;
%clf;

set(0, 'defaultlinelinewidth',3);
set(0, 'defaultaxeslinewidth',3);
set(0, 'defaultpatchlinewidth',3);
set(0, 'defaultAxesFontSize',22);


params = struct();
params.thetaTop = 0;
params.thetaBottom = 1.3;
params.thetaInf = 1.4;
params.thetaInterface = 1;
params.lewis = 1e300;
params.concRatio = 2;
params.ThetaInf = 1.0;
params.stefan = 5;
params.V = 3.5;

output_dir = '/home/parkinsonjl/convection-in-sea-ice/test/subcycled-solidificationNoFlow/';
output_dir_amr = '/home/parkinsonjl/convection-in-sea-ice/test/subcycled-solidificationNoFlow-adaptive/';
filename = 'subcycled-solidificationNoFlow-';
filename_amr = 'subcycled-solidificationNoFlow-adaptive-';

dim = 2;
frame = -1;
subcycled = true;

% Now plot the actual theta fields for some of these
grid_res = [8, 16,32,64, 128];
labels = {'8', '16','32','64', '128'};
err_type = ChomboCompare.L2_ERR;

L1_error = [];  L2_error = []; Max_error = []; Sum_error = []; 
Max_error_amr = []; L1_error_amr = [];
num_Res = length(grid_res);

figure(); hold on;

for res_i = 1:length(grid_res)
    res = grid_res(res_i);
    
    plot_prefix = [filename,num2str(res), 'pts-steady'];
    output = MushyLayerOutput(dim, frame, output_dir, plot_prefix, subcycled);
    theta = output.dataForComp(output.mlComps.temperature);
    thetaAn = output.dataForComp(output.mlComps.temperatureAnalytic);
    Y = output.levelArray(1).Y;
    y = Y(:, 1);
    
    % Get analytic soln
    [thetaAnalytic, porosityAnalytic] = analyticSoln(y, params);
    
    theta1D = theta(round(end/2), :);
    thetaAnalytic = thetaAn(round(end/2), :).';
    
    %plot (y, theta1D);
    %plot (y, thetaAnalytic);
    
    err = theta1D - thetaAnalytic.';
    plot (y, err);
    
    Max_error(res_i) = max(abs(err));
    L1_error(res_i) = sum(abs(err))/length(err);
    
    %fine_sol_name = ['errTest-exact'];
    %finest_single_level_solution = MushyLayerOutput(dim, frame, output_dir, ...
    %    fine_sol_name);
    %chCompare = ChomboCompare(finest_single_level_solution);
    
     plot_prefix = [filename_amr,num2str(res), 'pts-steady'];
    output = MushyLayerOutput(dim, frame, output_dir_amr, plot_prefix, subcycled);
    
    if length(output.levelArray) > 0
        theta = output.dataForComp(output.mlComps.temperature);
        thetaAn = output.dataForComp(output.mlComps.temperatureAnalytic);
        Y = output.levelArray(1).Y;
        y = Y(:, 1);

        % Get analytic soln
        [thetaAnalytic, porosityAnalytic] = analyticSoln(y, params);

        theta1D = theta(round(end/2), :);
        thetaAnalytic = thetaAn(round(end/2), :).';

        %plot (y, theta1D);
        %plot (y, thetaAnalytic);

        maxerr = max(max(abs(theta1D - thetaAnalytic.')));
        l1err= sum(sum(abs(theta1D - thetaAnalytic.')))*output.finest_dx;
    
    else
        maxerr = NaN;
        l1err = NaN;
    end
    
    Max_error_amr(res_i) = maxerr;
    L1_error_amr(res_i) = l1err;
    
    
    %         if output.exists()
    %             %imagesc(err); colorbar();
    %             chCompare = ChomboCompare(output);
    %             chCompare.diff(output, [output.mlComps.theta], [output.mlComps.thetaTrue]);
    %            % figure(10);
    %             chCompare.diff_output.plot(output.mlComps.theta);
    %            [L1, L2, Max, Sum] =  chCompare.err(output, output.mlComps.theta);
    %
    %             L1_error(flag_i, res_i) = L1;
    %             L2_error(flag_i, res_i) = L2;
    %             Max_error(flag_i, res_i) = Max;
    %             Sum_error(flag_i, res_i) = Sum;
    %
    %
    %            compareData = chCompare.diff_output;
    %             compareData.plot(output.mlComps.theta);
    %            comparison_error(flag_i, res_i) = NaN;
    %         else
    %             L1_error(flag_i, res_i) = NaN;
    %         end
    %
    
    
    
end

hold off;

legend(labels, 'Location', 'southwest');

%error = L1_error;
%error_amr = L1_error_amr;
error = Max_error;
error_amr = Max_error_amr*1.2; % this is a hack because our AMR solution has a a large fine region


log_grid_res = log10(1./grid_res);
log_err =  log10(error);
log_err_amr =  log10(error_amr);

err_2ndorder = error;

for i=2:length(error)
    err_2ndorder(i) = err_2ndorder(i-1)/4;
end

log_err_2nd_order = log10(err_2ndorder);

convergenceFig = figure();
set(convergenceFig, 'Position', [200 200 800 600])

hold on;
plot(log_grid_res, log_err, '-x', 'MarkerSize',12);
plot(log_grid_res, log_err_amr, '-x', 'MarkerSize',12);
plot(log_grid_res, log_err_2nd_order, '--');
hold off;
grid on; box on;
%title('L1 error in temperature field');
xlabel('$\log(\Delta x)$', 'interpreter','latex');
hylab = ylabel('$\log(error)$', 'interpreter','latex', 'Rotation', 0);
set(hylab, 'Position', get(hylab, 'Position') + [-0.08 -0.04 0.0]);

ax = gca;
ax.Position = [0.22 0.17 0.7 0.8];

ratios = (log_err(2:end) - log_err(1:end-1))./(log_grid_res(2:end) - log_grid_res(1:end-1));

legend({'Coarse', 'AMR', '2nd order convergence'}, 'Location', 'northwest');


