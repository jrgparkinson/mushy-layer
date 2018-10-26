%Load all output for Ra = 50 and compare convergence
%close all;
%clf;

%Define constants
NUSSELT_TEST = 'nusseltTest-';
BM2_FIXED = 'bm2-fixed-';
BM2_TGA_FIXED = 'bm2-TGA-fixed-';
BM2_NEWGRID_FIXED = 'bm2-newGrid-fixed-';
BM2_ONE_LEVEL = 'bm2-';
BM2_SIDEWALL_FIXED = 'bm2-sidewall-fixed-';
BM2_SIDEWALL_ONE_LEVEL = 'bm2-sidewall-';
BM2_SIDEWALL_PERIODIC_FIXED = 'bm2-sidewall-periodic-fixed-';
BM2_SIDEWALL_PERIODIC_ONE_LEVEL = 'bm2-sidewall-periodic-';

plot_title = containers.Map();

plot_title(NUSSELT_TEST) = 'nusselt test, 2 levels';
plot_title(BM2_ONE_LEVEL) = 'vertical heating, 1 level';
plot_title(BM2_FIXED) = 'vertical heating, 2 levels';
plot_title(BM2_NEWGRID_FIXED) = 'vertical heating, 2 levels (different grid)';
plot_title(BM2_TGA_FIXED) = 'vertical heating (TGA), 2 levels';
plot_title(BM2_SIDEWALL_FIXED) = 'sidewall, 2 levels';
plot_title(BM2_SIDEWALL_ONE_LEVEL) = 'sidewall, 1 level';
plot_title(BM2_SIDEWALL_PERIODIC_FIXED) = 'sidewall periodic, 2 levels';
plot_title(BM2_SIDEWALL_PERIODIC_ONE_LEVEL) = 'sidewall periodic, 1 level';

Nu_caltagirone = containers.Map({0,10,30,39,50,100,200,250}, [1 1 1 1  1.45 2.651 3.813 4.199]);
dim = 2;
frame = -1;

% Change these things!
output_dir = '/Users/parkinsonj/Dropbox/DPhil/Data/';
Ras = [50];
grid_res = [8, 16, 32];

dodgy_simulations = {{50, BM2_FIXED, 128, 'finished early'}, {100, BM2_SIDEWALL_ONE_LEVEL, 16, 'finished early'}};


flags = '';
%flags = 'adaptive-';
%flags = '2order-';
%flags = 'fixed-';
%flags = 'fixed-HOFluidInterp-';
%flags = 'HOFluidInterp-';


% Don't change anything else!
Nu_res_av = []; Nu_res_base = []; nusselt_error = [];


flags_list = {BM2_ONE_LEVEL};

figure_i = 1;

% Now plot the actual theta fields for some of these
Ra = Ras(1);
finest_res = grid_res(end);

% Here, the single level version MUST come first
flags_list = {BM2_ONE_LEVEL};
%flags_list = {BM2_ONE_LEVEL};
%flags_list = {BM2_SIDEWALL_ONE_LEVEL, BM2_SIDEWALL_FIXED};
%flags_list = {BM2_SIDEWALL_PERIODIC_ONE_LEVEL, BM2_SIDEWALL_PERIODIC_FIXED};
NUSSELT = -1; NUSSELT_ERROR = -2;
final_plot = NUSSELT_ERROR;
err_type = ChomboCompare.L2_ERR;

use_Nu_caltagirone = (strcmp(flags_list(1), BM2_FIXED) || strcmp(flags_list(1), BM2_ONE_LEVEL));
%use_Nu_caltagirone = false;

% determine if we have vertical or horizontal heating (for Nusselt
% calculation)
if (strcmp(flags_list(1), BM2_FIXED) || strcmp(flags_list(1), BM2_ONE_LEVEL) ...
        || strcmp(flags_list(1), NUSSELT_TEST))
    dir = 0;
else
    dir = 1;
end

fig = figure();
set(fig, 'Position', [500 500 800 1200])
%suptitle(['Steady state \theta error fields for Ra = ', num2str(Ra)]);
figure_i = 1;

nusselt_error = [];  nusselt = [];   L1_error = [];  L2_error = []; Max_error = []; Sum_error = [];

num_Res = length(grid_res); num_flags = length(flags_list);
num_Res = num_Res + 1; % this is a hack because I don't quite understand what's going on

fine_sol_name = [flags_list{1, 1}, 'Ra', num2str(Ra), '-',num2str(finest_res),'pts-steady'];
%fine_sol_name = ['errTest-exact'];
finest_single_level_solution = MushyLayerOutput(dim, frame, output_dir, ...
    fine_sol_name);
chCompare = ChomboCompare(finest_single_level_solution);

for flag_i=1:length(flags_list)
    figure_i = flag_i;
    
    for res_i = 1:length(grid_res)
        
        flags = flags_list{1, flag_i};
        res = grid_res(res_i);
        
        plot_prefix = [flags, 'Ra', num2str(Ra), '-',num2str(res),'pts-steady'];
        
        output = MushyLayerOutput(dim, frame, output_dir, plot_prefix);
              
        [Nu_av, Nu_base] = output.nusseltNumber(dir);
        nusselt(flag_i, res_i) = Nu_av;
        %error(res_i) = abs(Nu_av - Nu_caltagirone(Ra));
        
        subplot(num_Res+1, num_flags, figure_i);
        
%         output.plot(output.mlComps.theta);
        
        if output.exists()
            theta_comp = output.mlComps.theta;
            chCompare.diff(output);
            [L1, L2, Max, Sum] =  chCompare.err(output, output.mlComps.theta);
           %L1 = 0; L2=0; Max=0; Sum =0;
            L1_error(flag_i, res_i) = L1;
            L2_error(flag_i, res_i) = L2;
            Max_error(flag_i, res_i) = Max;
            Sum_error(flag_i, res_i) = Sum;
               
                
            compareData = chCompare.diff_output;
             compareData.plot(output.mlComps.theta);
           % comparison_error(flag_i, res_i) = NaN;
        else
            L1_error(flag_i, res_i) = NaN;
        end

       h = title({['Coarsest res = ', num2str(res)]; plot_title(flags)});
       P = get(h,'Position')    ;
       set(h,'Position',[P(1) P(2)+0.08 P(3)]);
       
        
        figure_i = figure_i + num_flags;
        
    end
    
    %Stick some convergence plots in the last row
    
    %     subplot(num_Res+1, num_flags, figure_i);
    %     plot( log(grid_res), log(error), 'x-');
    %     axis([2 6 -7 -1]);
    %     title('error convergence');
    
    
end

% Calculate error in Nusselt number, either against a value from the
% literature or against a high resolution single level run
% Call it Nu_true, even though it is only ever approximate
Nu_true = 0;
if use_Nu_caltagirone
    Nu_true = Nu_caltagirone(Ra);
else
    [Nu_true, ~] = finest_single_level_solution.nusseltNumber(dir);
end

for flag_i=1:length(flags_list)
    for res_i = 1:length(grid_res)
        nusselt_error(flag_i, res_i) = abs(nusselt(flag_i, res_i) - Nu_true);
    end
end


first_cell = (num_Res-1)*num_flags + 1;  last_cell = first_cell + num_flags - 1;

subplot(num_Res+1, num_flags, [first_cell last_cell]);
plotError(nusselt_error, NUSSELT_ERROR, plot_title, flags_list, grid_res, use_Nu_caltagirone)

%plotError(L1_error, ChomboCompare.L1_ERR, plot_title, flags_list, grid_res, use_Nu_caltagirone)

subplot(num_Res+1, num_flags, [first_cell+num_flags last_cell+num_flags]);
plotError(L2_error, err_type, plot_title, flags_list, grid_res, use_Nu_caltagirone)

