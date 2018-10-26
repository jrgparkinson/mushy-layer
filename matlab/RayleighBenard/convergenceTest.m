
close all;
clear all;

%Define constants
dim = 2;
frame = -1;

% Change these things!
%output_dir = '/Users/parkinsonj/Dropbox/DPhil/Data/'; %OS X
output_dir = '/home/parkinsonjl/convection-in-sea-ice/test/rayleighBenard/'; %Linux
prefix = 'rayleighBenard-';
grid_res = [4, 8, 16, 32];
finestRes = 64;
subcycled = true;

% Don't change anything else!
err_type = ChomboCompare.L2_ERR;

error = []; 
error_adaptive = [];


output_finest = MushyLayerOutput(dim, frame, output_dir, [prefix, num2str(finestRes), 'pts-steady'], subcycled);
data_finest = output_finest.dataForComp(output_finest.components.Temperature);

chComp = ChomboCompare(output_finest);

    for res_i = 1:(length(grid_res))
        
        res = grid_res(res_i);
        
        plot_prefix = [prefix, num2str(res),'pts-steady'];
        plot_prefix_adaptive = [prefix, 'adaptive-', num2str(res),'pts-steady'];
         
        output = MushyLayerOutput(dim, frame, output_dir, plot_prefix, subcycled);
        %output_adaptive = MushyLayerOutput(dim, frame, output_dir, plot_prefix_adaptive, subcycled);
        
        %data = output.dataForComp(output.components.Temperature);
        
        chComp.diff(output, [output.components.Temperature], [output.components.Temperature]);
        
        errField = chComp.diff_output.dataForComp(output.components.Temperature);
        %figure(); imagesc(errField); colorbar;
        [L1, L2, Max, Sum] = chComp.err(output, output.components.Temperature);
              
       % errField = output.dataForComp(output.components.yDarcyVelErr);
        %maxErr = max(max(abs(errField)));
        %nusselt(res_i) = Nu_av; 
        %error(res_i) = 1/res + rand()/100;
        
        error(res_i) = L1;
        
%         errField_adaptive = output_adaptive.dataForComp(output_adaptive.components.yDarcyVelErr);
%         if errField_adaptive
%         maxErr_adaptive = max(max(abs(errField_adaptive)));
%         %nusselt(res_i) = Nu_av; 
%         error_adaptive(res_i) = maxErr_adaptive;
%         else
%             error_adaptive(res_i) = NaN;
%         end
error_adaptive(res_i) = NaN;
        
    end
   
    
secondOrderSlope = error;
firstOrderSlope = error;
 for i=2:length(error)
     secondOrderSlope(i) = secondOrderSlope(i-1)/4; 
     firstOrderSlope(i) = firstOrderSlope(i-1)/2; 
 end

 
 
fig = figure();
set(fig, 'Position', [300 300 900 600]);

hold on;
plot(log(1./grid_res), log(error), 'x-');
plot(log(1./grid_res), log(error_adaptive), 'x-');
plot(log(1./grid_res), log(secondOrderSlope),'--');
plot(log(1./grid_res), log(firstOrderSlope),'--');

hold off;
xlabel('\Delta x'); 
ylabel('log(max err)');
title('Rayleigh-Benard convergence');
legend({'Max error', 'Max err adaptive','Second order comparison', '1st order'}, 'Location', 'southeast');

box on; grid on;
