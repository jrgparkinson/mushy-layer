
close all;
%clear all;

%Define constants
dim = 2;
frame = -1;

% Change these things!
%output_dir = '/Users/parkinsonj/Dropbox/DPhil/Data/'; %OS X
output_dir = '/home/parkinsonjl/gyre/mushyLayer-insulating/'; %Linux
prefix = 'mushyLayer-insulating-CR5.0RaC80.0Le150.0ChiCubedPermeabilitypts';
steadySuffix = '-steady';

doLimiting = true;
limitingString = '';

if doLimiting
output_dir = '/media/parkinsonjl/FREECOM HDD/mushyLayerConvTest-limiting--insulating/';
prefix = 'mushyLayerConvTest-limiting--insulating-CR5.0RaC100.0Le200ChiCubedPermeabilitypts';
limitingString  ='withLimiting';
else
output_dir = '/media/parkinsonjl/FREECOM HDD/mushyLayerConvTest-insulating/';
prefix = 'mushyLayerConvTest-insulating-CR5.0RaC100.0Le200ChiCubedPermeabilitypts';
limitingString  ='noLimiting';
end


steadySuffix = '--steady';
grid_res_all = [16, 32, 64, 128, 256];
finestRes = grid_res_all(end);
grid_res = grid_res_all(1:end-1);
subcycled = true;

% Don't change anything else!
err_type = ChomboCompare.L2_ERR;

error = []; 
error_adaptive = [];


output_finest = MushyLayerOutput(dim, frame, output_dir, [prefix, num2str(finestRes), steadySuffix], subcycled);
data_finest = output_finest.dataForComp(output_finest.components.Temperature);

chComp = ChomboCompare(output_finest);

    for res_i = 1:(length(grid_res))
        
        res = grid_res(res_i);
        
        plot_prefix = [prefix, num2str(res),steadySuffix];
        plot_prefix_adaptive = [prefix, 'adaptive-', num2str(res),steadySuffix];
         
        output = MushyLayerOutput(dim, frame, output_dir, plot_prefix, subcycled);
        %output_adaptive = MushyLayerOutput(dim, frame, output_dir, plot_prefix_adaptive, subcycled);
        
        %data = output.dataForComp(output.components.Temperature);
        comp = output.components.Temperature;
        comp = output.components.Liquidconcentration;
        %comp = output.components.Permeability;
        
        chComp.diff(output, [comp], [comp]);
        
        errField = chComp.diff_output.dataForComp(comp);
        figure(); imagesc(errField.'); set(gca,'YDir','normal'); colorbar;
        set(gca, 'XTickLabel', ''); set(gca, 'YTickLabel', '');
        print([num2str(res),limitingString], '-dpng');
        [L1, L2, Max, Sum] = chComp.err(output, comp);
              
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
plot(log10(1./grid_res), log10(error), 'x-');
plot(log10(1./grid_res), log10(error_adaptive), 'x-');
plot(log10(1./grid_res), log10(secondOrderSlope),'--');
plot(log10(1./grid_res), log10(firstOrderSlope),'--');

hold off;
xlabel('log(\Delta x)'); 
ylabel('log(Error)');
title('Mushy layer convergence');
legend({'Error', 'Max err adaptive','Second order comparison', '1st order'}, 'Location', 'southeast');



box on; grid on;

print(['ConvergenceTest', limitingString], '-dpng');