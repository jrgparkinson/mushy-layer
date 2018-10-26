%Load all output for Ra = 50 and compare convergence
close all;
%clear all;

%Define constants

Nu_caltagirone = containers.Map({0,10,30,39,50,100,200,250}, [1 1 1 1  1.45 2.651 3.813 4.199]);
dim = 2;
frame = -1;

% Change these things!
%output_dir = '/Users/parkinsonj/Dropbox/DPhil/Data/'; %OS X
output_dir = '/home/parkinsonjl/Dropbox/DPhil/Data/'; %Linux
prefix = 'bm2-';
Ra = 50;
grid_res = [8, 16, 32, 64, 128];
subcycled = false;

% Don't change anything else!
Nu_res_av = []; Nu_res_base = []; nusselt_error = [];

err_type = ChomboCompare.L2_ERR;


nusselt_error = []; 



    for res_i = 1:length(grid_res)
       
        res = grid_res(res_i);
        
        plot_prefix = [prefix, 'Ra', num2str(Ra), '-',num2str(res),'pts-steady'];
        
        output = MushyLayerOutput(dim, frame, output_dir, plot_prefix, subcycled);
              
        [Nu_av, Nu_base] = output.nusseltNumber();
        %nusselt(res_i) = Nu_av; 
        nusselt_error(res_i) = abs(Nu_base - Nu_caltagirone(Ra))/Nu_caltagirone(Ra);
        
    end
   
    
secondOrderSlope = nusselt_error;
firstOrderSlope = nusselt_error;
 for i=2:length(nusselt_error)
     secondOrderSlope(i) = secondOrderSlope(i-1)/4; 
     firstOrderSlope(i) = firstOrderSlope(i-1)/2; 
 end

 
 
fig = figure();

hold on;
plot(log(1./grid_res), log(nusselt_error), 'x-');
plot(log(1./grid_res), log(secondOrderSlope));
plot(log(1./grid_res), log(firstOrderSlope));

hold off;
xlabel('\Delta x'); 
ylabel('log(Nu - Nu_{calt})/Nu_{calt}');
legend({'Nusselt error','Second order comparison', '1st order'}, 'Location', 'southeast');

box on; grid on;
