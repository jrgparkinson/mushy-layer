% Compare AMR and uniform grid results for Rayleigh Benard convection
clear all; close all;


% 1) Horizontally averaged temperated vs y (vertical co ordinate)
% for Ra = 4000, Pr = 0.71
%output_dir = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/output-RB-fixed/';
%file = 'plt003431';


    
output_dir = '/home/parkinsonjl/convection-in-sea-ice/test/rayleighBenard-Ra4000/';
output_dir_amr = '/home/parkinsonjl/convection-in-sea-ice/test/rayleighBenard-adaptive-Ra4000/';
coarse_file = 'rayleighBenard-Ra4000-8pts-steady';
fine_file = 'rayleighBenard-Ra4000-16pts-steady';
amr_file = 'rayleighBenard-adaptive-Ra4000-8pts-steady';
%output_dir_amr = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/output-adaptive-full/';
%file = 'plt000800';

frame = -1;
dim = 2;
subcycled = true;

interpMethod = 'bicubic';
%interpMethod = 'nearest';

mlCoarse = MushyLayerOutput(dim, frame, output_dir, coarse_file, subcycled);
mlFine = MushyLayerOutput(dim, frame, output_dir, fine_file, subcycled);
mlAMR = MushyLayerOutput(dim, frame, output_dir_amr, amr_file, subcycled, interpMethod);

T_coarse = mlCoarse.dataForComp(mlCoarse.mlComps.temperature);
T_fine = mlFine.dataForComp(mlCoarse.mlComps.temperature);
T_amr = mlAMR.dataForComp(mlAMR.mlComps.temperature);


T_coarse_interpolated = resizem(T_coarse, 2, interpMethod);

T_err_coarse = T_coarse_interpolated - T_amr;
T_err_fine = T_fine - T_amr;

T_err_coarse_fine = T_coarse_interpolated - T_fine;

[X,Y] = mlAMR.grid();
y = Y(:,1);
x = X(1, :);

hfig = figure();
set(hfig, 'Position', [200 200 1400 1200]);

subplot(4, 1, 1);

imagesc(x, y, T_amr.');
colorbar();
set(gca,'YDir','normal')
axis normal;
title('AMR');

subplot(4, 1, 2);
imagesc(x, y, T_err_fine.');
colormap(bluewhitered), colorbar;
set(gca,'YDir','normal')
title('AMR - fine');
caxis([-0.1 0.1]);


subplot(4, 1, 3);
imagesc(x, y, T_err_coarse.');
colorbar();
set(gca,'YDir','normal')
title('AMR - coarse');
caxis([-0.1 0.1]);

subplot(4, 1, 4);
imagesc(x, y, T_err_coarse_fine.');
colorbar();
set(gca,'YDir','normal')
title('coarse - fine');
caxis([-0.1 0.1]);



T_uniform = mean(T_coarse, 1);
T_AMR = mean(T_amr, 1);




%Do averaging of coarse cells
% T_new = [];
% y_new = [];
% skip = false;
% for i=1:(length(T_uniform)-1)
%     if skip
%         % do nothing
%         skip = false;
%     else
%    if T_uniform(i) == T_uniform(i+1)
%        T_new(end+1) = (T_uniform(i)+T_uniform(i+1))/2;
%        y_new(end+1) = (y(i)+y(i+1))/2;
%        skip = true;
%    else
%        T_new(end+1) = T_uniform(i);
%        y_new(end+1) = y(i);
%    end
%    
%     end
%     
% end

% Get the last data point
% if skip == false
%     T_new(end+1) = T_uniform(i+1);
%     y_new(end+1) = y(i+1);
% end

%y_max = y(end) + (y(end) - y(end-1))/2;
%y = y/y_max;

% img = imread('toppaladoddi2.png');

% %axis limits, same for T and y
% Nu_min = 0;
% max = 1;
% 
% figure();
% %imagesc([Nu_min max], [Nu_min max], flipud(img));
% hold on;
% %plot(T_av, y, 'magenta');
% plot(T_uniform, y);
% plot(T_AMR, y);
% 
% set(gca, 'ydir', 'normal');
% 
% ylabel('y');
% xlabel('<T>');
% 
% box on; grid on;
% 
% hold off;
% 
% legend('Uniform', 'AMR');
% 
