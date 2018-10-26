% A variety of tests for Rayleigh-Benard convection following Toppaladoddi
% et al. (2015)
clear all; close all;

set(groot, 'DefaultTextInterpreter', 'latex')
set(groot, 'DefaultLegendInterpreter', 'latex')


doTest1 = true;
doTest2 = false;

base_dir = '/home/parkinsonjl/mnt/sharedStorage/AW002_PARKINSON_MUSH/AMRConvergenceTestRayleighBenard/';

% 1) Horizontally averaged temperated vs y (vertical co ordinate)
% for Ra = 4000, Pr = 0.71
%output_dir = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/output-RB-fixed/';
%file = 'plt003431';

if (doTest1)
    
    grid_res = [8, 16, 32, 64, 128];
    legend_str = {'8x48', '16x96', '32x192', '64x384', '128x768'};
    frames =    {'010096', '009600', '000993', '000600', '000200'};
    
    
    grid_res = [8, 16, 32];
    legend_str = {'Lipps (1976)', '8x48', '16x96', '32x192'};
    frames =    {'000714', '000769', '038000'};
    lineStyles = {':x', '--d','-.+'};

frame = -1;
dim = 2;
subcycled = true;

%y_max = y(end) + (y(end) - y(end-1))/2;
%y = y/y_max;

img = imread('toppaladoddi2.png');

%axis limits, same for T and y
Nu_min = 0;
max = 1;

hFig = figure();
set(hFig, 'Position', [300 400 1200 800])

%imagesc([Nu_min max], [Nu_min max], flipud(img));
lippsData = csvread('lipps.data', 1);

y_lipps = lippsData(:, 2);
T_lipps = lippsData(:, 1)+0.5;
plot(T_lipps, y_lipps , 's');

hold on;
%plot(T_av, y, 'magenta');


for res_i=1:length(grid_res)
    res = grid_res(res_i);
    frame = frames{res_i};
    output_dir = [base_dir, 'Uniform-rayleighBenard-', num2str(res), '-0/'];
file = ['AMRConvergenceTestRayleighBenard-Uniform-rayleighBenard-', num2str(res), '-', frame];

mlTest = MushyLayerOutput(dim, -1, output_dir, file, subcycled);


T = mlTest.dataForComp(mlTest.components.Temperature) - 1; % Have to subtract one

T_av = mean(T, 1);
[X,Y] = mlTest.grid();
y = Y(:,1);

%Do averaging of coarse cells
T_new = [];
y_new = [];
skip = false;
for i=1:(length(T_av)-1)
    if skip
        % do nothing
        skip = false;
    else
   if T_av(i) == T_av(i+1)
       T_new(end+1) = (T_av(i)+T_av(i+1))/2;
       y_new(end+1) = (y(i)+y(i+1))/2;
       skip = true;
   else
       T_new(end+1) = T_av(i);
       y_new(end+1) = y(i);
   end
   
    end
    
end

% Get the last data point
if skip == false
    T_new(end+1) = T_av(i+1);
    y_new(end+1) = y(i+1);
end

plot(T_new, y_new, lineStyles{res_i});

end

legend(legend_str, 'Location', 'eastoutside');


set(gca, 'ydir', 'normal');

ylabel('$y$');
xlabel('$< T >$');
%title('NB: toppaladoddi axes wrong way around?');

box on; grid on;

hold off;

end

if (doTest2)
% 2) Nu(Ra) scaling
%output_dir = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/output-NuRa/';
%file = 'plt';
output_dir = '/home/parkinsonjl/convection-in-sea-ice/test/';

frame = -1;
dim = 2;
subcycled = true;


Ra =         [2000,   2500,  3000,  5000,  10000 20000 30000 50000];
resolution = [8,      8,     8,     16,     16,    16,   32,   32];
Nu_clever =  [1.212, 1.475, 1.663,  2.116, 2.661 3.258 3.662 4.245] ; % 3.258, 3.662, 4.245];
Rac = 1708;
Nu_scaling = 1.56*(Ra/Rac).^0.296;

Nu = NaN * Ra;
for Ra_i = 1:length(Ra)
    file = ['rayleighBenard-Ra', num2str(Ra(Ra_i)), '/rayleighBenard-Ra', num2str(Ra(Ra_i)),'-',num2str(resolution(Ra_i)),'pts-steady' ];
    mlTest = MushyLayerOutput(dim, frame, output_dir, file, subcycled);
    
    % 1 means vertical heat transfer (dir = 1)
    [Nu_av, Nu_base] = mlTest.nusseltNumber(1);
    
    Nu(Ra_i) = Nu_base;
end


%axis limits, same for T and y
Nu_min = 0.8;
Nu_max = 6;

Ra_min = 10^3;
Ra_max = 10^5;

hfig = figure();
set(hfig, 'Position', [300 300 1000 600]);


semilogx(Ra, Nu, 'x', Ra, Nu_clever, 'o', Ra, Nu_scaling, '--');
axis([Ra_min, Ra_max, Nu_min, Nu_max]);
%plot(Ra, Nu, 'x');

legend('My results', 'Clever and Busse (1974)', 'Nu=1.56(Ra/Ra_c)^{0.296}', 'Location', 'southeast');

xlabel('Ra');
ylabel('Nu');
box on; grid on; 
%hold off;

percent_err = 100*(Nu-Nu_clever)./Nu_clever

end