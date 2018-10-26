%Load some output and calculate the Nusselt number
%clear all;
close all;

set(0, 'defaultlinelinewidth',2);
set(0, 'defaultaxeslinewidth',2);
set(0, 'defaultpatchlinewidth',2);
set(0, 'defaultAxesFontSize',20);
set(0, 'DefaultLineMarkerSize',10);

Nu_caltagirone = containers.Map({0,39,50,100,200,250, 300}, [1 1  1.45 2.651 3.813 4.199 4.523]);

Nu_parkinson_av = [];
Nu_parkinson_base = [];
% grid_res = containers.Map({0,39,50,100,200,250},...
%     [32,32,32,32,32,32]);

subcycled = false;
grid_res = containers.Map({0,39,50,100,200,250,300},...
    [32,32,128,128,128,128, 128]);

error = [];

k = keys(grid_res);
Ra_park = cell2mat(k);

dim = 2;


for Ra_i =1:length(Ra_park)
    %Frame = 7343;
    %plot_prefix = 'plotRa100-';
    Ra = Ra_park(Ra_i);
    
    frame = -1;
    res = grid_res(Ra);
    adaptive = '';
    if Ra > 300
        adaptive = 'adaptive-';
    end
    
    plot_prefix = ['bm2-', adaptive, 'Ra', num2str(Ra), '-',num2str(res),'pts-steady'];
    
    %Moving to dropbox as it syncs with Ubuntu
    
   % output_dir = '/Users/parkinsonj/Dropbox/DPhil/Data/'; %OS X
    output_dir = '/home/parkinsonjl/Dropbox/DPhil/Data/'; %Linux
    
    output = MushyLayerOutput(dim, frame, output_dir, plot_prefix, subcycled);
    
   
    
    [Nu_parkinson_av(Ra_i), Nu_parkinson_base(Ra_i)] = output.nusseltNumber();
    error(Ra_i) = abs(Nu_parkinson_base(Ra_i) - Nu_caltagirone(Ra));
    ratio = Nu_parkinson_base(Ra_i)/Nu_caltagirone(Ra)
end

hFig = figure();
set(hFig, 'Position', [400 400 600 400])

plot( cell2mat(keys(Nu_caltagirone)), cell2mat(values(Nu_caltagirone)), '-ro', ...
    Ra_park, Nu_parkinson_base, '--b+'); %, ...
   % Ra_park, Nu_parkinson_av, '-.g*');
%title('Nusselt number for HRL benchmark');
xlabel('Ra');
ylabel('Nu');
grid on; box on;
[hleg1, hobj1] = legend({' Caltagirone (1975)', ' This study'},'Location','southeast');
set(hleg1, 'Position', [0.5 0.2 0.3 0.15]);


axis([0 300 0 5]);

print('NusseltComparison','-dpng', '-r300');





