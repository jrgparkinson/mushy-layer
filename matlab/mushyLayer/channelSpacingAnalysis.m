% Plot evolution of channel spacing for wide domain simulations
clear all; close all;


leg = {};

base_folder ='/media/parkinsonjl/FREECOM HDD/channelSpacing-periodic/';
folders = {};

dirs = dir(base_folder);
for i = 1:length(dirs)
    % Ignore stupid entries like '.' and '..'
    if length(dirs(i).name) > 3
        folders{end+1} = dirs(i).name;
    end
end

times = NaN*ones(length(folders), 1);

spacings = NaN*ones(length(folders), 1);


h = figure();
set(h, 'Position', [100 200 1600 1200]);
title('Average channel spacing');
hold on;

%Ras = [200, 400, 800, 1600];
%for Ra_i = 1:length(Ras) 
   %folder = ['/media/parkinsonjl/FREECOM HDD/channelSpacing-periodic/CR1.25RaC', ...
   %    num2str(Ras(Ra_i)), 'Le200ChiCubedPermeabilitypts1024-0'];
for folder_i = 1:length(folders)
   folder = [base_folder, folders{folder_i}];
    
   diags = [folder, '/channelSpacing.out'];
   
   inputs = readInputs([folder, '/inputs']);
   
   Ra = str2double(inputs.rayleighComp);
   
   lineType = '--';
   symbol = 'x';
   if Ra == 800
       %lineType = '-';
       symbol = 'o';
   end
   
   
   fileID = fopen(diags);
   
   if exist(diags, 'file') == 0
       continue;
   end
   
   C = textscan(fileID,'%s | %s');
   
   file_times = C{1};
   file_spacings = C{2};
   
   t = []; spacings = [];
   for j = 2:length(file_times)
      t(end+1) = str2num(file_times{j}); 
      spacings(end+1) = str2num(file_spacings{j});
   end
   
   if sum(isinf(spacings)) ~= length(spacings)
   
       plot(t, spacings, [symbol, lineType]);

       if length(t) > 1
             leg{end+1} = folders{folder_i};
       end
   end

    
end
ylim([0 0.25]);
ylabel('Channel spacing, \lambda');




hold off;
xlim([0.8 2.5]);

lgd = legend(leg, 'Location','eastoutside');
box on;
xlabel('t');
