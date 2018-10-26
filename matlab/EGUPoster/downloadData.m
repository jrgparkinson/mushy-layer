function downloadData(local_dir, folder, remote_base)

endings = {'0', 'restart', 'restart2'};
%remote_base = '/home/parkinsonjl/mnt/sharedStorage/middleton/';

T = folder(2:end);

remote_dirs = {};
for i =1:length(endings)  
   remote_dirs{i} = [remote_base, 'CR1.179RaC80Le200KozenyPermeabilitypts256-S3.5-TB',T,'.0-R0.13-',endings{i},'/'];
end

remote_dirs{end+1} = [remote_base, folder, 'Outflow/'];


for i =1:length(remote_dirs)  
   remote_dir = remote_dirs{i};
    
   files = dir([remote_dir, 'toPlot*.mat']);
   for j=1:length(files)
       name = files(j).name;
       
       frame = strrep(name, 'toPlot', '');
       frame = str2num(strrep(frame, '.mat',''));
       
       if length(frame) == 0
           continue
       end
       
       localFileExists = (exist([local_dir, name], 'file') == 2);
       getFile = (mod(frame, 200) == 0);
       if  getFile && ~localFileExists
           
          fprintf('Copying %s \n', name);
          copyfile([remote_dir, name], local_dir);
          
       end
   end
end


end