close all;

%local_base_dir = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/matlab/EGUPoster/';
%remote_base = '/home/parkinsonjl/mnt/sharedStorage/middleton/' ; %CR1.179RaC80Le200KozenyPermeabilitypts256-S3.5-TB',folders{i},'.0-R0.13-', endings{j},'/diags.mat'];
 
local_base_dir = '/home/parkinsonj/convection-in-sea-ice/MushyLayer/matlab/EGUPoster/';
remote_base = '/network/group/aopp/oceans/AW002_PARKINSON_MUSH/middleton/';


folders = {'-10','-15','-20'};

fprintf('Concat diags \n');

for i=1:length(folders)
    
    concatDiagFile = [local_base_dir, 'T',folders{i}, '/diagConcat.mat'];
    
    concat_time = [];
    concat_Si = [];
    concat_Vi = [];
    concat_Sbox = [];
    concat_frames = [];
    
   
    theseFolders = {[remote_base, 'CR1.179RaC80Le200KozenyPermeabilitypts256-S3.5-TB',folders{i},'.0-R0.13-0/'], ...
        [remote_base, 'CR1.179RaC80Le200KozenyPermeabilitypts256-S3.5-TB',folders{i},'.0-R0.13-restart/'], ...
        [remote_base, 'CR1.179RaC80Le200KozenyPermeabilitypts256-S3.5-TB',folders{i},'.0-R0.13-restart2/'], ...
        [remote_base, 'T',folders{i},'Outflow']};
    
   for j=1:length(theseFolders)
      % endings = {'0','restart','restart2'};
      thisDiagFile = [theseFolders{j}, '/diags.mat'];
      %thisDiagFile  =  [data_dir, folders{i}, '/diags',num2str(j),'.mat'];
      
      clear frames;
      
      
      if exist(thisDiagFile, 'file') == 2
         load(thisDiagFile);
         
         %Print some info about what we've loaded
         fprintf('Loaded times, length %d \n', length(times));
         fprintf('Loaded Si, length %d \n', length(Si));
         fprintf('Loaded Vi, length %d \n', length(Vi));
         fprintf('Loaded Sbox, length %d \n', length(Sbox));
         fprintf('Loaded frames, length %d \n', length(frames));
         
          if exist('frames','var') == 1
              
          else
             frames = NaN*times; 
          end
         
         l = length(times);
          concat_time(end+1:end+l) = times;
    concat_Si(end+1:end+l) = Si;
    concat_Vi(end+1:end+l) = Vi;
    concat_Sbox(end+1:end+l) = Sbox;
    concat_frames(end+1:end+l) = frames;
         
      end
       
   end
   
   figure(); plot(concat_time, concat_Si, 'x'); title(['$T=',folders{i},'$']);
   
   save(concatDiagFile, 'concat_time','concat_Si', 'concat_Vi', 'concat_Sbox', 'concat_frames');
   
  
   
end


fprintf('Concat diags done \n');

