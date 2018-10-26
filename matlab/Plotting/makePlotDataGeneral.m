function makePlotDataGeneral(subfolders, base_dir,  plot_prefix)

fprintf('Loaded makePlotData\n');

if nargin < 3
    plot_prefix = '';
end

if nargin < 2
    base_dir = getDataDir('middleton/');
else
    base_dir = getDataDir(base_dir);
end

if nargin < 3
    
    % base_dir =  getDataDir('middleton/'); %'/media/parkinsonjl/FREECOM HDD/';
    
    ending3 = '0';
    ending2 = 'restart2';
    ending1 = 'restart';
    
    subfolders = {['CR1.179RaC80Le200KozenyPermeabilitypts256-S3.5-TB-10.0-R0.13-',ending1,'/'], ...
        [ 'CR1.179RaC80Le200KozenyPermeabilitypts256-S3.5-TB-15.0-R0.13-',ending1,'/'], ...
        [ 'CR1.179RaC80Le200KozenyPermeabilitypts256-S3.5-TB-20.0-R0.13-',ending1,'/'], ...
        [ 'CR1.179RaC80Le200KozenyPermeabilitypts256-S3.5-TB-10.0-R0.13-',ending2,'/'], ...
        [ 'CR1.179RaC80Le200KozenyPermeabilitypts256-S3.5-TB-15.0-R0.13-',ending2,'/'], ...
        [ 'CR1.179RaC80Le200KozenyPermeabilitypts256-S3.5-TB-20.0-R0.13-',ending2,'/'], ...
        [ 'CR1.179RaC80Le200KozenyPermeabilitypts256-S3.5-TB-10.0-R0.13-',ending3,'/'], ...
        [ 'CR1.179RaC80Le200KozenyPermeabilitypts256-S3.5-TB-15.0-R0.13-',ending3,'/'], ...
        [ 'CR1.179RaC80Le200KozenyPermeabilitypts256-S3.5-TB-20.0-R0.13-',ending3,'/'], ...
        [ 'T-10Outflow/'], ...
        [ 'T-15Outflow/'], ...
        [ 'T-20Outflow/']         };
    
end

%for i=1:length(subfolders)
%    full_folders{i} = fullfile(base_dir, subfolders{i});
%end
full_folders = fullfile(base_dir, subfolders);


% Three folders to process
for axi = 1:length(full_folders)
    
    % frame = frames(f_i);
    %frame = frames(axi);
    output_dir = full_folders{axi};
    
    fprintf('Output dir: %s\n', output_dir);
    
    
    
    
    
    % If diag file exists, check it hasn't been edited recently (this is the sign of a folder already being processed)
    lockFile = fullfile(output_dir, '/lock');
    lockFileExists = (exist(lockFile, 'file')==2);
    
    if lockFileExists
        fprintf('Directory locked, skipping\n');
        continue
    end
    
    try
        makePlotDataFolder(output_dir, plot_prefix);
    catch
        fprintf('Error, skipping\n');
    end
    
    
    
end % end loop over folders

end


