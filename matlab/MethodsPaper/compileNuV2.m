% Compile all

function compileNuV2(base_dir, chi, Da, Ra, res_uniform,  res_vm, NuLeBars)

if nargin == 0
    
    
    
   
    % u del u: 0, advection: 1, fixed code (hopefully) - good, clean up
    % some runs
    base_dir = '/home/parkinsonjl/mnt/sharedStorage/TestDiffusiveTimescale/zzConvectionDB-cfl0.23/';
    
    % New
    base_dir = '/home/parkinsonjl/mnt/sharedStorage/TestDiffusiveTimescale/ConvectionDB-cfl0.2/';
    
    
    % Also advection 1, but smaller cfl - quite good, needs finishing some
    % runs
   % base_dir = '/home/parkinsonjl/mnt/sharedStorage/TestDiffusiveTimescale/ConvectionDB-cfl0.09/';
    
   
    %base_dir = '/home/parkinsonjl/mnt/sharedStorage/TestDiffusiveTimescale/ConvectionDB-cfl0.1/';
    
   
    %  u del u: 0, advection: 0, Also OK but needs some work
   %base_dir = '/home/parkinsonjl/mnt/sharedStorage/TestDiffusiveTimescale/ConvectionDB-cfl0.24/';
   
   % Good too, but maybe old
   %base_dir = '/home/parkinsonjl/mnt/sharedStorage/TestDiffusiveTimescale/ConvectionDB-cfl0.15/';
    
    %  chi = '0.4';
    
    %  Da = '1e-06'; % '1.0e-02'
    %  res_uniform = 128;
    %  res_vm = 128;
    
    %  if strcmp(Da, '1e-06')
    %      Ra = {'1.0e+07','1.0e+08','1.0e+09'};
    %      NuLeBars = [1.08,3.07,12.9];
    %  else
    %      Ra = {'1.0e+03','1.0e+04','1.0e+05','5.0e+05'};
    %      NuLeBars =  [1.01,1.41,3.17,5.24];
    % end
    
    
end

UniformPrefix = 'Uniform-convectionDB-';
AMRPrefix = 'VariableMesh2SubcycleRefluxFreestream0.95-ref2-convectionDB-';

prefix = AMRPrefix;



%base_dir = getDataDir(['AMRConvergenceTest/ConvectionDB/']);

folders = dir(base_dir);

% Compile all data
field = '%10s';
fprintf(['%3s | %6s | %7s |',field,'|',field,'|',field,'|',field,'\n'], '\chi', 'Da', 'Ra', ...
   'Nu (32)', 'Nu (64)', 'Nu (128)', 'Nu (Le Bars)');
Nxs = [32, 64, 128];



for f_i=1:length(folders)
    fname = folders(f_i).name;
    
    % Get parameters for this folder
    scientific = '\d+\.\d+e[+-]\d+'; %'1.0e-06';
    
    regexStr = ['chi(\d\.\d)-Da(',scientific,')-Ra(',scientific,')'];
    [tokens,matches] = regexp(fname,regexStr,'tokens','match');
    
    %celldisp(tokens);
    
    if length(tokens) > 0
        chi = str2num(tokens{1}{1});
        Da = str2num(tokens{1}{2});
        Ra = str2num(tokens{1}{3});
        
        if chi == 1.0
            continue
        end

        fprintf('%1.1f | %1.1e | %1.1e |', chi, Da, Ra);

         % Now get Nu for different resolutions
         subfolders = dir(fullfile(base_dir, fname, [prefix, '*']));
         Nu_vals = num2cell(NaN*Nxs);
         for s_i = 1:length(subfolders)
             sname = subfolders(s_i).name;
             
             fullfolder = fullfile(base_dir, fname, sname);
             
             [tokens,matches] = regexp(sname,[prefix, '(\d+)--0'],'tokens','match');
             if length(tokens) > 0
                Nx = str2num(tokens{1}{1});
                
                [Nu_subfolder, stillRunning] = getNuFolder(fullfolder);
                
                Nu_i = find(Nxs == Nx);
                if isempty(Nu_i)
                    % Do nothing
                elseif stillRunning
                    Nu_vals{Nu_i} = sprintf('%8.2f*', Nu_subfolder);  
                else
                    Nu_vals{Nu_i} = sprintf('%9.2f', Nu_subfolder);
                end
             
             end
         end
         
         for Nu_i = 1:length(Nu_vals)
             fprintf('%10s|', Nu_vals{Nu_i}); 
         end
         
         
         % Print Nu (Le Bars)
         format = '%9.2f|';
         nu = getLeBarsNu(chi, Da, Ra);
         fprintf(format, nu);
         
         
       % fprintf(' ... | ... | ... ');

        fprintf('\n');
    end
end

end



function [Nu, stillRunning] = getNuFolder(folder)

%fprintf('Processing subfolder %s \n', folder);

Nu = NaN;
stillRunning = false;

diagFileExists = false;
diagFile = [folder,'/diagnosticsLatest.out'];
if exist(diagFile, 'file') == 2
    diagFileExists = true;
else
    diagFile = [folder,'/diagnostics.out'];
    if exist(diagFile, 'file') == 2
        diagFileExists = true;
    end
end

if diagFileExists
    
    diags = getDiagnostics(diagFile);
    
    %     recentMean = mean(diags.Nusselt(round(end-end/10):end));
    %     recentStd = std(diags.Nusselt(round(end-end/10):end));
    %
    %     recentDiff = abs(diags.Nusselt(end) - diags.Nusselt(end-1));
    %
    %     err = max(recentDiff, recentStd);
    
    %fprintf('Ra=%s, Nu=%1.5f +/- %1.5f \n', Ra, diags.Nusselt(end), err);
    
    if isfield(diags, 'Nusselt')
    
        Nu = diags.Nusselt(end);
    else
        Nu = NaN;
        
    end
    
        FileInfo = dir(diagFile);
        TimeStamp = FileInfo.datenum;

        ageDays = now - TimeStamp; % in days
        ageMinutes = ageDays*24*60;
        ageSeconds = ageDays*24*3600; % in seconds


        if ageMinutes > 60
            stillRunning = false;
        else
            stillRunning = true;
        end
        
   
    
    
end

% Some error checking
if Nu < 0 || Nu > 100
    Nu = NaN;
end
end

function folder = getOutputFolder(Ra, chi, Da, cfl)
if nargin < 4
    cfl = -1;
end

if cfl < 0
    folder = ['chi',chi,'-Da',Da,'-Ra',Ra];
else
    folder = ['chi',chi,'-Da',Da,'-Ra',Ra, '-cfl',cfl];
end

end

