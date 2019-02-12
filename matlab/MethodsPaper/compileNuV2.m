% Compile all

function compileNuV2(base_dir)

if nargin == 0
    % 0.18 has higher order bcs
    % 0.17 has higher order advection
    base_dir = '/home/parkinsonjl/mnt/sharedStorage/TestDiffusiveTimescale/ConvectionDB-cfl0.17/';
end

UniformPrefix = 'Uniform-convectionDB-';
AMRPrefix = 'VariableMesh2SubcycleRefluxFreestream0.95-ref2-convectionDB-';

prefix = UniformPrefix;

close all;
Nu_field_name = 'Nusselt';
%Nu_field_name = 'NusseltMiddle';
%Nu_field_name = 'NusseltLeft';

%base_dir = getDataDir(['AMRConvergenceTest/ConvectionDB/']);

folders = dir(base_dir);

% Compile all data
Nxs = [64, 128, 256, 512, 1024, 2048];

field = '%10s';

fields_list = [ field,'|'];
column_headers = {'\chi', 'Da', 'Ra'};
for i=1:length(Nxs)
    fields_list = [fields_list, field,'|'];
    column_headers{end+1} = ['Nu (', num2str(Nxs(i)) ,')'];
end
column_headers{end+1} = 'Nu (Le Bars)';

%fprintf(['%3s | %6s | %7s |',field,'|',field,'|',field,'|',field,'\n'], '\chi', 'Da', 'Ra', ...
%   'Nu (32)', 'Nu (64)', 'Nu (128)', 'Nu (Le Bars)');
fprintf(['%3s | %6s | %7s |',fields_list,'\n'], column_headers{:});


m = 3;
n= 4;
 figure('Position', [100 100 1600 1000]);
 subplot_i = 0;
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
         Nu_floats = NaN*ones(1, length(Nxs));
         for s_i = 1:length(subfolders)
             sname = subfolders(s_i).name;
             
             fullfolder = fullfile(base_dir, fname, sname);
             
             [tokens,matches] = regexp(sname,[prefix, '(\d+)--0'],'tokens','match');
             if length(tokens) > 0
                Nx = str2num(tokens{1}{1});
                
                [Nu_subfolder, stillRunning] = getNuFolder(fullfolder, Nu_field_name);
                
                Nu_i = find(Nxs == Nx);
                if isempty(Nu_i)
                    % Do nothing
                elseif stillRunning
                    Nu_vals{Nu_i} = sprintf('%8.2f*', Nu_subfolder);  
                else
                    Nu_vals{Nu_i} = sprintf('%9.2f', Nu_subfolder);
                end
                
                Nu_floats(Nu_i) = Nu_subfolder;
             
             end
         end
         
         for Nu_i = 1:length(Nu_vals)
             fprintf('%10s|', Nu_vals{Nu_i}); 
         end
         
         if f_i > 3
         
        subplot_i =subplot_i + 1;
         subplot(m, n, subplot_i);
         plot(1./Nxs, Nu_floats, 'x-');
         xlabel('$\Delta x$');
         ylabel('Nu');
         title(sprintf('Da=%.0e, Ra=%.1e', Da, Ra));
         
         %ax = gca;
         %ax.XScale = 'log';
         
         
         subplot_i =subplot_i + 1;
         subplot(m, n, subplot_i);
         plot((1./Nxs).^2, Nu_floats, 'x-');
         xlabel('$\Delta x^2$');
         ylabel('Nu');
         title(sprintf('Da=%.1e, Ra=%.1e', Da, Ra));
         
          %ax = gca;
         %ax.XScale = 'log';
         
         end
        %pause;
         
         
         % Print Nu (Le Bars)
         format = '%9.2f|';
         nu = getLeBarsNu(chi, Da, Ra);
         fprintf(format, nu);
         
         
       % fprintf(' ... | ... | ... ');

        fprintf('\n');
    end
end

end



function [Nu, stillRunning] = getNuFolder(folder, Nu_field_name)

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
    
    if isfield(diags, Nu_field_name)
        Nu = diags.(Nu_field_name)(end);
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

