% This is a collection of matlab scripts for analysing simulations of
% directional solidification at varying domain size

% There are three main functions:
% ** computeSimplifiedDiags(base_dir, recompute)
%    extracts the diagnostics of  interest from the simulations
% ** processEachFolder(base_dir, param_folders)
%    plots the time dependent fluxes for a single set of parameters and
%    various box sizes
% ** summariseFluxes(base_dir, param_folders)
%    plots the 'steady' fluxes for each set of parameters and box sizes


function channelSpacingAnalysis

close all;

base_dir = '/home/parkinsonjl/mnt/sharedStorage/channelSpacingV2/';

param_folders = dir([base_dir, 'CR*']);

global useWellsUnits;
useWellsUnits = true;

%computeSimplifiedDiags(base_dir, true);
computeSimplifiedDiags(base_dir, false);

useWellsUnits = false;
%processEachFolder(base_dir, param_folders)

useWellsUnits = true;
summariseFluxes(base_dir, param_folders);

end

function computeSimplifiedDiags(base_dir, recompute)

if nargin < 2 
    recompute = false;
end

param_folders = dir([base_dir, 'CR*']);

for p_folder_i = 1:length(param_folders)
    param_folder = param_folders(p_folder_i).name;
    this_dir = [base_dir, param_folder, '/'];
    
    data_file = [this_dir, 'data.mat'];
    if exist(data_file, 'file') == 2 && ~recompute
        fprintf('Skipping %s \n', param_folder);
        continue
    end
        
    fprintf('Processing %s \n', param_folder);
      
    
    folders = dir([this_dir, 'pts*']);
   
    folders_to_process = {};
    
    for folder_i=1:length(folders)       
        folder = folders(folder_i).name;
        folders_to_process{end+1} = folder;         
    end    
   
    clear('simplified_diags');
    
    for folder_i=1:length(folders_to_process)
        
        fname = folders_to_process{folder_i}; %strrep(folder, trim_prefix, '');
        fprintf('.. %s \n', fname);
        
        % Get cpp_diagnostics (WellsFlux) too
        
        diags = getDiagnostics([this_dir, fname, '/diagnostics.out']);
        
        wells_diags = getDiagnostics([this_dir, fname, '/cpp_diagnostics.out']);
        
        if isfield(diags, 'Fs_vertical_av')
            
          
            Fs = diags.Fs_vertical_av;
             
            simplified_diags(folder_i).name = fname;
            simplified_diags(folder_i).time = diags.time;
            simplified_diags(folder_i).flux = Fs;
            
        end
        
        if isfield(wells_diags, 'Time')
            
            simplified_diags(folder_i).wells_time = wells_diags.Time;
            simplified_diags(folder_i).wells_flux = -wells_diags.Fs_average;
            
        else
            simplified_diags(folder_i).wells_time = [NaN];
            simplified_diags(folder_i).wells_flux = [NaN];
            
        end
        
        
    end
    
    save(data_file, 'simplified_diags');
    
end

end

function summariseFluxes(base_dir, param_folders)

global useWellsUnits;

fixHeight = -2;
fixedResolution = 128;
fixedWidth = 1.0;

fixHeight = 2;
fixedWidth = -1.0;

h = figure();
h.Position = [50 400 1800 900];

hold on;

p = [];
leg  = {};

for p_folder_i = 1:length(param_folders)
    fname = param_folders(p_folder_i).name;
    
    clear('simplified_diags');
    load([base_dir, fname, '/data.mat']);
    
    disp(fname);
    
    % Compute width and steady flux
    
    width = [];
    height = [];
    steady_flux = [];
    error = [];
    
    simplified_diags = myFolderSort(simplified_diags);
    
    for i=1:length(simplified_diags)
        % Get the width from a string like pts128-width1.0-aspect0.125-0 
        thisName = simplified_diags(i).name;
        
        [tokens, matches] = regexp(thisName,'pts([^-]*)-width([^-]*)-aspect([^-]*)','tokens', 'Match');
        thisPts = str2num(tokens{1}{1});
        thisWidth = str2num(tokens{1}{2});
        thisAspect = str2num(tokens{1}{3});
        
        thisHeight = thisWidth/thisAspect;
        thisResolution = thisPts/thisWidth;
        if fixHeight > 0 && thisHeight ~= fixHeight ...
            || fixedResolution > 0 && thisResolution ~= fixedResolution ...
            || fixedWidth > 0 && thisWidth ~= fixedWidth
            continue;
            
        end
        
        
        width(end+1) = thisWidth;
        height(end+1) = thisHeight;
        %disp(thisName);
        
        
        % Get the steady flux by averaging time series then taking the last value
        
        if useWellsUnits
        Fs = simplified_diags(i).wells_flux;    
        else
        Fs = simplified_diags(i).flux;
        end
        
        [thisSteadyFlux, thisFluxError] = computeSteadyFlux(Fs);
        
        error(end+1) = thisFluxError;
        
        steady_flux(end+1) = thisSteadyFlux;
        
        
    end
    
    if sum(~isnan(steady_flux)) < 2
        continue;
    end
    
    
   if fixHeight > 0
       dep_var_name = 'Width';
       dependent_variable = width;
   
   else
       dep_var_name = 'Height';
       dependent_variable = height;
   end
   
    p(end+1) = errorbar(dependent_variable, steady_flux, error, '-x');
    xlabel(dep_var_name);
    ylabel('Steady flux');
    
    leg{end+1} = fname;
    
    legend(p, leg, 'Location', 'eastoutside');
    box on;
    drawnow;
    
    disp(steady_flux);
end

hold off;

end

function N = getNSmoothSteps(flux_timeseries)

N = 2e4;
end

function [thisSteadyFlux, thisFluxError] = computeSteadyFlux(Fs)
NSmoothSteps = getNSmoothSteps(Fs);

movingMeanFs = movmean(Fs, NSmoothSteps);

% Compute standard deviation over last NSmoothSteps
endi = length(Fs);
starti = max(1, endi-NSmoothSteps);
stdev = std(Fs(starti:endi));

thisFluxError = stdev;

thisSteadyFlux = movingMeanFs(end);
if abs(thisSteadyFlux) > 1e2
    thisSteadyFlux = NaN;
end

end

function processEachFolder(base_dir, param_folders)

global useWellsUnits;

for p_folder_i = 1:length(param_folders)
    param_folder = param_folders(p_folder_i).name;
    
    fprintf('Processing %s \n', param_folder);
    
    %param_folder =  'CR2.000RaC100Le100KozenyPermeabilityDa1.0e-02R1.0e-03';
    this_dir = [base_dir, param_folder, '/'];
    
    
    % Compare all files that start with this:
    
    % Decide what to plot
    fixedHeight = -1;
    fixedWidth = -1;
    resolution = -1;
    
    fixedHeight = 2.0;
    %fixedWidth = 2.0;
    %resolution=128;
    
    clear('simplified_diags');
    load([base_dir, param_folder, '/data.mat']);
    
    sorted_folders = myFolderSort(simplified_diags);
    
    h = figure();
    h.Position = [50 400 1800 900];
    leg = {};
    hold on;
        
    %diags_to_plot = [];
    clear('diags_to_plot');
    
    for folder_i=1:length(sorted_folders)
        
        folder = sorted_folders(folder_i).name;
        
      
        % Get details and check folder conforms to requirements
        folders_regex = 'CR(\d+\.?\d+?)RaC(\d+\.?\d+?)Le(\d+e?\d+)(\w+)pts(\d+)-(\w+)'; %Make sure it ends in 0
        folders_regex = '.*pts(\d+)-width(\d+\.\d+)-aspect(\d+\.\d+)-0';
        [~,matches,~]  = regexp(folder, folders_regex, 'match', 'tokens', 'tokenExtents');
        
        if isempty(matches)
            continue;
        else
            match = matches{1};
            pts = str2double(match{1});
            width = str2double(match{2});
            aspect = str2double(match{3});
            
            height = width/aspect;
            res = pts/width;
            
            if fixedHeight > -1 && height ~= fixedHeight
                continue;
            end
            
            if fixedWidth > -1 && width ~= fixedWidth
                continue;
            end
            
            if resolution > -1 && res ~= resolution
                continue
            end
            
            i = 1;
            
            if exist('diags_to_plot')
            i=length(diags_to_plot)+1;
            end
            
            diags_to_plot(i) = sorted_folders(folder_i);
            
        end
        
    end
    
    colors = distinguishable_colors(length(diags_to_plot));
    
    p = [];
    
    for folder_i=1:length(diags_to_plot)
         diags = diags_to_plot(folder_i);
         
        fname = diags.name; %strrep(folder, trim_prefix, '');
        %fprintf('.. %s \n', fname);
        
       % diags = getDiagnostics([this_dir, fname, '/diagnostics.out']);
      
        
        if isfield(diags, 'flux')
                
            % Try plotting moving mean
            if useWellsUnits
               Fs = diags.wells_flux;
               t = diags.wells_time;
            else
                Fs = diags.flux;
                t = diags.time;
            end
            NSmoothSteps = getNSmoothSteps(Fs); % What should this be??
            movingMeanFs = movmean(Fs, NSmoothSteps);
            
            c = colors(folder_i, :);
            p(folder_i) = plot(t, Fs,           'LineWidth', 1.0, 'Color', c);
            plot(t, movingMeanFs, 'LineWidth', 3.0, 'Color', c);
            
            leg{end+1} = fname;
            %leg{end+1} = sprintf('%s-Smoothed-%d', fname, NSmoothSteps);
            
            %theseDiags = struct('time', diags.time, 'flux', Fs);
            % simplified_diags(end+1) = theseDiags;
            simplified_diags(folder_i).name = fname;
            simplified_diags(folder_i).time = diags.time;
            simplified_diags(folder_i).flux = Fs;
            
            [thisSteadyFlux, thisFluxError] = computeSteadyFlux(Fs);
            
            % Make an extended set of times
            extended_time = diags.time;
            if length(extended_time) > 10
                dt = extended_time(end) - extended_time(end-1);
                for i=1:NSmoothSteps
                    extended_time(end+1) = extended_time(end) + dt;
                end
            end
            
            plot(extended_time, extended_time*0 + thisSteadyFlux, 'LineWidth', 1.0, 'LineStyle', '--', 'Color', c);
            %plot(diags.time, movingMeanFs, 'LineWidth', 3.0, 'Color', c);
            
        end
        
        legend(p, leg, 'Location', 'eastoutside');
        %xlim([0 1e5]);
        %  ylim([-2 -1]);
        xlabel('$t$');
        ylabel('$F_s$');
        title({param_folder, ''});
        box on;
        drawnow;
    end
    
    
    hold off;
    
    
    % Save image
    save_loc = [this_dir, 'fluxes.png'];
    
    print(h,save_loc,'-dpng','-r150')
    
    save_loc_2 = [base_dir, '/', param_folder, '.png'];
    print(h,save_loc_2,'-dpng','-r150')
    
%    save([this_dir, 'data.mat'], 'simplified_diags');
    
end

end



function sortedFolders = myFolderSort(folders)
fnames = {};

for i=1:length(folders)
    fnames{end+1} = folders(i).name;
end

q0 = regexp(fnames,'\d*','match');
q1 = str2double(cat(1,q0{:}));
[~,ii] = sortrows(q1,[2 1]);
out = fnames(ii);

%disp(out);

sortedFolders = folders(ii);

end

