% Object to analyse a pout file
% Input:
% a) either the location of a single pout file, or,
% b) for runs which have been stopped and restarted, a cell array
% containing the locations of multiple Pout files in the order in which
% they were run

classdef Pout
    properties
        file_loc
        times
        timesteps
        finest_dx
        dSdt
        dTdt
        dUdt
        minimumSl
        Nus
        fluxBottom
        fluxTop
        fluxSponge
        fails
        steadyState
        averageSalinity
        timescale
        dx
        pointsUpdated
        finalHeatFluxMismatch
    end
    
    methods
        
        function obj = Pout(files)
            obj.steadyState = false;
            obj.times = [];
            
            if nargin < 1
                return;
            end
            
            obj.file_loc = files;
            
            if ~iscell(files)
                files = {files};
            end
            
            for file_i = 1:length(files)
                file = files{file_i};
                
                %Options
                
                timestepInterval = 1;
                minTStep = 1;
                
                Ra = 100;
                Nu_caltagirone = containers.Map({0,39,50,100,200,250, 300}, [1 1  1.45 2.651 3.813 4.199 4.523]);
                
                % For Pr = 0.71, wavenumber = 3.117
                %For Ra = 4000 this is just an educated guess
                Ra = 4000;
                Nu_cleverBusse = containers.Map({2000, 2500, 3000, 4000, 5000, 10000, 20000, 30000, 50000}, ...
                    [1.212 1.475 1.663 1.985 2.116 2.661 3.258 3.662 4.245]);
                
                %Nu_analytic = Nu_caltagirone(Ra);
                Nu_analytic = Nu_cleverBusse(Ra);
                
                fileID = fopen(file,'r');
                
                test = exist(file, 'file');
                
                if exist(file, 'file') == 2
                    %fprintf(['Loading', file, '\n']);
                    file_contents = fileread(file);
                else
                    return;
                end
                
                % First need to determine if diagnostics are in pout or if
                % we've put them in a separate diagnostics.out file
                
                diag_file_contents = '';
                
                diagFile = strrep(file, 'pout.0', 'diagnostics.out');
                if exist(diagFile, 'file') == 2
                    %fprintf(['Also loading', diagFile, '\n']);
                    diag_file_contents = fileread(diagFile);
                    %fprintf('... done\n');
                end
                
                
                realNumberRegex = '-?\d+\.?\d*e?[+-]?\d*';
                
                
                
                pattern3 = ['Nusselt number \(level 0\) = (', realNumberRegex,')'];
                [mat3,tok3,ext3]  = regexp(file_contents, pattern3, 'match', ...
                    'tokens', 'tokenExtents');
                
                %          'AMRLevelMushyLayer::postTimeStep 0[\n]',...
                %'Solute flux \(level 0\) = ',realNumberRegex, '[\n]',...
                containsFs = false;
                if (strfind(file_contents, 'Solute flux'))
                    containsFs = true;
                end
                
                if (strfind(file_contents, 'all amr levels said steady state was reached'))
                    obj.steadyState = true;
                end
                
                containsFh = false;
                if (strfind(file_contents, 'Heat flux'))
                    containsFh = true;
                end
                
                containsNu = false;
                if (strfind(file_contents, 'Nusselt number'))
                    containsNu = true;
                end
                
                pattern4 = ['coarse time step\s+(\d+)\s+old time = (',realNumberRegex ,')[^\n]+[\n]'...
                    'd\/dt \(Enthalpy\s*\) = (',realNumberRegex,')\s*[\n]' ...
                    'd\/dt \(Bulk concentration\s*\) = (',realNumberRegex,')\s*[\n]' ...
                    'd\/dt \(Darcy velocity\s*\) = (',realNumberRegex,')\s*[\n]'];
                
                %Remove this from pout - not currently interested in it
                optionalSyncLine = 'Sum\(RHS\) for sync solve = [^\n]+[\n]';
                file_contents = regexprep(file_contents, optionalSyncLine, '');
                
                
                
                % Search for diagnostics in pout.0
                if length(diag_file_contents) == 0
                    
                    
                    %pattern4 = [optionalSyncLine, pattern4];
                    if containsFh
                        pattern4 = ['Heat flux: top = ',realNumberRegex, ', bottom = ',realNumberRegex, '[\n]',...
                            'Relative heat flux mismatch = ',realNumberRegex, '\(flux diff = ',realNumberRegex, ', H diff = ',realNumberRegex, '\)[\n]',...
                            pattern4];
                        
                    end
                    
                    if containsFs
                        pattern4 = ['Solute flux: top = (',realNumberRegex, '), bottom = (',realNumberRegex, ')[\n]',...
                            'Solute flux: left = ',realNumberRegex, ', right = ',realNumberRegex, '[\n]',...
                            'Relative solute flux mismatch = ',realNumberRegex, ' \(flux diff = ',realNumberRegex, ', salt diff = ',realNumberRegex, ' \)[\n]',...
                            pattern4];
                        
                    end
                    
                    
                else
                    % Search for diags in diagnostics.out file
                    %fprintf('Regex search diagnostics file');
                    split = '=====================';
                    bits = strsplit(diag_file_contents, split);
                    
                    %fprintf('Timesteps: %d \n', length(bits));
                    
                    diagPattern = ['Time = (',realNumberRegex, ')[\n]',...
                        'Solute flux: top = (',realNumberRegex, '), bottom = (',realNumberRegex, ')[\n]',...
                        'Solute flux: left = (',realNumberRegex, '), right = (',realNumberRegex, ')[\n]',...
                        'Heat flux: top = (',realNumberRegex, '), bottom = (',realNumberRegex, ')[\n]'];
                        %'Horizontally averaged salinity at \(0, 20\%, 40\%, 60\%\):  (',realNumberRegex, '), (',realNumberRegex, '), (',realNumberRegex, '), (',realNumberRegex, '),\s+[\n]',...
                       % 'Mushy layer depth = (',realNumberRegex, ')'];
                    
                    diagMatches = {};
                    for i=1:length(bits)
                        [~,thisMatch,~]  = regexp(bits{i}, diagPattern, 'match', ...
                         'tokens', 'tokenExtents');
                     %fprintf('%d, ', i);
                     if length(thisMatch) > 0
                        diagMatches{end+1} = thisMatch{1};
                     end
                    end
                end
                
                [mat4,tok4,ext4]  = regexp(file_contents, pattern4, 'match', ...
                    'tokens', 'tokenExtents');
                
                pattern5 =  'Solver failed -- restarting((?!coarse).)*[\n]coarse time step\s+(\d+)'; %[.\n]*AMR([.]*)\n
                [mat5,tok5,ext5]  = regexp(file_contents, pattern5, 'match', ...
                    'tokens', 'tokenExtents');
                
                [~, timescaleMatch, ~] = regexp(file_contents, ['timescale: (',realNumberRegex,')'], 'match', ...
                    'tokens', 'tokenExtents');
                if length(timescaleMatch) > 0
                obj.timescale = str2double(timescaleMatch{1}{1});
                end
                
                numTimesteps = length(tok4);
                
                dTdt = NaN*ones(numTimesteps, 1);
                timestep = 1;
                
                obj.averageSalinity = NaN*ones(numTimesteps, 4);
                
                
                fails = [];
                
                %Get dx:
                [~, dxMatch, ~] = regexp(file_contents, ['read dx = (',realNumberRegex,')'], 'match', ...
                    'tokens', 'tokenExtents');
                if length(dxMatch) > 0
                   dx = str2num(dxMatch{1}{1});
                   obj.dx = dx;
                end
                
                % Get num points updated
                [~, NPointsMatch, ~] = regexp(file_contents, ['total number of points updated = (\d+)'], 'match', ...
                    'tokens', 'tokenExtents');
                if length(NPointsMatch) > 0
                   temp = str2num(NPointsMatch{1}{1});
                   obj.pointsUpdated = temp;
                end
                
                 % Get final heat flux mismatch
                [~, fluxMismatch, ~] = regexp(file_contents, ['Heat flux mismatch = (',realNumberRegex,') \(abs\)'], 'match', ...
                    'tokens', 'tokenExtents');
                if length(fluxMismatch) > 0
                   temp = str2num(fluxMismatch{end}{1});
                   obj.finalHeatFluxMismatch = temp;
                end
                
                
                
                %                 obj.timesteps = NaN*ones(numTimesteps, 1);
                %                 obj.times = [];
                %
                %                 obj.Nus = NaN*ones(numTimesteps, 1);
                %                 obj.fluxBottom  = []; obj.fluxTop = [];
                %                 obj.dSdt =  NaN*ones(numTimesteps, 1); obj.dUdt =  NaN*ones(numTimesteps, 1);
                %                 obj.minimumSl = NaN*ones(numTimesteps);
                
                t_i = length(obj.times) + 1;
                
                for i = 1:length(tok4)
                    
                    thisTok = tok4{i};
                    
                   
                    
                    ind = 1;
                    
                    if length(diag_file_contents) == 0
                        
                        if containsFs
                            
                            obj.fluxTop(t_i) = str2double( tok4{i}{ind});
                            obj.fluxBottom(t_i) = str2double( tok4{i}{ind+1});
                            %obj.fluxSponge(t_i) = str2double( tok4{i}{ind+2});
                            ind = ind + 2;
                            
                        end
                        
                    elseif length(diagMatches) >= i
                        diagTok = diagMatches{i};
                         
                        obj.fluxTop(t_i) = str2double( diagTok{2});
                        obj.fluxBottom(t_i) = str2double( diagTok{3});
                        
                        if length(diagTok) > 7
                        temp =  [str2double(diagTok{8}) str2double(diagTok{9}) ...
                            str2double(diagTok{10}) str2double(diagTok{11})];
                        obj.averageSalinity(t_i, :) = temp;
                        end
                    end
                    
                    
                    
                    tstep =  str2double(tok4{i}{ind});
                    if mod(tstep, timestepInterval) ~= 0 || tstep < minTStep
                        continue;
                    end
                    
                    obj.timesteps(t_i) = tstep;
                    obj.times(t_i) =str2double( tok4{i}{ind+1});
                    obj.dTdt(t_i) = str2double( tok4{i}{ind+2});
                    obj.dSdt(t_i) = str2double( tok4{i}{ind+3});
                    obj.dUdt(t_i) = str2double( tok4{i}{ind+4});
                    
                    if (length(tok3) >= i)
                        obj.Nus(t_i) =  str2double(tok3{i}{1});
                    else
                        obj.Nus(t_i) = NaN;
                    end
                    
                    Nu_err = 100*(obj.Nus(t_i) - Nu_analytic)./Nu_analytic;
                    
                    %fprintf('timestep %d, rate = %s, Nusselt = %s, Nusselt error (percent) = %s \n', ...
                    %    obj.timesteps(t_i), obj.dTdt(t_i), obj.Nus(t_i), Nu_err);
                    
                    t_i = t_i + 1;
                end
                
                
                
                %identify solver fails
                for i = 1:length(tok5)
                    obj.fails(i) = str2double( tok5{i}{2});
                    %fprintf('failed timestep: %d \n', fails(i));
                end
                
            end
            
        end
        
        
        function doPlots(obj)
            movingAverageFlux = [];
            averageLength = 200;
            timeSpan = 0.1;
            movingAverageTimesteps = []; movingAverageTimes = [];
            
            for i = 1:length(obj.fluxBottom)
                % j = ceil(i/averageLength);
                %
                %if mod(i, averageLength) == 0
                %
                %    movingAverageFlux(j) = sum(fluxBottom(i+1-averageLength:i))/averageLength;
                %    movingAverageTimesteps(j) = timesteps(i);
                %    movingAverageTimes(j) = times(i);
                % end
                
                %Calculating moving average considering all fluxes within some specified
                %previous time
                j = i;
                if (obj.times(i) - obj.times(1) ) > timeSpan
                    fluxSum = 0;
                    N = 0;
                    for k = 1:i
                        if (obj.times(i) - obj.times(k)) < timeSpan
                            fluxSum = fluxSum + obj.fluxBottom(k);
                            N=N+1;
                        end
                    end
                    movingAverageFlux(j) = fluxSum/N;
                    movingAverageTimes(j) = obj.times(i);
                end
            end
            
            %Calculate change in moving averages
            deltaAverageFlux = (movingAverageFlux(2:end) - movingAverageFlux(1:end-1))./ movingAverageFlux(1:end-1);
            
            
            
            % log_iters = log10(1:10);
            %
            hFig = figure(1);
            set(hFig, 'Position', [300 300 1100 900])
            
            subplot(3, 1, 1);
            hold on;
            
            maxRate = max(1.2, max(obj.dTdt)*1.2);
            minRate = -1;
            hold on;
            
            yyaxis left;
            
            if plotRestarts
                for i=1:length(obj.fails)
                    plot([obj.fails(i) obj.fails(i)], [ minRate maxRate], 'g--');
                end
            end
            
            test = max(obj.dTdt);
            
            
            h1 = plot(obj.timesteps, obj.dTdt, '-b', 'LineWidth', 4);
            h2 = plot(obj.timesteps, obj.dSdt, '-y',  'LineWidth', 4);
            %    ylim([-1 max(max(dTdt), max(dSdt))*1.2]);
            hold off;
            
            yyaxis right;
            h3 = plot(obj.timesteps, obj.dUdt, 'LineWidth', 4);
            
            
            %plot(timesteps, dSdt./max(dSdt));
            %plot(timesteps, dUdt./max(dUdt));
            
            
            
            %     for i=1:length(spikes)
            %        plot([spikes(i) spikes(i)], [ minRate maxRate], 'r--');
            %     end
            
            legend([h1 h2 h3], {'H (left)', 'S (left)', 'U (right)'});
            
            %axis([timesteps(1) timesteps(end)+(timesteps(end)-timesteps(1))*0.2 minRate maxRate]);
            xlim([obj.timesteps(1) obj.timesteps(end)+(obj.timesteps(end)-obj.timesteps(1))*0.4]);
            grid on; box on;
            xlabel('timestep');
            ylabel('rate');
            
            %figure(2);
            subplot(3, 1, 2);
            %yyaxis left;
            hold on;
            plot(obj.timesteps, log10(obj.dTdt), '-b', 'LineWidth', 4);
            plot(obj.timesteps, log10(obj.dSdt), '-y', 'LineWidth', 4);
            
            % yyaxis right;
            
            plot(obj.timesteps, log10(obj.dUdt), '-r',  'LineWidth', 4);
            hold off;
            legend({'H', 'S', 'U'});
            
            xlim([obj.timesteps(1) obj.timesteps(end)+(obj.timesteps(end)-obj.timesteps(1))*0.4]);
            grid on; box on;
            xlabel('timestep');
            ylabel('log10(rate)');
            
            
            subplot(3, 1, 3);
            
            
            if containsFs
                
                yyaxis left;
                hold on;
                plot(obj.times, obj.fluxTop,'-', 'LineWidth', 4);
                plot(obj.times, obj.fluxBottom,'--', 'LineWidth', 4);
                %ylim([-7 7]);
                
                
                hold off;
                
                yyaxis right;
                
                fractionalFluxDiff = abs((obj.fluxTop - obj.fluxBottom)./obj.fluxTop);
                plot(obj.times, log10(fractionalFluxDiff));
                %hold on;
                
                %plot(movingAverageTimes, movingAverageFlux, '-m',  'LineWidth', 3);
                %hold off;
                
                
                legend('Fs top', 'Fs bottom', 'log(fractional diff)');
                xlim([obj.times(1) obj.times(end)+(obj.times(end)-obj.times(1))*0.4]);
                grid on; box on;
                xlabel('time');
                
                
            elseif containsNu
                
                plot(obj.times, obj.Nus);
                grid on; box on;
                ylabel('Nu');
            end
            
            
            figure();
            plot(obj.times, obj.minimumSl);
            box on; grid on;
            
            %       if plotRestarts
            % for i=1:length(fails)
            %       plot([fails(i) fails(i)], [ minRate maxRate], 'g--');
            % end
            %    end
            
        end
        
        function plot(obj, var)
            figure();
            plot(obj.timesteps, var);
            grid on; box on;
        end
        
    end
    
end









