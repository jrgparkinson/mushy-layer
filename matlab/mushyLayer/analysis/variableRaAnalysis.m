%clear all;
close all;

set(0, 'defaultlinelinewidth',2);
set(0, 'defaultaxeslinewidth',2);
set(0, 'defaultpatchlinewidth',2);
set(0, 'defaultAxesFontSize',22);

% Length scale (metres)
L = 2;

baseDir = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/';

dataFolder = '/media/parkinsonjl/FREECOM HDD/';

RaCs = [75, 100, 150, 200, 250, 300, 400, 500, 600, 700, 800];
%RaCs = [600];

doChannels = true;

if exist('optimalFluxes', 'var') == 0
    optimalFluxes = NaN*RaCs;
end

if exist('nonOptimalFluxes', 'var') == 0
    nonOptimalFluxes = NaN*RaCs;
end

if exist('error', 'var') == 0
    error = NaN*RaCs;
end

if exist('optimalWidths', 'var') == 0
    optimalWidths = NaN*RaCs;
end

if exist('optimalChannelWidth', 'var') == 0
    optimalChannelWidth = NaN*RaCs;
end

if exist('optimalChannelDepth', 'var') == 0
    optimalChannelDepth = NaN*RaCs;
end

if exist('channelErr', 'var') == 0
    channelErr = NaN*RaCs;
end

%hStates = figure();

for i=1:length(RaCs)
    if isnan(optimalFluxes(i))
        
        
        
        % Try and find the optimal solute flux
        RaC = RaCs(i);
        
        % Find all files with this RaC in this folder:
        
        insulating_folder = [dataFolder , 'mushyLayerLowC-upperBranch-insulating'];
        all_files = dir(insulating_folder);
        files = {};
        fluxes = [];
        widths = [];
        output_dir = {};
        plot_prefix = {};
        frame = [];
        
        chanWidths = [];
        chanDepths = [];
        
        for file_i=1:length(all_files)
            folder = all_files(file_i).name;
            if findstr(folder, ['CR1.25RaC',num2str(RaCs(i)),'Le200ChiCubedPermeabilitypts'])
                % files(end+1) = filename;
                
                parts = strsplit(folder, 'pts');
                
                if length(parts) > 1
                    parts2 = strsplit(parts{2}, '-');
                    width = str2num(parts2{1});
                    
                    poutFile = [insulating_folder, '/', folder, '/pout.0'];
                    thisPout = Pout(poutFile);
                    if thisPout.steadyState && length(thisPout.fluxBottom) > 0
                        thisFlux =  thisPout.fluxBottom(end);
                        widths(end+1) = width;
                        fluxes(end+1) = thisFlux;
                        
                        if doChannels
                       
                        % Get the steady state file
                        thisFullFolder = [insulating_folder, '/', folder];
                        files = dir([thisFullFolder, '/*.2d.hdf5']);
                        finalFile = files(end);
                        
                        [mat,tok,ext]  = regexp(finalFile.name, '(.*-)(\d+)\.2d\.hdf5', 'match', ...
                            'tokens', 'tokenExtents');
                        
                        
                        thisTok = tok{1};
                        
                        %steadyStateFile(end+1) = finalFile;
                        
                        frame  = str2num(thisTok{2});
                        output_dir = [thisFullFolder, '/']; %thisFullFolder;
                        plot_prefix = thisTok{1};
                        
                         %Calculate channel geometry for all states with
                         %this RaC
                        dim = 2; subcycled=true;
                        output =  MushyLayerOutput(dim, frame, output_dir,...
                        plot_prefix, subcycled);
        
                        [chanWidths(end+1), chanDepths(end+1)] = output.channelGeometry();
                        
                        end
                        
                    end
                    
                end
                
            end
        end
        
        if length(fluxes) > 0
            
            % Fit quadratic to fluxes near the maximum
            maxFluxi = find(fluxes == max(fluxes));
            
            singleOptimali = maxFluxi;
            if length(maxFluxi) > 1
                singleOptimali = maxFluxi(1);
            end
            
            % use min() and max() because maxFluxi may be a vector
            mini = max(1,              min(maxFluxi) - 2);
            maxi = min(length(fluxes), max(maxFluxi) + 2);
            
            while isnan(fluxes(mini))
                mini = mini + 1;
            end
            
            while isnan(fluxes(maxi))
                maxi = maxi-1;
            end
            
            
            % fit quadratic polynomial
            dx = 0.2/128;
            p = polyfit(widths(mini:maxi)*dx,fluxes(mini:maxi),2);
            
            %  c-b^2/(4*a)
            maximumFlux =  p(3) - p(2)^2 / (4*p(1));
            optimalFluxes(i) = maximumFlux;
            
            %  -b/(2*a)
            optimalWidths(i) = -p(2)/(2*p(1));
            
            % These should be decreasing
            if i > 1 && (optimalWidths(i) > optimalWidths(i-1))
                optimalWidths(i) = optimalWidths(i) / 2;
            end
            
            if doChannels
           
            % Fit for channel width
            %pWidth = polyfit(widths(mini:maxi)*dx,chanWidths(mini:maxi),2);
            optimalChannelWidth(i) =  interp1(widths(mini:maxi)*dx, chanWidths(mini:maxi),optimalWidths(i)); 
            
            % Fit for channel depth
            %pDepth = polyfit(widths(mini:maxi)*dx,chanDepths(mini:maxi),2);
            optimalChannelDepth(i) =   interp1(widths(mini:maxi)*dx, chanDepths(mini:maxi),optimalWidths(i)); 
            
            
            if singleOptimali > 1 && singleOptimali < length(chanWidths)
               
            channelErr(i) = max( abs(chanWidths(singleOptimali) - chanWidths(singleOptimali-1)), ...
                abs(chanWidths(singleOptimali) - chanWidths(singleOptimali+1)));
            elseif singleOptimali > 1
                     channelErr(i) =abs(chanWidths(singleOptimali) - chanWidths(singleOptimali-1));
            elseif  singleOptimali < length(chanWidths)
                 channelErr(i) = abs(chanWidths(singleOptimali) - chanWidths(singleOptimali+1));
            else
                
            end
            
            %dim = 2; subcycled=true;
            %output =  MushyLayerOutput(dim, frame(maxFluxi), output_dir{maxFluxi},...
            %plot_prefix{maxFluxi}, subcycled);
        
            % Calculate optimal channel geometry by fitting quadratic
            %[optimalChannelWidth(i), optimalChannelDepth(i)] = output.channelGeometry();
            
            end
            
        end
        
        
    end
    
    
end

% hack because something has gone wrong
optimalWidths(5:7) = optimalWidths(5:7)*2;

% poutNames = {'Ra200', 'Ra250', 'Ra300'};
lineStyles = {'-', '--',':', '-.'};
lineColours = {'r', 'b', 'g', 'm', 'black'};
lineColours = {'m', 'c', 'r', 'g', 'b', 'k'};
h = figure();
set(h, 'Position', [200 200 600 400]);


legendStr = {};



hold on;
%e = errorbar(RaCs, optimalFluxes-1, error, 'x-', 'LineWidth', 2.0);
%e.MarkerSize = 10;
%hold on;
plot(RaCs, optimalFluxes-1, 'x-');
%plot(RaCs, nonOptimalFluxes-1, 'o-');

%hold off;
axis([0 850 0 0.2]);

%legend({'Upper branch', 'Lower branch'});

%title('Steady state vertical solute flux');

xlabel('$Ra_s$', 'Interpreter','latex');
ylabel('$F_0$', 'Interpreter','latex');

set(get(gca,'YLabel'),'Rotation',0, ...
    'Units', 'Normalized', 'Position', [-0.2, 0.65, 0])

set(gca, 'position', [0.25 0.22 0.7 0.73]);

ax = gca;
ax.YTick = [0,0.05, 0.1, 0.15, 0.2];
ax.YTickLabel ={'0.0','', '0.1','', '0.2'};

%legend({'Optimal flux', 'Largest flux found so far'});

%grid on;
box on;

%For saving figure as PDF
% h =gcf;
% set(h,'Units','Inches');
% pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(h,'optimalFlux.pdf','-dpdf','-r0')



if doChannels
    
inverseChanErr = channelErr./(optimalChannelWidth.^2);

figure();
yyaxis left;
%errorbar(RaCs, optimalChannelWidth, channelErr); ylabel('a');
errorbar(RaCs, 1./optimalChannelWidth, inverseChanErr);
ylabel('1/a');

%yyaxis right;
%plot(RaCs, optimalChannelDepth);
%ylabel('d');

box on;
xlabel('Ra_S'); 

end


 %set(hStates,'Units','Inches');
 %pos = get(hStates,'Position');
 %set(hStates,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
 %print(hStates,['optimalStates.pdf'],'-dpdf','-r0')
 
 
 %
% figure(); plot(log(RaCs), optimalFluxes-1); xlabel('log(RaC)'); ylabel('F_0 - V');
%figure(); plot(log(RaCs), optimalFluxes-1, 'x-'); xlabel('log(RaC)'); ylabel('F_0 - V');
%figure(); plot(1./optimalWidths, optimalFluxes-1, 'x-'); xlabel('log(RaC)'); ylabel('F_0 - V');