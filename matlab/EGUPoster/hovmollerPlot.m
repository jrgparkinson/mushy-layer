%function hovmollerPlot(folder)
close all;
clear all;

folder = 'T-10';

centralPeak = [0 2.0];
%centralPeak = [0.75 1.25]; % -15, -10
%centralPeak = [0.6 1.0]; % -20
%numPeaks = 1;

numPeaks = 4;

centralPeak = round(centralPeak/2.0*256);
centralPeak = max(centralPeak, 1);

thisDir = strrep(mfilename('fullpath'), 'hovmollerPlot', '');

%local_base_dir = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/matlab/EGUPoster/';
remote_base = '/home/parkinsonjl/mnt/sharedStorage/middleton/';

local_base_dir = thisDir; %'/home/parkinsonj/convection-in-sea-ice/MushyLayer/matlab/EGUPoster/';
%remote_base = '/network/group/aopp/oceans/AW002_PARKINSON_MUSH/middleton/';

saveFiles = false;

TESTING = false;



output_dir = [local_base_dir, folder, '/'];
hovmollerFile = [output_dir, 'hovmoller.mat'];

% Download data
%downloadData(output_dir, folder, remote_base);

% Diagnostics for this simulation:
load([output_dir, 'diagConcat.mat']);

if exist(hovmollerFile, 'file') == 2
    load(hovmollerFile)
    freezingFront = max(freezingFront) - freezingFront;
else
    
    
    % Make frames from plotdata
    [files, frames] = getFiles(output_dir, 'toPlot');
    frames = sort(frames);
    %frames = getFiles(output_dir, 'videoFrame');
    
    numFiles = length(files);
    
    if TESTING
        numFiles = 2;
    else
        set(0, 'DefaultFigureVisible', 'on');
    end
    
    
    
    domWidth = 256;
    
    %numFiles = 50;
    
    step = 2;
    numDataPoints = floor(numFiles/step);
    
    % This contains a single horizontal strip for each timestep
    hovmoller = NaN*ones(numDataPoints, domWidth);
    times = NaN*ones(numDataPoints, 1);
    freezingFront = NaN*ones(numDataPoints, 1);
    
    t_min = 0;
    load([output_dir, 'toPlot',num2str(frames(end)) , '.mat']);
    t_max = t;
    
    h = figure();
    set(h, 'Position', [200 200 400 800]);
    
    for i=1:numDataPoints
        
        frame = frames(i*step);
        
        file = ['toPlot', num2str(frame), '.mat'];
        
        
        fprintf('Processing %s \n', file);
        
        
        
        dataFile = [output_dir, file];
        load(dataFile);
        
        fprintf('t =  %1.2f \n', t);
        
        if i > 2 && t < times(i-1)
            fprintf('Time seems to be going backwards - skip file \n');
            continue
        end
        
        % Find vertical position we want
        
        
        [mlx, mly] = find(porosity < 1);
        if mly
            mlMin = min(min(mly));
        else
            mlMin = size(Z, 1);
        end
        
        
        
        
        slice_y_i = mlMin - 10;
        
        sl_slice = squeeze(Sl(:, slice_y_i));
        
        hovmoller(i, :) = sl_slice;
        times(i) = t;
        freezingFront(i) = mlMin;
        avMushyPorosity(i) = mean(mean(porosity(porosity<1)));
        avMushySl(i) = mean(mean(Sl(porosity<1)));
        
        x = squeeze(X(1,:));
        
        [X, T] = meshgrid(x,times);
        
        pcolor(X, getTime(T), getSalnity(hovmoller));
        xlabel('x'); ylabel('t');
        
        c = colorbar('Location', 'eastoutside');
        c.Label.String = 'S_l';
        caxis([-1 -0.95]);
        
        ylim([t_min, t_max]);
        
        drawnow;
        %pause(0.1);
        
        
        
    end
    
    
    save(hovmollerFile, 'X', 'T', 'hovmoller', 'freezingFront',...
        'avMushyPorosity', 'avMushySl');
    
    
end


h2 = figure();
set(h2, 'Position', [200 200 600 800]);
pcolor(X, getTime(T), getSalinity(hovmoller));
axHov = gca;
axHov.Position = [0.15 0.13 0.65 0.83];
xlabel('x (cm)'); ylabel('t (s)');
%set(gca, 'Ydir', 'reverse');

c = colorbar('Location', 'eastoutside');
c.Label.String = 'S_l (g/kg)';


if saveFiles
h2.InvertHardcopy = 'off';
h2.Color = 'white';
  set(h2,'Units','Inches');
  pos = get(h2,'Position');
  set(h2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
  
  filename = ['/home/parkinsonjl/Documents/Talks/TransientSaltFluxes/hovmuller', folder];
  print(h2,[filename, '.pdf'],'-dpdf','-r0')
print(h2,[filename, '.png'],'-dpng','-r400')
end

%caxis([-1 -0.95]);


% Now do some analysis;
%Y = fft2(hovmoller);
%h3 = figure();
%imagesc(log10(abs(fftshift(Y))));
%imagesc((abs(fftshift(Y))));
%cfft = colorbar();
%title('log10(Fourier transform)');

% Peak finding to track plumes


h4=figure()

numTimeSteps = size(hovmoller, 1);
peakPositions = NaN*ones(numTimeSteps, numPeaks);

x = X(1, :);
dx = x(2)-x(1);

for t_i = 1:numTimeSteps
   Sl_slice = squeeze(hovmoller(t_i, :));
   
   % Find peak positions
   [pks, locs] = findpeaks(Sl_slice(centralPeak(1):centralPeak(2)));
   
   if length(locs) == numPeaks
       peakPositions(t_i, :) = locs*dx;
   end
   
   %plot(Sl_slice);
   %ylim([-1 max(max(hovmoller))*0.9]);
   %drawnow;
   
end

clf;

times = squeeze(T(:, 1));

hold on;
for i = 1:size(peakPositions, 2)
plot(peakPositions(:, i), T(:, 1));
end
hold off;

ylabel('t');
xlabel('x');


% Now do some fourier analysis
% Find frequency and amplitude of peaks

multiplepeakPositions = peakPositions;


peakPositions = sum(peakPositions, 2)/size(peakPositions, 2);

numSignals = size(peakPositions, 2);
for i=1:numSignals
   plumePos = squeeze(peakPositions(:, i)); 
   
   meanPos = nanmean(plumePos);
   
   [pks, locs] = findpeaks(plumePos, 'MinPeakHeight', meanPos);
   
   % Also need to find minima
   inverseData = max(plumePos)-plumePos;
   [pksMin, locsMin] = findpeaks(inverseData, 'MinPeakHeight', nanmean(inverseData));
   pksMin = max(plumePos)-pksMin;
   
   figure();
   findpeaks(plumePos, 'MinPeakHeight', meanPos);
%    
   figure();
   findpeaks(inverseData, 'MinPeakHeight', nanmean(inverseData));
   
   % Extract periods from these locations
   periodMax = [];
   periodMaxTimes = [];
   for j=1:length(locs)
       if j < length(locs)
        periodMax(j) = times(locs(j+1)) - times(locs(j));
        periodMaxTimes(j) = times(locs(j));
       end
       
   end
   
    periodMin = [];
   periodMinTimes = [];
   for j=1:length(locsMin)
       if j < length(locsMin)
        periodMin(j) = times(locsMin(j+1)) - times(locsMin(j));
        periodMinTimes(j) = times(locsMin(j));
       end
       
   end
   
   % combine  period calculations from maxima and minima
   periodTimes = [periodMaxTimes periodMinTimes];
   periodTimes = sort(periodTimes);
   periods = NaN*periodTimes;
   for j=1:length(periodTimes)
       thisTime = periodTimes(j);
       max_i = find(periodMaxTimes == thisTime);
       if max_i > 0
           
           periods(j) = periodMax(max_i);
       else
           min_i = find(periodMinTimes == thisTime);
           periods(j) = periodMin(min_i);
       end
   end
   
   amplitudes = [];
   amplitudeTimes = [];
   % Now process peaks. For each max peak, find the following min peak
   for j=1:length(pks)
       maxPeakPos = locs(j);
       
       for k=1:length(pksMin)
           if locsMin(k) > maxPeakPos
               amplitudes(end+1) = pks(j) - pksMin(k);
               amplitudeTimes(end+1) = times(locs(j));
               break;
           end
       end
       
   end
   
   % remove anomalies
   periods(periods < nanmean(periods)/2) =  NaN;
   
   hSummary = figure();
   set(hSummary, 'Position', [100 100 1200 800]);
   
        m = 1; n = 5;
        subplot(m, n, 1);
        plot(plumePos, times);
        title('x position');
        ylabel('t');

        subplot(m, n, 2);
        plot(periods, periodTimes, '-x');
        title('Period');
        ylabel('t');

        subplot(m, n, 3);
        plot(amplitudes, amplitudeTimes, '-x');
        title('Amplitude');
        ylabel('t');
        
        subplot(m, n, 4);
        plot(freezingFront, times, '-');
        title('Freezing front');
        ylabel('t');
        
         subplot(m, n, 5);
         
         growthRate = 1./freezingFront;
         %Interpolate growth rate onto period data
         
         interpGrowthrate = interp1(times, growthRate, periodTimes);
         scaledPeriod = periods.*interpGrowthrate;
         
            plot(scaledPeriod, periodTimes, '-');
            title('Period/Freezing front');
            ylabel('t');
   
            
            
            fprintf('Mean amplitude: %1.4f \n', nanmean(amplitudes));
            
end





figure();

plot(concat_Si, concat_time);
xlabel('$S_{ice}$')
ylabel('t');



