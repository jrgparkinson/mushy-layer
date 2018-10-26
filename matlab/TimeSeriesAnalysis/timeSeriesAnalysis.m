function timeSeriesAnalysis

close all; clear all;

figure();
hold on;

temps = {'-10','-15','-20'};


for temp_i = 1:length(temps)
    
output_dir = ['/home/parkinsonjl/convection-in-sea-ice/MushyLayer/',...
    'matlab/EGUPoster/T',temps{temp_i},'/']; 
load([output_dir, 'diagConcat.mat']);


t = getTime(concat_time);
mi = concat_Vi*(0.08)*1000; % mass of water/m^2
deltaS =  (35- getSalinity(concat_Si)).*mi; %g/m^2

saltRelease = deltaS/1000;
dtOrig = t(end)-t(end-1);
%resampleSpacing = 100/ ( (t(end)-t(end-1)));
%saltReleaseEven = resample(saltRelease, t, resampleSpacing);

saltReleaseEven = resample(saltRelease, t);
tEven = linspace(min(t), max(t), length(saltReleaseEven));

%tEven = linspace(min(t), max(t), length(t)*10);
%saltReleaseEven = interp1(t,saltRelease,tEven,'spline');


saltFlux = gradient(saltRelease,t);
saltFluxEven = gradient(saltReleaseEven, tEven);

saltFluxEven = resample(saltFlux, t);
%saltFluxEven = mean(saltFlux) + sin(2*pi*tEven); % + 0*saltFluxEven;
%tEven = linspace(min(t), max(t), length(saltReleaseEven));

leg= {'Uneven', 'Even'};

% figure();
% m = 2; n=1;
% subplot(m, n, 1);
% hold on;
% plot(t, saltRelease);
% plot(tEven, saltReleaseEven);
% hold off;
% 
% legend(leg, 'Location', 'eastoutside');
% 
% subplot(m,n,2);
% hold on;
% plot(t, saltFlux);
% plot(tEven, saltFluxEven);
% 
% hold off;
% legend(leg, 'Location', 'eastoutside');


Fs = 1/(tEven(2)-tEven(1));

Y = fft(saltFluxEven);
L = length(saltFluxEven);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);       
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;



plot(f,P1) 


end

title('Single-Sided Amplitude Spectrum')
xlabel('f (Hz)')
ylabel('|P1(f)|')

legend(temps);


end

function t = getTime(dimlessTime)
timescale = 800; %seconds
t = dimlessTime*timescale/60; % minutes
%t = dimlessTime;
end


function S = getSalinity(dimlessSalinity)
Se = 233;
deltaS = 233-35;

S = dimlessSalinity*deltaS + Se;
end
