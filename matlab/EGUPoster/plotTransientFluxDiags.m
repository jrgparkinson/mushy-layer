function plotTransientFluxDiags
fontsize = 16;
set(0, 'defaultAxesFontSize',fontsize);

clear all;
%close all;
data_dir = getDataDir('middleton/');
data_dir = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/matlab/EGUPoster/';
%run = 'CR1.179RaC80Le200KozenyPermeabilitypts256-S3.5-TB-20.0-R0.13-domWidth2.5-stoppedEarlyWavey';
%run = 'CR1.179RaC80Le200KozenyPermeabilitypts256-S3.5-TB-15.0-R0.13-0';

%runs = {'CR1.179RaC800Le200KozenyPermeabilitypts128-S3.5-TB-15.0-R0.013-domWidth1.0-0',...
%    'CR1.179RaC800Le200KozenyPermeabilitypts128-S3.5-TB-20.0-R0.013-domWidth1.0-0',};
%leg = {'T=-15','T=-20'};

timeSlices = [3 20 30 50 60];

colors = [0.635 0.078 0.184; % warmest
    %0.85 0.325 0.098; % medium
    0.929 0.694 0.125;
    0 0.447 0.7410 % coldest
    ];

runs = {'T-10',...
    'T-15',...
    'T-20'};
leg = {'$T_a=-10^\circ$C', '$T_a=-15^\circ$C','$T_a=-20^\circ$C'};


% Need some new vectors so we can combine diag files
time = []; %NaN*linspace(3,1e5);
FsAv = []; %NaN*linspace(3,1e5);
FsMax = []; % NaN*linspace(3,1e5);

times2 = [];
Si2 = [];
Vi2 = [];
Sbox2 = [];

for i=1:length(runs)
    diagFile = [data_dir, runs{i} '/diagConcat.mat'];
    if exist(diagFile, 'file') == 2
        load(diagFile);
        l = length(concat_time);
        times2(i, 1:l) = concat_time;
        Si2(i, 1:l) = concat_Si;
        Vi2(i, 1:l) = concat_Vi;
        Sbox2(i, 1:l) = concat_Sbox;
    else
        times2(i, :) = [NaN];
        Si2(i, :) = [NaN];
        Vi2(i, :) = [NaN];
        Sbox2(i, :) = [NaN];
    end
end


loadSavedData = true;
diagFile = '/diagnostics.out';

maxTime = 70;

if loadSavedData
    load([data_dir, 'timeseries.mat']);
else
    
    FsField = 'L1FsVertFluid'; %'L0FsVertFluid'
    %FsField = 'L0FsVertFluid';
    
    
    diag1 = getDiagnostics([data_dir, runs{1}, diagFile]);
    d(1:length(runs)) = diag1;
    
    time(1, 1:length(diag1.time)) = diag1.time;
    FsAv(1, 1:length(diag1.time)) = diag1.averageLiquidSalinity;
    FsMax(1, 1:length(diag1.time)) = diag1.(FsField);
    
    for i=2:length(runs)
        diags = getDiagnostics([data_dir, runs{i}, diagFile]);
        d(i) = diags;
        
        time(i, 1:length(diags.time)) = diags.time;
        FsAv(i, 1:length(diags.time)) = diags.averageLiquidSalinity;
        FsMax(i,1:length(diags.time)) = diags.(FsField);
    end
    
    for j=3:4
        
        for i=1:length(runs)
            diagFile = [data_dir, runs{i}, '/diagnostics',num2str(j),'.out'];
            fprintf('%s \n', diagFile);
            diag1 = getDiagnostics(diagFile);
            
            if isfield(diag1, 'averageLiquidSalinity')
                
                yi = find(time(i,:)==0);
                
                start = min(find(time(i,:)==0));
                
                if length(start) == 0
                    start = length(time(i, :));
                end
                
                k = 1;
                prevSl =  FsAv(i, start-k);
                while isnan(prevSl)
                    prevSl = FsAv(i, start-k);
                    k = k + 1;
                    
                end
                
                time(i, start:start+length(diag1.time)-1) = diag1.time;
                offset =  - diag1.averageLiquidSalinity(1) + prevSl;
                FsAv(i,  start:start+length(diag1.time)-1) = diag1.averageLiquidSalinity + offset;
                FsMax(i,  start:start+length(diag1.time)-1) = diag1.(FsField);
                
                fprintf('Offset = %1.5f (old Sl = %1.5f, new Sl = %1.5f) \n', offset, FsAv(i, start-1), diag1.averageLiquidSalinity(1));
                
            end
        end
        
    end
    
    % Trying to remove dodgy valu
    time(time==0) = NaN;
    
    save([data_dir, 'timeseries.mat'], 'time', 'FsAv', 'FsMax');
    
    
end



%folder = [data_dir, run, '/'];
%diagFile = [folder, 'diagnostics.out'];
%d = getDiagnostics(diagFile);

%timescale = 800;
%time = d.time*timescale/60;


posLower = [0.12 0.1 0.85 0.44];
posUpper = [0.12 0.54 0.85 0.44];

h =figure();
set(h, 'Position', [200 200 800 400]);
m = 2; n=1;

%subtractOff = [80000 40000 000];
subtractOff = [80000 40000 5000];

subplot(m, n, 1);
hold on
for i=1:length(runs)
    %t = getTime(squeeze(time(i, 1:end-subtractOff(i))));
    %Sl = getSalinity(squeeze(FsAv(i, 1:end-subtractOff(i))));
    %plot( t(~isnan(Sl)), Sl(~isnan(Sl)), '-', 'Color', colors(i, :))
    
    
    t = getTime(squeeze(times2(i, :)));
    Si = getSalinity(squeeze(Si2(i, :)));
    mi = squeeze(Vi2(i, :))*(0.08)*1000; % mass of water/m^2
    %Delta S (g/m^2)= Salinity (g/kg) * mass of water (kg/m^2)
    deltaS =  (35- Si).*mi; %g/m^2
    deltaS = deltaS/1000;
    
     dSdt = gradient(deltaS,t);
     
     doPlot = (~isnan(deltaS)).*(t ~= 0);
     
     dSdt = min(dSdt, 0.018);
    
   % plot( t(~isnan(F)), BSF(~isnan(F)), '-', 'Color', colors(i, :))
    plot( t(doPlot==1), (max(dSdt(doPlot==1),0)), '-', 'Color', colors(i, :))
    
    
    %plot( getTime(time.'), getSalinity(FsAv.'), '-')
end
ax = gca;
yl = ax.YLim;
yl = [-0.001 0.02];
plotTimeSlices(timeSlices, yl);
ax.YLim = yl;

%subtractOff = [80000 40000 5000];

hold off
box on;
%xlabel('$t$ (mins)');
ax.XTick = 0:10:70;
ax.XTickLabels={'', '', '', '', '', '', '', '', '', '', ''};

ax.YTick = [0 0.01 0.02];
ax.YTickLabels = {'0', '0.01', '0.02'};

ylabel('$F_s$ (kg m$^{-2}$ s$^{-1}$)');
%ax.YLim = [-0.0 0.05];

%ylabel('$\bar{S}_l$ (g/kg)');
%ax.YTick = [34 35 36 37 38 39 40 41 42 43];
%ax.YTickLabels = {'', '35', '', '', '', '', '40','', '', ''};
%ylabel('$\bar{S}_i$ (g/kg)');
legend(leg, 'Location', 'northeast');

ax.YLim = yl;


%ax.YTick = [34 35 36 37 38 39 40 41 42 43];
%ax.YTickLabels = {'', '35', '', '', '', '', '40','', '', ''};




ax.XLim = [0 maxTime];

ax.Position = posUpper;

subplot(m, n, 2);
hold on
for i=1:length(runs)
%     F = getFlux(squeeze(FsMax(i, 1:end-subtractOff(i)))); %(m/s) (g/kg)
%     t = getTime(squeeze(time(i, 1:end-subtractOff(i))));
%     
%     A = 0.02^2; %m^2
%     rho = 1000; %g/kg
%     BSF = F*rho; %g/m^2/s
%     BSF = BSF/1000; %kg/m^2/s

 t = getTime(squeeze(times2(i, :)));
    Si = getSalinity(squeeze(Si2(i, :)));
    mi = squeeze(Vi2(i, :))*(0.08)*1000; % mass of water/m^2
    %Delta S (g/m^2)= Salinity (g/kg) * mass of water (kg/m^2)
    deltaS =  (35- Si).*mi; %g/m^2
    deltaS = deltaS/1000;
    
    doPlot = (~isnan(deltaS)).*(t ~= 0);
    
    
    plot(t(doPlot==1),deltaS(doPlot==1), '-', 'Color', colors(i, :))
    
   % dSdt = gradient(deltaS,t);
    
   % plot( t(~isnan(F)), BSF(~isnan(F)), '-', 'Color', colors(i, :))
   % plot( t(~isnan(dSdt)), dSdt(~isnan(dSdt)), '-', 'Color', colors(i, :))
    %plot( getTime(time.'), getFlux(FsMax.')*100, '-')
end
ax = gca;
yl = ax.YLim;
yl(1) = 0;
yl(2) = 0.24;
plotTimeSlices(timeSlices, yl);
ax.YLim = yl;

ax.Position = posLower;

hold off
box on;
xl = xlabel('$t$ (mins)'); 

ylabel('$\Delta S$ (kg m$^{-2}$)');


ax.XLim = [0 maxTime];
ax.XTick = 0:10:70;
ax.XTickLabels = {'0', '10', '20', '30', '40', '', '60', '70'};

xlPos = xl.Position;
xl.Position = [50 -0.01]; %;[xlPos(1)-30.0 xlPos(2)+0.017];


% h =gcf;
% h.InvertHardcopy = 'off';
% h.Color = 'white';
%   set(h,'Units','Inches');
%   pos = get(h,'Position');
%   set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%   filename = ['/home/parkinsonjl/convection-in-sea-ice/figures', '/transientDiagnostics'];
%   print(h,[filename, '.pdf'],'-dpdf','-r0')
% print(h,[filename, '.png'],'-dpng','-r400')





end

function F = getFlux(dimlessFlux)
fluxScale = 0.0029;
%fluxScale= 1.2500e-05;
F = dimlessFlux*fluxScale;
end

function S = getSalinity(dimlessSalinity)
Se = 233;
deltaS = 233-35;

S = dimlessSalinity*deltaS + Se;
end

function t = getTime(dimlessTime)
timescale = 800; %seconds
t = dimlessTime*timescale/60; % minutes
%t = dimlessTime;
end

function plotTimeSlices(times, ylim)

for i = 1:length(times)
    plot([times(i) times(i)], ylim, '--k');
    
end
end

