% Extract data from pout file to look at how we approach (or don't
% approach) steady state. Mainly useful for seeing if we oscillate about a
% steady state.

clear all;
close all;

%Options
plotRestarts = false;
timestepInterval = 1;
minTStep = 1;


pout_file = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/pout-16.0';
pout_file = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/pout-32.0';
%pout_file = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/pout-64.0';

%pout_file = '/home/parkinsonjl/convection-in-sea-ice/test/pout.0';
%pout_file = '/home/parkinsonjl/convection-in-sea-ice/test/mushyLayer/128-0/pout.0';
%pout_file = '/home/parkinsonjl/gyre/mushyLayer2/64-0/pout.0';
%pout_file = '/home/parkinsonjl/gyre/mushyLayer-insulating/128-0/pout.0';
%pout_file = '/home/parkinsonjl/gyre/rayleighBenardAMR-adaptive/16-1/pout.0';

Ra = 100;
Nu_caltagirone = containers.Map({0,39,50,100,200,250, 300}, [1 1  1.45 2.651 3.813 4.199 4.523]);

% For Pr = 0.71, wavenumber = 3.117
%For Ra = 4000 this is just an educated guess
Ra = 4000;
Nu_cleverBusse = containers.Map({2000, 2500, 3000, 4000, 5000, 10000, 20000, 30000, 50000}, ...
                                [1.212 1.475 1.663 1.985 2.116 2.661 3.258 3.662 4.245]);

%Nu_analytic = Nu_caltagirone(Ra);

Nu_analytic = Nu_cleverBusse(Ra);

fileID = fopen(pout_file,'r');

file_contents = fileread(pout_file);

%Options:


realNumberRegex = '-?\d+.\d+e[+-]\d+';

%Darcy velocity
%Temperature
pattern = ['d\/dt \(Temperature\) = (', realNumberRegex,')'];
%pattern = 'd\/dt \(Temperature\)';
[mat,tok,ext]  = regexp(file_contents, pattern, 'match', ...
               'tokens', 'tokenExtents');
           
% pattern2 = 'coarse time step\s+(\d+)';           
% [mat2,tok2,ext2]  = regexp(file_contents, pattern2, 'match', ...
%                'tokens', 'tokenExtents');   
           
pattern3 = ['Nusselt number \(level 0\) = (', realNumberRegex,')'];
[mat3,tok3,ext3]  = regexp(file_contents, pattern3, 'match', ...
               'tokens', 'tokenExtents');  
           
 %          'AMRLevelMushyLayer::postTimeStep 0[\n]',...
%'Solute flux \(level 0\) = ',realNumberRegex, '[\n]',...
containsFs = false;
if (strfind(file_contents, 'Solute flux'))
    containsFs = true;
end

containsNu = false;
if (strfind(file_contents, 'Nusselt number'))
    containsNu = true;
end

%'Time =((?!coarse).)*[\n]',...

pattern4 = ['coarse time step\s+(\d+)\s+old time = (',realNumberRegex ,')[^\n]+[\n]'...
    'd\/dt \(Enthalpy\) = (',realNumberRegex,')[\n]' ... 
    'd\/dt \(Bulk concentration\) = (',realNumberRegex,')[\n]' ...
    'd\/dt \(Darcy velocity\) = (',realNumberRegex,')[\n]'];  


if containsFs
    pattern4 = ['Solute flux \(level 0\). Bottom = (',realNumberRegex, '), top = (',realNumberRegex,')[\n]',...
        pattern4];
end

pattern4 = ['Minimum liquid salinity:\s+(', realNumberRegex, ')[\n]', ...
    pattern4];

[mat4,tok4,ext4]  = regexp(file_contents, pattern4, 'match', ...
               'tokens', 'tokenExtents');
           
pattern5 =  'BaseLevelTGA:: WARNING: solver exitStatus == 4\n((?!coarse).)*[\n]coarse time step\s+(\d+)'; %[.\n]*AMR([.]*)\n
[mat5,tok5,ext5]  = regexp(file_contents, pattern5, 'match', ...
               'tokens', 'tokenExtents');
           
           pattern5 =  'Solver failed -- restarting((?!coarse).)*[\n]coarse time step\s+(\d+)'; %[.\n]*AMR([.]*)\n
[mat5,tok5,ext5]  = regexp(file_contents, pattern5, 'match', ...
               'tokens', 'tokenExtents');
           
           
%tok

dTdt = NaN*ones(length(tok), 1);
timestep = 1;


fails = [];

timesteps = NaN*ones(length(tok), 1);
times = [];

numTsteps = length(timesteps);

Nus = NaN*ones(length(tok), 1);
fluxBottom  = []; fluxTop = [];
dSdt =  NaN*ones(length(tok), 1); dUdt =  NaN*ones(length(tok), 1);
minimumSl = NaN*ones(numTsteps);

t_i = 1;

for i = 1:length(tok4)
  
    ind = 1;
    
    minimumSl(t_i) = str2double( tok4{i}{ind});
    ind = ind+1;
    
    if containsFs
        
    fluxBottom(t_i) = str2double( tok4{i}{ind});
    fluxTop(t_i) = str2double( tok4{i}{ind+1});
    ind = ind + 2;
    
    end
    
    tstep =  str2double(tok4{i}{ind});
    if mod(tstep, timestepInterval) ~= 0 || tstep < minTStep
        continue;
    end
    
    timesteps(t_i) = tstep;
    times(t_i) =str2double( tok4{i}{ind+1}); 
    dTdt(t_i) = str2double( tok4{i}{ind+2});
    dSdt(t_i) = str2double( tok4{i}{ind+3});
    dUdt(t_i) = str2double( tok4{i}{ind+4});
          
    
     
     if (length(tok3) >= i) 
     Nus(t_i) =  str2double(tok3{i}{1});
     else
         Nus(t_i) = NaN;
     end
     
     Nu_err = 100*(Nus(t_i) - Nu_analytic)./Nu_analytic;
          
     fprintf('timestep %d, rate = %s, Nusselt = %s, Nusselt error (percent) = %s \n', timesteps(t_i), dTdt(t_i), Nus(t_i), Nu_err);
        
    t_i = t_i + 1;
end



movingAverageFlux = [];
averageLength = 200;
timeSpan = 0.1;
movingAverageTimesteps = []; movingAverageTimes = [];

for i = 1:length(fluxBottom)
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
   if (times(i) - times(1) ) > timeSpan
       fluxSum = 0;
       N = 0;
   for k = 1:i
      if (times(i) - times(k)) < timeSpan
          fluxSum = fluxSum + fluxBottom(k);
          N=N+1;
      end
   end
     movingAverageFlux(j) = fluxSum/N;
     movingAverageTimes(j) = times(i);
   end
end

%Calculate change in moving averages
deltaAverageFlux = (movingAverageFlux(2:end) - movingAverageFlux(1:end-1))./ movingAverageFlux(1:end-1);


%identify solver fails
for i = 1:length(tok5)
    fails(i) = str2double( tok5{i}{2});
    %fprintf('failed timestep: %d \n', fails(i));
end

% log_iters = log10(1:10);
% 
hFig = figure(1);
set(hFig, 'Position', [300 300 1100 900])

subplot(3, 1, 1);
hold on;
    
maxRate = max(1.2, max(dTdt)*1.2);
minRate = -1;
hold on;

yyaxis left;

if plotRestarts
 for i=1:length(fails)
       plot([fails(i) fails(i)], [ minRate maxRate], 'g--');  
 end
end
    
 test = max(dTdt);
 
    
    h1 = plot(timesteps, dTdt, '-b', 'LineWidth', 4);
    h2 = plot(timesteps, dSdt, '-y',  'LineWidth', 4);
%    ylim([-1 max(max(dTdt), max(dSdt))*1.2]);
    hold off;
    
    yyaxis right;
    h3 = plot(timesteps, dUdt, 'LineWidth', 4);
    
    
    %plot(timesteps, dSdt./max(dSdt));
    %plot(timesteps, dUdt./max(dUdt));
    
    
   
%     for i=1:length(spikes)
%        plot([spikes(i) spikes(i)], [ minRate maxRate], 'r--');  
%     end
  
     legend([h1 h2 h3], {'H (left)', 'S (left)', 'U (right)'});
    
    %axis([timesteps(1) timesteps(end)+(timesteps(end)-timesteps(1))*0.2 minRate maxRate]);
    xlim([timesteps(1) timesteps(end)+(timesteps(end)-timesteps(1))*0.4]);
    grid on; box on;
xlabel('timestep');
ylabel('rate');

%figure(2);
subplot(3, 1, 2);
%yyaxis left;
hold on;
    plot(timesteps, log10(dTdt), '-b', 'LineWidth', 4);
    plot(timesteps, log10(dSdt), '-y', 'LineWidth', 4);
   
   % yyaxis right;
    
    plot(timesteps, log10(dUdt), '-r',  'LineWidth', 4);
     hold off;
     legend({'H', 'S', 'U'});
    
    xlim([timesteps(1) timesteps(end)+(timesteps(end)-timesteps(1))*0.4]);
    grid on; box on;
xlabel('timestep');
ylabel('log10(rate)');


subplot(3, 1, 3);


if containsFs

    yyaxis left;
hold on;
    plot(times, fluxTop,'-', 'LineWidth', 4);
     plot(times, fluxBottom,'--', 'LineWidth', 4);
    %ylim([-7 7]);
    

    hold off;
    
    yyaxis right;
    
    fractionalFluxDiff = abs((fluxTop - fluxBottom)./fluxTop);
    plot(times, log10(fractionalFluxDiff));
    %hold on;
    
    %plot(movingAverageTimes, movingAverageFlux, '-m',  'LineWidth', 3);
    %hold off;
    
    
    legend('Fs top', 'Fs bottom', 'log(fractional diff)');
     xlim([times(1) times(end)+(times(end)-times(1))*0.4]);
    grid on; box on;
xlabel('time');

    
elseif containsNu
    
    plot(times, Nus);
     grid on; box on;
     ylabel('Nu');
end


figure();
plot(times, minimumSl);
box on; grid on;
    
%       if plotRestarts
% for i=1:length(fails)
%       plot([fails(i) fails(i)], [ minRate maxRate], 'g--');  
% end
%    end     


   



