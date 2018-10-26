close all;

folders = {
    %'CR2.0RaC60Le100.0ChiCubedPermeabilitypts128-steady', ...
   % 'CR2.0RaC60Le200ChiCubedPermeabilitypts128-steady', ...
   % 'CR2.0RaC60Le300.0ChiCubedPermeabilitypts128-steady', ...
    'CR2.0RaC60Le800.0ChiCubedPermeabilitypts128-steady'
    };
Learr = [100, 200, 300, 800];

Learr = [800];

for f_i = 1:length(folders)
    
data_dir = ['/home/parkinsonjl/mnt/sharedStorage/saltDiffusionChannelWidth/'...
   folders{f_i}];
Le = Learr(f_i);
Ra = 60;

%Le = 200; Ra = 60;

getData = false;
if getData
ml = getFinalPlotFile(data_dir);


T = ml.dataForComp(ml.components.Temperature);
Sl = ml.dataForComp(ml.components.Liquidconcentration);
psi = ml.dataForComp(ml.components.streamfunction);
chi = ml.dataForComp(ml.components.Porosity);

end

[X,Y] = ml.grid();

dz = X(1,2) - X(1,1);
s = size(X);
Nz = s(1);
Nx = s(2);

channel_j = 150;

channel_side = 1;
if channel_side < 0
  channel_search_vector = 1:1:Nx;
else
    channel_search_vector = Nx:-1:1;
end

for n=1:length(channel_search_vector)
   i = channel_search_vector(n);
   
   if chi(i, channel_j) < 1
       channel_i = i;
       break;
   end
end


if channel_side < 0
    channel_pos = [1 channel_i];
else
    channel_pos = [channel_i Nx];
end
%channel_search_vector = 
%for i = channel_side


%figure();
%hold on;
%pcolor(chi.');
%plot(channel_pos, [channel_j, channel_j],'--');


Tchan = squeeze(T(channel_pos(1):channel_pos(2), channel_j));
Slchan = squeeze(Sl(channel_pos(1):channel_pos(2), channel_j));
psichan = squeeze(psi(channel_pos(1):channel_pos(2), channel_j));
xchan = X(channel_j, channel_pos(1):channel_pos(2));

% if channel_side > 0
%    dSdz =  (Sl(channel_pos(2), channel_j+1) - Sl(channel_pos(2), channel_j-1))/dz;
%    dTdz =  (T(channel_pos(2), channel_j+1) - T(channel_pos(2), channel_j-1))/dz;
% else
%    dSdz =  (Sl(channel_pos(1), channel_j+1) - Sl(channel_pos(1), channel_j-1))/dz;
%    dTdz =  (T(channel_pos(1), channel_j+1) - T(channel_pos(1), channel_j-1))/dz;
% end

dSdz =  ((Sl(channel_pos(1), channel_j+1) - Sl(channel_pos(1), channel_j-1))/dz + ...
   (Sl(channel_pos(2), channel_j+1) - Sl(channel_pos(2), channel_j-1))/dz )/2;
dTdz =  ((T(channel_pos(1), channel_j+1) - T(channel_pos(1), channel_j-1))/dz + ...
    (T(channel_pos(2), channel_j+1) - T(channel_pos(2), channel_j-1))/dz)/2;

%dSdz = dSdz/5;

if channel_side > 0
    Tchan = flip(Tchan);
    Slchan = flip(Slchan);
    psichan = flip(psichan);
    %xchan = fliplr(xchan);
end

xchan = xchan-xchan(1);



% Now get analytic soln

T0 = Tchan(1);
S0 = Slchan(1);
%dTdz = -5.0;
%dSdz = -0.4*dTdz;
doPlots = false;
[Tmodel,Smodel,psimodel,xmodel] = channelSolutionSin(Le,Ra,T0,S0,dSdz,dTdz,doPlots);


figure(f_i);

hold on;

plot(xchan,Tchan-T0, 'b-');
plot(xchan,Slchan-Smodel(end), 'r-');
plot(xchan,psichan-psichan(1), 'k-');

plot(xmodel,Tmodel-T0, 'b--');
plot(xmodel,Smodel-Smodel(end), 'r--');
plot(xmodel,psimodel-psimodel(1), 'k--');

plot([xchan(end) xchan(end)], [0 0.1], ':');

hold off;

xlabel('x');

title(['Ra = ',num2str(Ra),', Le = ',num2str(Le)]);

legend('T', 'Sl','psi', 'T model', 'Sl model', 'psi model', 'Location', 'eastoutside');

drawnow;

end