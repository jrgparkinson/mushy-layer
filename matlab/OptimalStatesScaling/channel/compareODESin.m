% Compare ODE and sinusoidal solutions
% It is NOT possible to predict features that vary along the height of the
% channel with this method!

close all;

Ra = 100; %60;
Le = 200;
St = 5;

H = 1/St;
dTdz = -St;
dSdz = -0.3*dTdz;

h = H/2;
%z = h/2; % z = 0 is bottom of channel, z=h is top
z = 0:0.001:h;


TTop = 0.5;
TBottom = 0.9;
dTdz = (TTop-TBottom)/h + 0*z;

STop = -TTop;
SBottom = STop-0.1;
dSdz = (STop-SBottom)/h + 0*z;

% Need dTdz(h) = Le*dSdz;
%dSdz = -dTdz/Le + ((STop-SBottom)/h)*(h-z);


T0 = TBottom+z.*dTdz;
S0 = SBottom+z.*dSdz;

% These should actually be at least quadratic I think
%S0 = SBottom + dSdz.*(z);
%T0 = TBottom + dTdz.*(z);

z_i = round(length(z)/2);

%dSdz = 5.0; dTdz = -4.0;
%T0 = 0.738; S0 = T0 + abs(0.565-T0);

doPlots = true;
%ainit = 0.0113; stepSize = 5e-5;
%[Tode,Sode,psiode,areturn,xode] = channelSolutionODE(Le,Ra,T0,S0,dSdz,dTdz,doPlots,ainit, stepSize);


a=channelWidth(Le,Ra,T0(z_i),S0(z_i),dSdz(z_i),dTdz(z_i))

[Tsin,Ssin,psisin,xsin] = channelSolutionSin(Le,Ra,T0(z_i),S0(z_i),dSdz(z_i),dTdz(z_i),doPlots,a);

close all;

figure(1);
plot(z, T0, '-', z, S0, '-', z,T0+S0,'-');
xlabel('z');
legend('T0', 'S0','T0+S0');


hFig = figure(2);
set(hFig, 'Position', [200 200 1000 900]);

%m = 1; n = 2;

%subplot(m,n, 1);

hold on;

%plot(xode, Tode,'r-');
%plot(xode,Sode,'g-');
%plot(xode,psiode,'b-');

plot(xsin,Tsin,'r--');
plot(xsin,Ssin,'g--');
plot(xsin,psisin,'b--');


hold off

xlim([0 a]);

legend('T','-S ','psi', 'Location', 'northoutside'); %'Tode','Sode','psiode',


% subplot(m, n, 2);
% 
% 
% a = T0*NaN;
% 
% for i=1:length(T0)    
%     a(i) = channelWidth(Le,Ra,T0(i),S0(i),dSdz(i),dTdz(i));
% end
% 
% hold on;
% plot(a,z, '-');
% hold off;
% xlabel('a');
% ylabel('z');

%daspect([1 10 1]);





% Le = logspace(0,10,100);
% a = Le*NaN;
% 
% for i=1:length(Le)
%     a(i) = channelWidth(Le(i),Ra,T0,S0,dSdz,dTdz);
% end
% 
% figure();
% 
% xPlot = log10(1./Le);
% yPlot = log10(a);
% 
% a_half = xPlot*0.5-1;
% 
% hold on;
% plot(xPlot, yPlot, '-');
% plot(xPlot, a_half, '--');
% 
% hold off;
% 
% xlabel('log10(1/Le)'); ylabel('a');


%Le = logspace(0,10,100);
