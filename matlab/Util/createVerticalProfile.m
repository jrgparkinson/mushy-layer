close all;
clear all;


z = linspace(0,2,100); %cm

H = NaN*z;
S = NaN*z;

p.St = 5;
p.Le = 200;
p.C = 1.15;
p.thetaInf = 1.001;
p.cp = 1;
p.pc = 1e-5;
p.thetaEutectic = 0;
p.ThetaEutectic = -p.thetaEutectic;

HTop = 2.5; %6.01;
Hbottom = 6.1;

Hmin = 2.0;
Hav = (Hbottom + Hmin)/2;

Tbottom = Hbottom - 5;
Tice = 1.01;
Tav = (Tbottom+Tice)/2;

blSize = 1;
blsizeChi = 0.2;

% ice depth
h =  z(end)/2;
z_h = min(find(z>h));
z_top = round(0.9*length(z));
htop = z(end) - z(z_top);
% 
T = Tav -  (Tbottom-Tav)*tanh((z-h)/blSize);
%T(z_top:end) = p.thetaInf - (p.thetaInf)*(z(z_top:end)-z(z_top))/htop;

Sice = -1.12; Sbottom = -1; Sav = (Sice+Sbottom)/2;
S = Sav + (Sice-Sav)*tanh((z-h)/blsizeChi);

chiMin = 0.2; chiAv = (1+chiMin)/2;
%chi = 1 + 0*z;
chi = chiAv  - (1-chiAv)*tanh((z-h)/blsizeChi);

%S = -1 + 0*z;
%S(z_h:end) = Sl(z_h:end).*chi(z_h:end);
H = p.St*chi + T;


%H(1:end) = p.St+p.thetaInf;
%H(z_h:end) = H(z_h) -3.5 ; %max(0.8,H(z_h) -3.5 - sqrt((z(z_h:end)-h))*0.005);
%H = Hav - (Hbottom-Hav)*tanh(z-h);
%H(z_top:end) = Hmin + (HTop-Hmin)* (z(z_top:end)-z(z_top))/(z(end)-z(z_top));
%
%S = -1+0*z;
figure();
hold on;
%plot(H, z, '-');
plot(S+2, z, '-');
plot(T, z, '-');
plot(chi, z, '-');

hold off
ylabel('z (cm)');

legend({'S+2', 'T', 'chi'}, 'location', 'eastoutside');


% Compute resulting temperature and porosity
[T,Sl, Ss, chi, HS, HE, HL] = computeEnthalpyVariables(H,S,p);

figure();
hold on;
plot(H, z, '-');
plot(HS, z, '-');
plot(HE, z, '-');
plot(HL, z, '-');

hold off
ylabel('z (cm)');

legend({'H','$H_S$', '$H_E$', '$H_L$'}, 'location', 'eastoutside');


figure();
hold on;
%plot(H, z, '-');
plot(S+2, z, '-');
plot(T, z, '-');
plot(chi, z, '-');
plot(Sl+2, z, '-');

hold off
ylabel('z (cm)');

legend({'S+2', 'T', 'chi', 'Sl+2'}, 'location', 'eastoutside');