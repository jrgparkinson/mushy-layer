%Load an output file
%close all;

set(groot, 'DefaultTextInterpreter', 'latex')
set(groot, 'DefaultLegendInterpreter', 'latex')



dim = 2;
frame = -1;
subcycled = true;

% typical value
dt = 1.0;

% Should be valid:
data_dir = ['/home/parkinsonjl/mnt/sharedStorage/AW002_PARKINSON_MUSH/', ...
    'optimalStates-highRes-newSoluteFlux/CR1.25RaC200Le200ChiCubedPermeabilitypts56-0/'];
output = MushyLayerOutput(dim, frame, data_dir, ...
    'CR1.25RaC200Le200ChiCubedPermeabilitypts56-274737', ...
    subcycled); %, 'bilinear')



data_dir = ['/home/parkinsonjl/convection-in-sea-ice/MushyLayer/', ...
    'run/refluxTest/Uniform-HRLSteady-16-0/'];
output = MushyLayerOutput(dim, frame, data_dir, ...
    'refluxTest-Uniform-HRLSteady-16-000759', ...
    subcycled); %, 'bilinear')


% Should not be valid:
% data_dir = ['/home/parkinsonjl/mnt/sharedStorage/AW002_PARKINSON_MUSH/', ...
%     'optimalStates-highRes-new/CR1.5RaC100Le200ChiCubedPermeabilitypts84-0/'];
% output = MushyLayerOutput(dim, frame, data_dir, ...
%     'CR1.5RaC100Le200ChiCubedPermeabilitypts84-358727', ...
%     subcycled); %, 'bilinear')

Le = 200;
V = 0;

%Get theta and velocity
T = output.dataForComp(output.components.Temperature);
H = output.dataForComp(output.components.Enthalpy);
Sl = output.dataForComp(output.components.Liquidconcentration);
chi = output.dataForComp(output.components.Porosity);
S = output.dataForComp(output.components.Bulkconcentration);
Ux = output.dataForComp(output.components.xAdvectionvelocity);
Uz = output.dataForComp(output.components.yAdvectionvelocity);
[X,Y] = output.grid();
dx = output.finest_dx();

X = X.'; Y = Y.';
x = X(:, 1);
y = Y(1, :);

H = H.';
S = S.';
T = T.';
Sl = Sl.';
Ux = Ux.';
Uz = Uz.';
chi = chi.';

% contourf(X, Y, theta);
[dH_dx, dH_dz] = gradient2order(H, dx, dx);
[dS_dx, dS_dz] = gradient2order(S, dx, dx);

[dtheta_dx, dtheta_dz] = gradient2order(T, dx, dx);
[dSl_dx, dSl_dz] = gradient2order(Sl, dx, dx);

T_advection = (Ux.*dtheta_dx + Uz.*dtheta_dz);
salt_advection = Ux.*dSl_dx + Uz.*dSl_dz;

lap_theta = 4*del2(T, dx);

[diff1, ~] = gradient2order(chi.*dSl_dx, dx, dx);
[~, diff2] = gradient2order(chi.*dSl_dz, dx, dx);

%[diff1, ~] = gradient2order(dSl_dx, dx, dx);
%[~, diff2] = gradient2order(dSl_dz, dx, dx);

salt_diffusion = (diff1 + diff2)/Le;

steady_state_heat = (lap_theta - V*dH_dz - T_advection)*dt;
steady_state_heat_half_advection = (lap_theta - 0.5*V*dH_dz - T_advection)*dt;
steady_state_salt = (salt_diffusion - V*dS_dz - salt_advection)*dt;




hFig = figure(1);
set(hFig, 'Position', [100 100 1400 800]);

m = 2; n = 2;
subplot(m, n, 1);
H = pcolor(x, y, steady_state_heat);
maxVal = max(max(abs(steady_state_heat)))/10;
set(H,'edgecolor','none');
caxis([-maxVal maxVal]);
 colormap(bluewhitered);
colorbar();
title('$\Delta t [\nabla^2 T - V dH/dz - U \cdot \nabla T]$');
%set(gca, 'YDir', 'normal');

subplot(m, n, 2);
H = pcolor(x, y, steady_state_salt);
set(H,'edgecolor','none');
maxVal = max(max(abs(steady_state_salt)));
caxis([-maxVal maxVal]);
colormap(bluewhitered);
colorbar();
title('$\Delta t [(1/Le) \nabla \cdot \chi \nabla S_l - V dS/dz - U \cdot \nabla S_l]$');
%set(gca, 'YDir', 'normal');

subplot(m, n, 3);
H = pcolor(x, y, steady_state_heat_half_advection);
maxVal = max(max(abs(steady_state_heat_half_advection)));
set(H,'edgecolor','none');
caxis([-maxVal maxVal]);
 colormap(bluewhitered);
colorbar();
title('$\Delta t [\nabla^2 T - V dH/dz - U \cdot \nabla T]$');
