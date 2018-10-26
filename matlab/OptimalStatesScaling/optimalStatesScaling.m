%clear all;
close all;

set(0, 'defaultlinelinewidth',2);
set(0, 'defaultaxeslinewidth',2);
set(0, 'defaultpatchlinewidth',2);
set(0, 'defaultAxesFontSize',22);

% Length scale (metres)
L = 2;

axisLabels = true;

baseDir = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/';

dataFolder = '/media/parkinsonjl/FREECOM HDD/mushyLayerLowC-upperBranch-insulating/';

files = []; legendStr = {};
y_wcrit = [];
y_chicrit = [];
y_psicrit = [];
chicrit = 0.4; wcrit = 0.8; psicrit = 0.5;
getData = true;
makeAllThePlots = true;
interiorPlots = false;


i = 0;

% List all optimal states here:
i = i+1;
files(i).output_dir = 'CR1.25RaC75Le200ChiCubedPermeabilitypts48-0/';
files(i).plot_prefix = 'mushyLayerLowC-upperBranch-insulating-CR1.25RaC75Le200ChiCubedPermeabilitypts48-';
files(i).frame = 59592;
files(i).legend = 'Ra = 75';
files(i).Ra = 75;

%
i = i+1;
files(i).output_dir = 'CR1.25RaC100Le200ChiCubedPermeabilitypts40-0/';
files(i).plot_prefix = 'mushyLayerLowC-upperBranch-insulating-CR1.25RaC100Le200ChiCubedPermeabilitypts40-';
files(i).frame = 50771;
files(i).legend = 'Ra = 100';
files(i).Ra = 100;

% Optimal state
i=i+1;
 files(i).output_dir = 'CR1.25RaC150Le200ChiCubedPermeabilitypts40-0/';
 files(i).plot_prefix = 'mushyLayerLowC-upperBranch-insulating-CR1.25RaC150Le200ChiCubedPermeabilitypts40-';
 files(i).frame = 45631;
 files(i).legend = 'Ra = 150';
 files(i).Ra = 150;

% Optimal state
i = i+1;
files(i).output_dir = 'CR1.25RaC200Le200ChiCubedPermeabilitypts32-0/';
files(i).plot_prefix = 'mushyLayerLowC-upperBranch-insulating-CR1.25RaC200Le200ChiCubedPermeabilitypts32-';
files(i).frame = 43334;
files(i).legend = 'Ra = 200';
files(i).Ra = 200;

% Optimal state
i = i+1;
files(i).output_dir = 'CR1.25RaC300Le200ChiCubedPermeabilitypts28-0/';
files(i).plot_prefix = 'mushyLayerLowC-upperBranch-insulating-CR1.25RaC300Le200ChiCubedPermeabilitypts28-';
files(i).frame = 44705;
files(i).legend = 'Ra = 300';
files(i).Ra = 300;

% Optimal state
i = i+1;
files(i).output_dir = 'CR1.25RaC400Le200ChiCubedPermeabilitypts28-0/';
files(i).plot_prefix = 'mushyLayerLowC-upperBranch-insulating-CR1.25RaC400Le200ChiCubedPermeabilitypts28-';
files(i).frame = 141902;
files(i).legend = 'Ra = 400';
files(i).Ra = 400;

% files(2).output_dir = 'CR1.25RaC250Le200ChiCubedPermeabilitypts32-0/';
% files(2).plot_prefix = 'mushyLayerLowC-insulating-CR1.25RaC250Le200ChiCubedPermeabilitypts32-';
% files(2).frame = 37322;
% files(2).legend = 'Ra = 250';

% Optimal state
i = i+1;
  files(i).output_dir = 'CR1.25RaC500Le200ChiCubedPermeabilitypts28-0/';
  files(i).plot_prefix = 'mushyLayerLowC-upperBranch-insulating-CR1.25RaC500Le200ChiCubedPermeabilitypts28-';
  files(i).frame = 142094;
  files(i).legend = 'Ra = 500';
  files(i).Ra = 500;



% Optimal state
i = i+1;
files(i).output_dir = 'CR1.25RaC600Le200ChiCubedPermeabilitypts56-0/';
files(i).plot_prefix = 'mushyLayerLowC-upperBranch-insulating-CR1.25RaC600Le200ChiCubedPermeabilitypts56-';
files(i).frame = 082852;
files(i).legend = 'Ra = 600';
files(i).Ra = 600;


%Optimal state
 i = i+1;
   files(i).output_dir = 'CR1.25RaC700Le200ChiCubedPermeabilitypts48-0/';
   files(i).plot_prefix = 'mushyLayerLowC-upperBranch-insulating-CR1.25RaC700Le200ChiCubedPermeabilitypts48-';
   files(i).frame = 117900;
   files(i).legend = 'Ra = 700';
   files(i).Ra = 700;

%optimal state
% i = i+1;
% files(i).output_dir = 'CR1.25RaC800Le200ChiCubedPermeabilitypts48-0/';
% files(i).plot_prefix = 'mushyLayerLowC-upperBranch-insulating-CR1.25RaC800Le200ChiCubedPermeabilitypts48-';
% files(i).frame = 122778;
% files(i).legend = 'Ra = 800';
% files(i).Ra = 800;

i = i+1;
files(i).output_dir = 'CR1.25RaC800Le200ChiCubedPermeabilitypts48-0/';
files(i).plot_prefix = 'mushyLayerLowC-upperBranch-insulating-CR1.25RaC800Le200ChiCubedPermeabilitypts48-';
files(i).frame = 160877;
files(i).legend = 'Ra = 800';
files(i).Ra = 800;


%colors = {'b', 'r', 'k', 'm', 'g', 'c', 'y', 'b', 'r', 'k'};
%colors =  linspecer(length(files)+5);
colors = [   0    0.4470    0.7410
    0.8500    0.3250    0.0980
    
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840
    0.1     0.1     0.1
    0.9290    0.6940    0.1250];
%axes('NextPlot','replacechildren', 'ColorOrder',C);

if getData
    
    RaCs = NaN*ones(1, length(files));
fluxes  = NaN*RaCs;
widths = NaN*fluxes;
channelDepths = NaN*widths;
maxMushyU = NaN*widths;
maxLiqU = NaN*widths;
maxLiqV = NaN*widths;
maxMushyV = NaN*widths;
maxSpeed = NaN*widths;
    plots = [];
logPsiMaxVec = NaN*ones(length(files), 128);

channelDepthPorosity = NaN*widths;
averagePorosity = NaN*widths;
channelDepthPerm = NaN*widths;
averagePerm = NaN*widths;

maxMushySaltAdv = NaN*widths;

psiMax = NaN*widths;
psiMaxMush = NaN*widths;

halfPorosityDepth = NaN*widths;
mushyHeights = NaN*widths;

initial_file_i = 1; % should be 1 to do all files

for file_i = initial_file_i:length(files)
    
    
    %If possible, get optimal state properties by interpolating from the
    %three states near the optimal state.
    
    
    pout = Pout([dataFolder, files(file_i).output_dir, '/pout.0']);
    inputs = readInputs([dataFolder, files(file_i).output_dir, '/inputs']);
    
    widths(file_i) = str2double(inputs.domain_length);
    RaCs(file_i) = str2double(inputs.rayleighComp);
    fluxes(file_i) = pout.fluxBottom(end);
    
    dim = 2; subcycled = true;
    
    output = MushyLayerOutput(dim, files(file_i).frame, [dataFolder, files(file_i).output_dir],...
        files(file_i).plot_prefix, subcycled);
    
    
    
    [chanWidth, chanDepth, numChannels, mushyHeights(file_i)] = output.channelGeometry();
    [maxMushyU(file_i), maxMushyV(file_i), maxSpeed(file_i)] = output.maxMushyVel();
     [maxLiqU(file_i), maxLiqV(file_i), ~] = output.maxLiqVel();
    [channelDepthPorosity(file_i), averagePorosity(file_i),...
        halfPorosityDepth(file_i)] = output.porosityMetrics(-0.5);
    [channelDepthPerm(file_i), averagePerm(file_i)] = output.permMetrics();
    
    maxMushySaltAdv(file_i) = output.maxMushyAdvectionSrcTerm();

    channelDepths(file_i) = chanDepth;
    
    if numChannels > 0
        widths(file_i) = widths(file_i)/numChannels;
    end
    
%     zPorosity = output.horizontallyAverage(output.components.Porosity, false);
%     V_averaged = output.horizontallyAverage(output.components.yAdvectionvelocity, true);
%     
%     Vfield = output.dataForComp(output.components.yAdvectionvelocity);
%     V_max = output.horizontallyMax(Vfield, true);
%     
     porosity = output.dataForComp(output.components.Porosity);
     psi = output.getStreamfunction(2000, 1).';
     psiMax(file_i) = max(max(abs(psi)));
     
     psiMush = psi;
     psiMush(porosity > 0.98) = NaN;
     psiMush(porosity < 0.001) = NaN;
     psiMaxMush(file_i) = max(max(abs(psiMush)));
     
     
     if makeAllThePlots
         
     %Calculate each term in the governing steady state equations and see
     %how they scale (in the mushy layer)
     H = output.dataForComp(output.components.Enthalpy);
     U = output.dataForComp(output.components.xAdvectionvelocity);
     V = output.dataForComp(output.components.yAdvectionvelocity);
     T = output.dataForComp(output.components.Temperature);
     S = output.dataForComp(output.components.Bulkconcentration);
     Sl = output.dataForComp(output.components.Liquidconcentration);
     Ss = output.dataForComp(output.components.Solidconcentration);
     Permeability = output.dataForComp(output.components.Permeability);
     psi = output.getStreamfunction();
     [X,Z] = output.grid();
     dx = output.problemDomain.dxCoarse;
     Le = str2double(inputs.lewis);
     Ra = RaCs(file_i);
     concratio = str2double(inputs.compositionRatio);
     stefan = str2double(inputs.stefan);
     
     % All future operations expect these dimensions to be the other way
     % around
     
     % These two possibly alright now
     %X = X.'; Z = Z.';
     %vorticity = vorticity.';
     
     porosity = porosity.';
     H = H.';
     T = T.';
     Sl = Sl.';
     Ss = Ss.';
     Permeability = Permeability.';
     U = U.'; V = V.';
     S = S.';
     
     %Convert to Wells dimensionless units
     Ss = Ss + 1;
     Sl = Sl + 1;
     
     [dChidx, dChidz] = gradient(porosity, dx);
     [~, dthetaChidz] = gradient(porosity.*T, dx);
     
     %Heat equation:
     [~, dHdz] = gradient(H, dx);
     lapT = 4*del2(T, dx);
     [dTdx, dTdz] = gradient(T, dx);
     heatAdvection = U.*dTdx + V.*dTdz;
     
     %Salt equation:
     [~, dSdz] = gradient(S, dx);
     [dSldx, dSldz] = gradient(Sl, dx);
     [~, dChiSldz] = gradient(Sl.*porosity, dx);
     [~, dChiSsdz] = gradient(Ss.*(1-porosity), dx);
     saltAdvection = U.*dSldx + V.*dSldz;
     saltDiffusion = (1/Le)*divergence(X, Z, porosity.*dSldx, porosity.*dSldz);
     
     %Momentum:
     vorticityDissipation = 4*del2(psi, dx);
     baroclinicTorque = Ra*Permeability.*dSldx;
     [dPsidx, dPsidz] = gradient(psi, dx);
     [dPermdx, dPermdz] = gradient(Permeability, dx);
     lastTerm = (dPsidx.*dPermdx + dPsidz.*dPermdz)./Permeability;
     
     % Work out approximate mushylayer 
     %isMush = porosity;
     minPorosity = min(min(porosity));
     eutecticPorosity = (1-1/concratio);
     %eutecticPorosity = 0.0;
     if interiorPlots
     maxPoros = 0.4; % 0.97
     else
         maxPoros = 0.97; 
     end
     isMush = (porosity < maxPoros) + (porosity > 0.01) - 1;
     [y_i, x_i] = find(isMush == 1);
     %isMush(max(y_i), :) = 0; % Trying to remove solid-mush boundary
     %isMush(max(y_i)-1, :) = 0; % Trying to remove solid-mush boundary
     % isMush(max(y_i)-2, :) = 0; % Trying to remove solid-mush boundary
      % isMush(max(y_i)-3, :) = 0; % Trying to remove solid-mush boundary
      
      y_max = max(y_i)-2; % Get rid of solid-mush boundary layer
      y_min = min(y_i);
       
       %Remove left hand channel to make comparison easier
       x_min = 1;
       [yLength, xLength]  = size(isMush);
       if numChannels > 1
           
          x_min = round(xLength/2);
       end
       
       x_max = xLength;
       if interiorPlots
        x_max = xLength - 5; % Get rid of chimney - interior plots
       end
       
     dHdz(~isMush) = NaN;                dHdz = dHdz(y_min:y_max, x_min:x_max);
     heatAdvection(~isMush) = NaN;       heatAdvection = heatAdvection(y_min:y_max, x_min:x_max);
     lapT(~isMush) = NaN;                lapT = lapT(y_min:y_max, x_min:x_max);
     
     dSdz(~isMush) = NaN;                   dSdz = dSdz(y_min:y_max, x_min:x_max);
     saltAdvection(~isMush) = NaN;          saltAdvection = saltAdvection(y_min:y_max, x_min:x_max);
     saltDiffusion(~isMush) = NaN;          saltDiffusion = saltDiffusion(y_min:y_max, x_min:x_max);
     dChiSldz(~isMush) = NaN;               dChiSldz = dChiSldz(y_min:y_max, x_min:x_max);
     dChiSsdz(~isMush) = NaN;               dChiSsdz = dChiSsdz(y_min:y_max, x_min:x_max);
     
     vorticityDissipation(~isMush) = NaN;   vorticityDissipation =vorticityDissipation(y_min:y_max, x_min:x_max);
     baroclinicTorque(~isMush) = NaN;       baroclinicTorque = baroclinicTorque(y_min:y_max, x_min:x_max);
     lastTerm(~isMush) = NaN;               lastTerm = lastTerm(y_min:y_max, x_min:x_max);
     
     psi(~isMush) = NaN;                    psi = psi(y_min:y_max, x_min:x_max);
     U(~isMush) = NaN;                      U = U(y_min:y_max, x_min:x_max);
     V(~isMush) = NaN;                      V = V(y_min:y_max, x_min:x_max);
     Sl(~isMush) = NaN;                     Sl = Sl(y_min:y_max, x_min:x_max);
     T(~isMush) = NaN;                     T = T(y_min:y_max, x_min:x_max);
     porosity(~isMush) = NaN;                    porosity = porosity(y_min:y_max, x_min:x_max);
     
     dChidz(~isMush) = NaN;                 dChidz = dChidz(y_min:y_max, x_min:x_max);
     dTdz(~isMush) = NaN;               dTdz = dTdz(y_min:y_max, x_min:x_max);
     dthetaChidz(~isMush) = NaN;        dthetaChidz = dthetaChidz(y_min:y_max, x_min:x_max);
     dSldz(~isMush) = NaN;                  dSldz = dSldz(y_min:y_max, x_min:x_max);
     
     
     X = X(y_min:y_max, x_min:x_max);
     Z = Z(y_min:y_max, x_min:x_max);
     
     
     x = X(1, :);
     z = Z(:, 1);
     
     %saltDiffusion = 0.0*saltAdvection;
     
     hPoros = figure(1);
      set(hPoros, 'Position', [200 200 800 800]);
     pcolor(x, z, porosity); colorbar();
      title('porosity');
     set(gca,'YDir','normal');      axis equal; shading flat; 
     
     
     h = figure(2);
     set(h, 'Position', [50 50 1600 1200]);
     m = 3;  n=4;
     
     colormap jet;
     
     
     plot_i = 1;
     subplot(m, n, plot_i);
     pcolor(x, z, abs(dTdz)); colorbar();
      title('dT/dz');
     set(gca,'YDir','normal');      axis equal; shading flat; 
     
     plot_i = plot_i + 1;
     subplot(m, n, plot_i);
     pcolor(x, z, stefan*abs(dChidz)); colorbar();
      title('St*d(\chi)/dz');
     set(gca,'YDir','normal');      axis equal; shading flat; 
  
     
     plot_i = plot_i + 1;
     subplot(m, n, plot_i);
     pcolor(x, z, abs(heatAdvection)); colorbar();
      title('-U.grad(T)');
     set(gca,'YDir','normal');      axis equal; shading flat; 
     
     plot_i = plot_i + 1;
     subplot(m, n, plot_i);
     pcolor(x, z, abs(lapT)); colorbar();
     title('grad^2 T');
     set(gca,'YDir','normal');      axis equal; shading flat; 
     
     
     if n>4
         plot_i = plot_i + 1;
         subplot(m, n, plot_i);
         pcolor(x, z, abs(dHdz  + heatAdvection - lapT)); colorbar();
         title('sum');
         set(gca,'YDir','normal');      axis equal; shading flat; 
     end
     
%      plot_i = plot_i + 1;
%      subplot(m, n, plot_i);
%      pcolor(x, z, abs(dChiSldz)); colorbar();
%      title('d(\chi S_l)/dz');
%      set(gca,'YDir','normal');      axis equal; shading flat; 
 plot_i = plot_i + 1;
     subplot(m, n, plot_i);
     pcolor(x, z, abs(Sl.*dChidz)); colorbar();
     title('S_l d(\chi)/dz');
     set(gca,'YDir','normal');      axis equal; shading flat; 
     
      plot_i = plot_i + 1;
     subplot(m, n, plot_i);
     pcolor(x, z, abs(porosity.*dSldz)); colorbar();
     title('\chi d(S_l)/dz');
     set(gca,'YDir','normal');      axis equal; shading flat; 
     
      plot_i = plot_i + 1;
     subplot(m, n, plot_i);
     pcolor(x, z, (dChiSldz)); colorbar();
     title('d(S_l \chi)/dz');
     set(gca,'YDir','normal');      axis equal; shading flat; 
     
     plot_i = plot_i + 1;
     subplot(m, n, plot_i);
     pcolor(x, z, (dChiSsdz)); colorbar();
     title('d((1-\chi) S_s)/dz');
     set(gca,'YDir','normal');      axis equal; shading flat; 
     
%      subplot(m, n, plot_i);
%      pcolor(x, z, abs(saltAdvection)); colorbar();
%      title('U.grad(S_l)');
%      set(gca,'YDir','normal');      axis equal; shading flat; 
     
plot_i = plot_i + 1;
     subplot(m, n, plot_i);
     pcolor(x, z, dSdz); colorbar();
     title('dS/dz');
     set(gca,'YDir','normal');      axis equal; shading flat; 
     
     
%      plot_i = plot_i + 1;
%      subplot(m, n, plot_i);
%      pcolor(x, z, saltDiffusion); colorbar();
%      title('(1/Le)div(chi grad(S_l))');
%      set(gca,'YDir','normal');      axis equal; shading flat; 
     
     if n>4
     plot_i = plot_i + 1;
     subplot(m, n, plot_i);
     pcolor(x, z, abs(dSdz + saltAdvection + saltDiffusion)); colorbar();
     title('sum');
     set(gca,'YDir','normal');      axis equal; shading flat; 
     end
     
     
     % Blank plot
%      plot_i = plot_i + 1;
%          subplot(m, n, plot_i);
         
     
     plot_i = plot_i + 1;
     subplot(m, n, plot_i);
     
     pcolor(x, z, abs(vorticityDissipation)); colorbar();
     title('grad^2 \psi');
     set(gca,'YDir','normal');      axis equal; shading flat; 
     
     plot_i = plot_i + 1;
     subplot(m, n, plot_i);
     
      pcolor(x, z, abs(baroclinicTorque)); colorbar();
      title('Ra \Pi dT/dx');
      set(gca,'YDir','normal');      axis equal; shading flat; 
      
      plot_i = plot_i + 1;
     subplot(m, n, plot_i);
     
      pcolor(x, z, abs(lastTerm)); colorbar();
      title('(grad \psi . grad \Pi) / \Pi');
    set(gca,'YDir','normal');      axis equal; shading flat; 
    
    if n> 4
    plot_i = plot_i + 1;
    subplot(m, n, plot_i);
      pcolor(x, z, abs(vorticityDissipation - baroclinicTorque - lastTerm)); colorbar();
      title('sum');
    set(gca,'YDir','normal');      axis equal; shading flat; 
    end
    
    interiorPlotsStr = '';
    if interiorPlots
        interiorPlotsStr = 'Interior';
    end
    
     print(h,['equationScaling/', interiorPlotsStr, 'Rac', num2str(Ra),'equationScaling.png'],'-dpng','-r150')
    
     
     hFields = figure(3);
     colormap jet;
     set(hFields, 'Position', [50 50 1600 1200]);
     
     m = 2; n = 3;
     plot_i = 1;
     subplot(m, n, plot_i);
     pcolor(x, z, abs(porosity)); colorbar();
      title('\chi');
    set(gca,'YDir','normal');      axis equal; shading flat; 
    
    plot_i = plot_i+1;
    subplot(m, n, plot_i);
     pcolor(x, z, Sl); colorbar();
      title('S_l');
    set(gca,'YDir','normal');      axis equal; shading flat; 
    
     plot_i = plot_i+1;
    subplot(m, n, plot_i);
     pcolor(x, z, abs(psi)); colorbar();
      title('\psi');
    set(gca,'YDir','normal');      axis equal; shading flat; 
     
     plot_i = plot_i+1;
    subplot(m, n, plot_i);
     pcolor(x, z, abs(U)); colorbar();
      title('U');
    set(gca,'YDir','normal');      axis equal; shading flat; 
    
    plot_i = plot_i+1;
    subplot(m, n, plot_i);
     pcolor(x, z, abs(V)); colorbar();
      title('V');
    set(gca,'YDir','normal');      axis equal; shading flat; 
    
       print(hFields,['fields/', interiorPlotsStr, 'Rac', num2str(Ra),'fields.png'],'-dpng','-r150')
    
    
     % pause;
      
     end
      
     %psiMaxMush(file_i) = nanmedian(nanmedian(abs(psiMush)));
     %psiMaxMush(file_i) = nanmean(nanmean(abs(psiMush)));
     
%     psi_max = output.horizontallyMax(psi, true);
%     
%     probDomain = output.problemDomain;
%     dx = probDomain.dxCoarse;
%     
%     numx = (probDomain.domainExtent.hi_i - probDomain.domainExtent.lo_i) + 1;
%     numy = (probDomain.domainExtent.hi_j - probDomain.domainExtent.lo_j) + 1;
%     
%     x = linspace(double(probDomain.domainExtent.lo_i), double(probDomain.domainExtent.hi_i), numx)*dx;
%     y = linspace(double(probDomain.domainExtent.lo_j), double(probDomain.domainExtent.hi_j), numy)*dx;
%     
%     
%     xlo = double(probDomain.domainExtent.lo_i)*dx;
%     xhi = double(probDomain.domainExtent.hi_i)*dx;
%     width = xhi-xlo;
    
    %[X, Y] = meshgrid(x, y);
    
%     y_wcrit(file_i) = NaN;
%     y_chicrit(file_i)  = NaN;
%     y_psicrit(file_i)= NaN;
%     for j = 1:numy
%         
%         if V_max(j) < wcrit && isnan(y_wcrit(file_i))
%             y_wcrit(file_i) = y(j);
%         end
%         
%         if zPorosity(j) < chicrit && isnan(y_chicrit(file_i))
%             y_chicrit(file_i) =  y(j);
%         end
%         
%         if psi_max(j) < psicrit && isnan(y_psicrit(file_i)) && y(j) > 0.05
%             y_psicrit(file_i) = y(j);
%         end
%     end
%     
%    
    
    legendStr{file_i} = files(file_i).legend;
    
end

end % end if get data

h = figure();
h.Position = [100 100 1000 1000];
m = 3; %rows
n = 3; % cols
subplot(m, n, 1)
plot(log(RaCs), fluxes-1, 'x');
xlabel('log(Ra_c)'); ylabel('F-V'); 
 
subplot(m, n, 2)
plot(log(RaCs).^0.5, 1./widths, 'x');
xlabel('log(Ra_c)^{0.5}'); ylabel('1/\lambda_0'); 

subplot(m, n, 3)
%plot(log(RaCs).^0.5, channelDepths, 'x');
plot(1./log(1./RaCs), channelDepths, 'x');
xlabel('-1/log(Ra_c)'); ylabel('h'); 

subplot(m, n, 4)
plot(RaCs, widths./channelDepths, 'x');
xlabel('Ra_C'); ylabel('\lambda_0 / h');


% subplot(m, n, 6)
% plot(RaCs, widths.*channelDepths, 'x');
% xlabel('Ra_c'); ylabel('h \lambda_0'); 

subplot(m, n, 5)
plot(RaCs.^(1/3) .* log(RaCs), maxMushyU, 'x');
xlabel('Ra_c^{1/3} log(Ra_c)'); ylabel('u_{max}'); 

subplot(m, n, 6)
plot((RaCs.^(4/3)), maxMushyV, 'x');
xlabel('Ra_c^{4/3}'); ylabel('v_{max}'); 

subplot(m, n, 7)
plot(RaCs./log(RaCs), maxMushyV./maxMushyU, 'x');
xlabel('Ra_C/log(Ra_c)'); ylabel('v_{max}/u_{max}');


subplot(m, n, 8)
plot(RaCs.^0.3 .* log(RaCs), maxSpeed, 'x');
xlabel('Ra_C^{1/3} log(Ra_C)'); ylabel('|u|_{max}');

subplot(m, n, 9)
plot(RaCs.^0.3, maxLiqU, 'x');
xlabel('Ra_C^{1/3}'); ylabel('u_{max, liquid}');


% subplot(m, n, 9)
% plot(log(RaCs), channelDepths./widths, 'x');
% xlabel('log(Ra_C)'); ylabel('h / \lambda_0');


% Need another figure as there's so many plots
m = 3; n = 3;
h2 = figure();
h2.Position = [400 100 1000 1000];
subplot(m, n, 1);
plot(RaCs.^(-1/2), averagePerm, 'x'); xlabel('1/Ra_c^{1/2}'); ylabel('<\Pi> ');

%subplot(m, n, 2);
%plot(RaCs.^(-1/2), averagePerm, 'x'); xlabel('1/Ra_c^{1/2}'); ylabel('<\Pi>');
subplot(m, n, 2);
plot((1/stefan)*(1+averagePerm.^(1/3)), widths, 'x'); xlabel('St^{-1}(1+<\Pi>^{1/3})'); ylabel('\lambda');


subplot(m, n, 3);
plot(averagePorosity.^3, averagePerm, 'x'); 
%xlabel('log(1/Ra_C)^3'); 
xlabel('\chi^3'); 
ylabel('<\Pi>');

subplot(m, n, 4);
plot((1/stefan)*(1+averagePerm.^(1/3)), mushyHeights, 'x'); 
%xlabel('log(1/Ra_C)^3'); 
xlabel('St^{-1}(1+<\Pi>^{1/3})'); 
ylabel('H_{mush}');

subplot(m, n, 5);
plot(RaCs, psiMaxMush./(widths.*averagePerm.^2), 'x'); 
xlabel('Ra_c'); 
ylabel('\psi_{mush}/(\lambda <\Pi>^2');

subplot(m, n, 6);
plot(1-averagePorosity, psiMaxMush./widths, 'x'); 
xlabel('1-<\chi>');
ylabel('\psi_{mush}/\lambda');

%subplot(m, n, 6);
%plot(log(RaCs).*(RaCs.^(-1/3)), halfPorosityDepth, 'x'); 
%xlabel('log(1/Ra_C)^3'); 
%subplot(m, n, 6);
%xlabel('log(Ra_C)/Ra_C^{1/3}'); 
%plot(RaCs, halfPorosityDepth, 'x'); 
%xlabel('Ra_C');
%ylabel('y (\chi = <\chi>) ');


%subplot(m, n, 7)
%plot(RaCs, (psiMax.*averagePerm)./(psiMaxMush), 'x');
%xlabel('Ra_C'); ylabel('\Pi \psi_{max liq} / (\psi_{max mush})', 'interpreter','tex');

subplot(m, n, 7)
%plot((1-1./RaCs.^(1/6)), fluxes-1, 'x');
plot((1-averagePorosity).*(averagePorosity + concratio*(1-averagePorosity)), fluxes-1, 'x');
xlabel('\Delta \chi (\chi + C \Delta \chi)'); ylabel('F-V', 'interpreter','tex');


subplot(m, n, 8)
plot((1-1./RaCs.^(1/6)), fluxes-1, 'x');
xlabel('1-1/Ra_c^{1/6}'); ylabel('F-V', 'interpreter','tex');

subplot(m, n, 9)
plot(RaCs.^(1/6), fluxes-1, 'x');
xlabel('Ra_c^{1/6}'); ylabel('F-V', 'interpreter','tex');



%subplot(m, n, 8)
%plot(averagePerm, widths, 'x');
%xlabel('<\Pi>'); ylabel('\lambda');


%subplot(m, n, 9)
%plot(RaCs, psiMaxMush, 'x');
%xlabel('Ra_C'); ylabel('\psi_{max} (mush)');


% subplot(m, n, 9)
% plot(averagePerm.^2, widths./channelDepths, 'x');
% xlabel('<\Pi>^2'); ylabel('\lambda / H');



iceDepth = 0.1;
%scaledRa = RaCs.*exp(-iceDepth./(channelDepthPerm));
scaledRa = RaCs.*averagePerm.*channelDepths.^2;

h = figure();
h.Position = [100 100 1000 1000];
m = 3; %rows
n = 3; % cols
subplot(m, n, 1)
plot(scaledRa, fluxes-1, 'x');
xlabel('scaled Ra_c'); ylabel('F-V'); 
 
subplot(m, n, 2)
plot(scaledRa, 1./widths, 'x');
xlabel('scaled Ra_c'); ylabel('1/\lambda_0'); 

subplot(m, n, 3)
plot(scaledRa, channelDepths, 'x');
xlabel('scaled Ra_c'); ylabel('h'); 

subplot(m, n, 4)
plot(scaledRa, widths./channelDepths, 'x');
xlabel('scaled Ra_c'); ylabel('\lambda_0 / h');


% subplot(m, n, 6)
% plot(RaCs, widths.*channelDepths, 'x');
% xlabel('Ra_c'); ylabel('h \lambda_0'); 

% subplot(m, n, 5)
% plot(RaCs.^(1/3) .* log(RaCs), maxMushyU, 'x');
% xlabel('Ra_c^{1/3} log(Ra_c)'); ylabel('u_{max}'); 
% 
% subplot(m, n, 6)
% plot((RaCs.^(4/3)), maxMushyV, 'x');
% xlabel('Ra_c^{4/3}'); ylabel('v_{max}'); 
% 
% subplot(m, n, 7)
% plot(RaCs./log(RaCs), maxMushyV./maxMushyU, 'x');
% xlabel('Ra_C/log(Ra_c)'); ylabel('v_{max}/u_{max}');
% 

% subplot(m, n, 8)
% plot(RaCs.^0.3 .* log(RaCs), maxSpeed, 'x');
% xlabel('Ra_C^{1/3} log(Ra_C)'); ylabel('|u|_{max}');
% 
% subplot(m, n, 9)
% plot(RaCs.^0.3, maxLiqU, 'x');
% xlabel('Ra_C^{1/3}'); ylabel('u_{max, liquid}');


% subplot(m, n, 9)
% plot(log(RaCs), channelDepths./widths, 'x');
% xlabel('log(Ra_C)'); ylabel('h / \lambda_0');

%
% ylabel('y');
% xlabel('Ra');


 