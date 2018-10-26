% Plot optimal values

clear all;
close all;

output_folder = 'optimalStates-Brinkman/';
output_folder = 'optimalStates-highRes-new/';

%set(groot, 'DefaultAxesColorOrder', [1 0 0; 0 1 0; 0 0 1], ...
%    'DefaultAxesLineStyleOrder', 'x-|x--|x:');
set(groot, 'DefaultAxesLineStyleOrder', 'x-|x--|x:');
%set(groot, 'DefaultAxesLineStyleOrder',{'--',':'})

data_dir = getDataDir(output_folder);

figure_output_dir = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/matlab/OptimalStatesScaling/';

% load data:
load([data_dir, 'optimalVals.mat']);

printQuality = '-r100';
makeStandardPlots = false;
makeOceanPlots  = false;
makeDetailedPlots = false;
makeStabilityPlots = true;
makeFluxDecompositionPlots = false;

scaleDetailedPlotsWithRmCrit = true;

if strcmp(output_folder,'optimalStates-highRes-new/')
    RM_MAX = 600;
    CR_MAX = 20.0;
    Da = 1.0;
else
    RM_MAX = 1e5;
    CR_MAX = 40.0;
    Da = 5e-3;
end


% Get Rayleigh numbers and concentration ratios

% Get all concentration ratios and rayleigh numbers
k = keys(optimalVals);
CR_arr = [];
Rm_arr = [];
for CR_i = 1:length(k)
    keyPair = k(CR_i);
    keyPair = keyPair{1};
    temp = keyPair(1);  thisC  = temp{1};
    temp = keyPair(2);  thisRa = temp{1};
    
    
    if  ~ismember(thisC, CR_arr)
        CR_arr(end+1) = thisC;
    end
    
    if ~ismember(thisRa, Rm_arr)
        Rm_arr(end+1) = thisRa;
    end
    
end

Rm_arr = sort(Rm_arr); CR_arr = sort(CR_arr);

%Stefan number
St = 5;
Le = 200;
Pr = 10;


% Now, for each Ra, C pair, get matrices of optimal values

optimalFluxMat = NaN*ones(length(CR_arr), length(Rm_arr));
optimalWidthMat = optimalFluxMat;
optimalCritWidthMat= optimalFluxMat;
optimalPorosMat = optimalFluxMat;
optimalPermMat = optimalFluxMat;
optimalPsiMat = optimalFluxMat;

optimalL2DiffusionFluxMat = optimalFluxMat;
optimalL2FrameFluxMat = optimalFluxMat;
optimalL2FluidFluxMat = optimalFluxMat;

optimalL0DiffusionFluxMat = optimalFluxMat;
optimalL0FrameFluxMat = optimalFluxMat;
optimalL0FluidFluxMat = optimalFluxMat;

optimallMat = optimalFluxMat;
optimalHMat = optimalFluxMat;
optimalhMat = optimalFluxMat;
optimalhPorosityMat = optimalFluxMat;
optimalChannelHMat = optimalFluxMat;

averagePorosMat = optimalFluxMat;
averagePermMat = optimalFluxMat;

optimalHeatAdvection = optimalFluxMat;
optimalHeatDiffusion = optimalFluxMat;
optimalLatentHeat = optimalFluxMat;
optimalTempFrame = optimalFluxMat;

optimalSaltAdvection = optimalFluxMat;
optimalSaltDiffusion = optimalFluxMat;
optimalSolidSaltFrame = optimalFluxMat;
optimalLiquidSaltFrame = optimalFluxMat;

optimalVorticityDissipation = optimalFluxMat;
optimalBaroclinicTorque = optimalFluxMat;
optimalVorticityPermeability = optimalFluxMat;

optimalExtrapolationFlag = optimalFluxMat;
optimalExtrapolationFlag(isnan(optimalExtrapolationFlag)) = false;
optimalFlagCritWidth = optimalExtrapolationFlag;



for CR_i = 1:length(k)
    key = k(CR_i); key = key{1};
    CR = key{1}; Ra = key{2};
    state = optimalVals(CR, Ra);
    
    fprintf('CR: %1.2f, Ra: %1.3f, flux: %1.5f \n', CR, Ra, state.flux);
    
    
    if strcmp(output_folder,'optimalStates-Brinkman/')
        Rm = Da*Ra;
    elseif strcmp(output_folder,'optimalStates-Brinkman-lowRes/')
         Rm = Ra ; %Da*Ra*Pr*Pr;
    else
        Rm = Ra;
    end
    
    Rm_i = find(Rm_arr == Rm);
    CR_i = find(CR_arr == CR);
    
    % Let's skip some things for now
    if Rm > RM_MAX || CR > CR_MAX
        continue
    end
    
    optimalFluxMat(CR_i, Rm_i) = state.flux;
    optimalWidthMat(CR_i, Rm_i) = state.fullWidth;
    
    if isfield(state, 'criticalWidth')
        optimalCritWidthMat(CR_i, Rm_i) = state.criticalWidth;
    else
        optimalCritWidthMat(CR_i, Rm_i) = NaN;
    end
    
    if isfield(state, 'flagCritWidth')
        optimalFlagCritWidth(CR_i, Rm_i) = state.flagCritWidth;
    else
        optimalFlagCritWidth(CR_i, Rm_i) = false;
    end
    
    optimallMat(CR_i, Rm_i) = state.smallLScale;
    optimalHMat(CR_i, Rm_i) = state.H;
    optimalhMat(CR_i, Rm_i) = state.h;
    optimalChannelHMat(CR_i, Rm_i) = state.channelHeight;
    optimalChannelWidthMat(CR_i, Rm_i) = state.channelWidth;
    
    optimalhPorosityMat(CR_i, Rm_i) = state.hporosity;
    
    % Interior:
    %     optimalHeatAdvection(CR_i, Ra_i) = state.heatAdvectionInterior;
    %     optimalHeatDiffusion(CR_i, Ra_i) = state.heatDiffusionInterior;
    %     optimalLatentHeat(CR_i, Ra_i) = state.latentHeatInterior;
    %     optimalTempFrame(CR_i, Ra_i) = state.TFrameAdvectionInterior;
    %
    %     optimalSaltAdvection(CR_i, Ra_i) = state.saltAdvectionInterior;
    %     optimalSaltDiffusion(CR_i, Ra_i) = state.saltDiffusionInterior;
    %     optimalSolidSaltFrame(CR_i, Ra_i) = state.solidSalinityFrameInterior;
    %     optimalLiquidSaltFrame(CR_i, Ra_i) = state.liquidSalinityFrameInterior;
    %
    %     optimalVorticityDissipation(CR_i, Ra_i) = state.vorticityDiffusionInterior;
    %     optimalBaroclinicTorque(CR_i, Ra_i) = state.baroclinicTorqueInterior;
    %     optimalVorticityPermeability(CR_i, Ra_i) = state.vorticityPermeabilityInterior;
    
    optimalPorosMat(CR_i, Rm_i) = state.maxPorosityInterior;
    optimalPermMat(CR_i, Rm_i) = state.maxPermeabilityInterior;
    %optimalPsiMat(CR_i, Rm_i) = state.maxStreamfunctionInterior;
    optimalPsiMat(CR_i, Rm_i) = state.maxStreamfunctionMush;
    averagePorosMat(CR_i, Rm_i) = state.avPorosityInterior;
    averagePermMat(CR_i, Rm_i) = state.avPermeabilityMush;
   
    
    if isfield(state, 'L2FsVertDiffusion')
        optimalL2DiffusionFluxMat(CR_i, Rm_i) = state.L2FsVertDiffusion;
        optimalL2FrameFluxMat(CR_i, Rm_i) = state.L2FsVertFrame;
        optimalL2FluidFluxMat(CR_i, Rm_i) = state.L2FsVertFluid;
        
        optimalL0DiffusionFluxMat(CR_i, Rm_i) = state.L0FsVertDiffusion;
        optimalL0FrameFluxMat(CR_i, Rm_i) = state.L0FsVertFrame;
        optimalL0FluidFluxMat(CR_i, Rm_i) = state.L0FsVertFluid;
    else
        optimalL2DiffusionFluxMat(CR_i, Rm_i) = NaN;
        optimalL2FrameFluxMat(CR_i, Rm_i) = NaN;
        optimalL2FluidFluxMat(CR_i, Rm_i) = NaN;
        
        optimalL0DiffusionFluxMat(CR_i, Rm_i) = NaN;
        optimalL0FrameFluxMat(CR_i, Rm_i) = NaN;
        optimalL0FluidFluxMat(CR_i, Rm_i) = NaN;
    end
    
    % Mush:
    
    optimalSaltAdvection(CR_i, Rm_i) = state.saltAdvectionMush;
    optimalSaltDiffusion(CR_i, Rm_i) = state.saltDiffusionMush;
    optimalSolidSaltFrame(CR_i, Rm_i) = state.solidSalinityFrameMush;
    optimalLiquidSaltFrame(CR_i, Rm_i) = state.liquidSalinityFrameMush;
    
    optimalVorticityDissipation(CR_i, Rm_i) = state.vorticityDiffusionMush;
    optimalBaroclinicTorque(CR_i, Rm_i) = state.baroclinicTorqueMush;
    optimalVorticityPermeability(CR_i, Rm_i) = state.vorticityPermeabilityMush;
    
    optimalHeatAdvection(CR_i, Rm_i) = state.heatAdvectionMush;
    optimalHeatDiffusion(CR_i, Rm_i) = state.heatDiffusionMush;
    optimalLatentHeat(CR_i, Rm_i) = state.latentHeatMush;
    optimalTempFrame(CR_i, Rm_i) = state.TFrameAdvectionMush;
    
    if isfield(state, 'extrapolatedMax') && ~isnan(state.extrapolatedMax)
        optimalExtrapolationFlag(CR_i, Rm_i) = state.extrapolatedMax;
    else
        optimalExtrapolationFlag(CR_i, Rm_i) = false;
    end
    
    
    
end

% Let's calculate critical rayleigh numbers properly, by extrapolating
% fluxes to zero

Rm_crit = CR_arr;
LOpt_crit = CR_arr;
Pi_crit = CR_arr;
h_crit = CR_arr;

for CR_i = 1:length(CR_arr)
    CR = CR_arr(CR_i);
    
    fluxForCR = optimalFluxMat(CR_i, :);
    LForRm = optimalWidthMat(CR_i, :);
    actualFluxes= ~isnan(fluxForCR);
    %PiForRm = (optimalhPorosityMat(CR_i, :)).^3; % TODO check this is right
    PiForRm = (averagePermMat(CR_i, :)); % TODO check this is right
    hForRm = optimalhMat(CR_i, :); % TODO check this is right
    
    %hForRm = optimalChannelHMat(CR_i, :);
    
    % Take smallest two fluxes and fit a straight line
    fluxesToFit = [];
    RmToFit = []; LToFit = []; PiToFit = []; hToFit = [];
    for i=1:length(fluxForCR)
        if ~isnan(fluxForCR(i)) && fluxForCR(i)/CR > 1e-5
            fluxesToFit(end+1) = fluxForCR(i);
            RmToFit(end+1) = Rm_arr(i);
            LToFit(end+1) = LForRm(i);
            PiToFit(end+1) = PiForRm(i);
            hToFit(end+1) = hForRm(i);
        end
    end
    
    if length(RmToFit) > 1
        Rm_crit(CR_i) = RmToFit(1) - fluxesToFit(1)*(RmToFit(2)-RmToFit(1))/(fluxesToFit(2)-fluxesToFit(1));
        LOpt_crit(CR_i) = LToFit(1) - fluxesToFit(1)*(LToFit(2)-LToFit(1))/(fluxesToFit(2)-fluxesToFit(1));
        Pi_crit(CR_i) = PiToFit(1) - fluxesToFit(1)*(PiToFit(2)-PiToFit(1))/(fluxesToFit(2)-fluxesToFit(1));
        h_crit(CR_i) = hToFit(1) - fluxesToFit(1)*(hToFit(2)-hToFit(1))/(fluxesToFit(2)-fluxesToFit(1));
    elseif ~isempty(RmToFit)
        Rm_crit(CR_i) = RmToFit(1);
        LOpt_crit(CR_i) = LToFit(1);
        Pi_crit(CR_i) = PiToFit(1);
        h_crit(CR_i) = hToFit(1);
    else
        
        defaultVal = 0.0; %NaN
        Rm_crit(CR_i) = defaultVal;
        LOpt_crit(CR_i) = defaultVal;
        Pi_crit(CR_i) = defaultVal;
        h_crit(CR_i) = defaultVal;
        
        
    end
    
    fprintf('CR = %1.3f, Rm_crit = %1.3f, L_crit = %1.3f\n', CR, Rm_crit(CR_i), LOpt_crit(CR_i));
    
    
end

%figure(); semilogx(CR_arr, 1./L_crit); xlabel('CR'); ylabel('$a_c$');

%T = table(CR_arr, Rm_crit)




if makeStabilityPlots
    
    hStability = figure();
    set(hStability, 'Position', [100 100 1400 900]);
    
    m = 2; % Num rows
    n = 3; % Num cols
    subplot(m, n, 1);
    plot(CR_arr, 1./LOpt_crit);
    xlabel('$\mathcal{C}$');
    ylabel('$\alpha_{Opt, critical}$');
    
    
    
    subplot(m, n, 2);
    plot((CR_arr), (h_crit));
    xlabel('$\mathcal{C}$');
    ylabel('$h_{critical}$');
    
    subplot(m, n, 3);
    plot((CR_arr), (Pi_crit));
    xlabel('$\mathcal{C}$');
    ylabel('$\Pi_{critical}$');
    
    subplot(m, n, 4);
    hold on;
    plot((1./CR_arr), (Rm_crit));
    plot((1./CR_arr), (Rm_crit.^2));
    
    hold off;
    legend('$Rm_{crit}$', '$Rm_{crit}^2$');
    xlabel('$1/\mathcal{C}$');
    %ylabel('$Rm_{critical}$');
    
    subplot(m, n, 5);
    plot((1./(St*(h_crit).^2)), (Rm_crit));
    xlabel('$1/(St h^2)$');
    ylabel('$Rm_{critical}$');
    
    
    x = (CR_arr.*LOpt_crit.^2./(h_crit.*Pi_crit));
    x = (CR_arr.*LOpt_crit.^2/(h_crit.*Pi_crit));
    
    x = CR_arr; %(CR_arr./Rm_crit).^(1/2);
    
    y = Pi_crit./h_crit;
    
    y = Rm_crit.*h_crit.^2.*Pi_crit;
    
    x = 1./CR_arr.^(0.5);
    y = Rm_crit;
    
    x = log10(CR_arr);
    %y = St*Rm_crit.*h_crit.^2;
    y = St*Rm_crit.*h_crit.^2;
    
    %y = Rm_crit.*h_crit.^2;
    
    hErr = 0.4/256;
    errBar = Rm_crit.*(hErr);
    
    subplot(m, n, 6);
    maxC_i = 17;
    maxC_i = min(maxC_i, length(y));
    
    %plot(x(1:maxC_i), y(1:maxC_i));
    errorbar(x(1:maxC_i), y(1:maxC_i), errBar(1:maxC_i));
    
    
    %xlabel('$\mathcal{C} L /(\Pi h)$');
    %ylabel('$Rm_{critical}$');
    
    %xlabel('$(\mathcal{C} / Rm_{critical})^{1/2}$');
    %ylabel('$ \Pi$');
    
    %xlabel('$1/\mathcal{C}^{1/2}$');
    %ylabel('$Rm_c$');
    
    xlabel('log$_{10}(\mathcal{C})$');
    ylabel('$\St \, \hat{h}^2 \, \hat{Rm_S}$');
    
end



hScaledScalings = figure();
set(hScaledScalings, 'Position', [100 100 1600 1000]);
m = 2;
n = 3;

subplot(m, n, 1);

hold on;
CR_leg = {};
for CR_i = 1:length(CR_arr(1:8))
    
    fluxForCR = optimalFluxMat(CR_i, :);
    actualFluxes= ~isnan(fluxForCR);
    
    psiForCR = optimalPsiMat(CR_i, :);
    
    hForCR = optimalhMat(CR_i, :);
    LForCR = optimalWidthMat(CR_i, :);
    
    l= optimallMat(CR_i, :);
    
    
    %RmStar = (Rm_arr-Rm_crit(CR_i))./Rm_crit(CR_i);
    RmStar = (Rm_arr)./Rm_crit(CR_i);
    CR = CR_arr(CR_i);
    
    %RmStar = Rm_arr .* hForCR.^2;
    
    
    x = RmStar.*hForCR./LForCR;
    x = Rm_arr.*hForCR.^2 ./ LForCR;
    x = RmStar.^0.5;
    %x = CR_arr(CR_i)*LForCR;
    
    x = Rm_arr.*hForCR.^2 ./ LForCR;
    
    x = 5*RmStar.*hForCR.^3 ./ (l*CR);
    
    if sum(actualFluxes) > 0
        
        plot(x(actualFluxes), psiForCR(actualFluxes));
        CR_leg{end+1} = ['$\mathcal{C}=', num2str(CR_arr(CR_i)),'$'];
    end
    
end
hold off;
box on;
xlabel('$ (Rm^*)^{1/2} $'); 
xlabel('$ (Rm) h^2 / L $'); 
xlabel('$h^3 Rm^* St/(l \mathcal{C})$');
ylabel('$\psi$');


subplot(m, n, 2);

hold on;
leg = {};
for Ra_i = 1:length(Rm_arr)
    
    F = optimalFluxMat(:, Ra_i);
    actualFluxes= ~isnan(F).*((CR_arr < 0.25).');
    
    psi = optimalPsiMat(:, Ra_i);
    
    %h = optimalhMat(:, Ra_i);
    h = optimalChannelHMat(:, Ra_i);
    
    L = optimalWidthMat(:, Ra_i);
    
    l= optimallMat(:, Ra_i);
    
   % Rm = Rm_arr - Rm_crit(CR_i);
    
    %x = Rm.*hForCR./LForCR;
   % x = Rm_arr.*hForCR.^2 ./ LForCR;
   
   Rm = Rm_arr(Ra_i);
   RmStar = (Rm./Rm_crit).';
   
   x = CR_arr.^(0.5);
    
   x = 1./(h.^(0.5));
   x = (CR_arr.').^(0.5)./(L.*RmStar.^(1/2));
   
    %x = CR_arr(CR_i)*LForCR;
    
    if sum(actualFluxes) > 0
        
        plot(x(actualFluxes==1), ...
            psi(actualFluxes==1));
        leg{end+1} = ['$Rm=', num2str(Rm_arr(Ra_i)),'$'];
    end
    
end
hold off;
box on;
xlabel('$C^{1/2}$'); 
xlabel('$C^{1/2}/(L Rm^{1/2})$');
ylabel('$\psi$');

legend(leg, 'Location', 'eastoutside')









subplot(m, n, 3);

hold on;
CR_leg = {};
for CR_i = 1:length(CR_arr(1:7))
      %RmStar = (Rm_arr-Rm_crit(CR_i))./Rm_crit(CR_i);
    RmStar = (Rm_arr)./Rm_crit(CR_i);
    CR = CR_arr(CR_i);
    
    %RmStar = Rm_arr .* hForCR.^2;
    
    fluxForCR = optimalFluxMat(CR_i, :);
    actualFluxes= ~isnan(fluxForCR).*~isnan(RmStar);
    
    psi = optimalPsiMat(CR_i, :);
    
    h= optimalhMat(CR_i, :);
    L = optimalWidthMat(CR_i, :);
    
    H = optimalHMat(CR_i, :);
    
  l= optimallMat(CR_i, :);
    
    
    x = RmStar.*hForCR./LForCR;
    x = Rm_arr.*hForCR.^2 ./ LForCR;
    x = RmStar;
    
    x = 1./hForCR;
    
    x = hForCR;
    
    %x = CR./(5*H.^2);
    
    %x = 1./H.^2;
    
    %x = CR_arr(CR_i)*LForCR;
    
    y = psi;
    x = (CR)*hForCR./LForCR;
    
    if sum(actualFluxes) > 0
        
        plot(x(actualFluxes==1), y(actualFluxes==1));
        CR_leg{end+1} = ['$\mathcal{C}=', num2str(CR_arr(CR_i)),'$'];
    end
    
end
hold off;
box on;
%xlabel('$ 1/H^2 $'); 
xlabel('$ Rm^*$');
xlabel('$C h/L $');
ylabel('$\psi$');



subplot(m, n, 4);

hold on;
CR_leg = {};
for CR_i = 1:length(CR_arr(1:9))
     %RmStar = (Rm_arr-Rm_crit(CR_i))./Rm_crit(CR_i);
    RmStar = (Rm_arr)./Rm_crit(CR_i);
    CR = CR_arr(CR_i);
    
    %RmStar = Rm_arr .* hForCR.^2;
    
    fluxForCR = optimalFluxMat(CR_i, :);
    actualFluxes= ~isnan(fluxForCR).*~isnan(RmStar);
    
    psi = optimalPsiMat(CR_i, :);
    
    h= optimalhMat(CR_i, :);
    L = optimalWidthMat(CR_i, :);
    
    
    
    chi_h = optimalhPorosityMat(CR_i, :);
    deltaChi = 1-chi_h;
    
    
    x = RmStar.*hForCR./LForCR;
    x = Rm_arr.*hForCR.^2 ./ LForCR;
    x = (CR./RmStar).^(1/6);
    x = RmStar;
    
    %x = (RmStar/CR).^(1/6);
    
    %x = CR_arr(CR_i)*LForCR;
    
    y = h./LForCR;
    x = sqrt(RmStar);
   % y = h;
   
   y = fluxForCR;
   x = St*psi./(LForCR*CR^(1/4));
   x = St*psi./(LForCR);
    
    if sum(actualFluxes) > 0
        
        plot(x(actualFluxes==1), y(actualFluxes==1));
        CR_leg{end+1} = ['$\mathcal{C}=', num2str(CR_arr(CR_i)),'$'];
    end
    
end
hold off;
box on;
%xlabel('$ (\mathcal{C}/Rm^*)^{1/6} $'); 
%xlabel('$St \psi  / (L \mathcal{C}^{1/4})$');
xlabel('$St \psi / (L )$');
ylabel('$F_O$');




subplot(m, n, 5);

hold on;
CR_leg = {};
for CR_i = 1:length(CR_arr(1:7))
     %RmStar = (Rm_arr-Rm_crit(CR_i))./Rm_crit(CR_i);
    RmStar = (Rm_arr)./Rm_crit(CR_i);
    CR = CR_arr(CR_i);
    
    %RmStar = Rm_arr .* hForCR.^2;
    
    fluxForCR = optimalFluxMat(CR_i, :);
    actualFluxes= ~isnan(fluxForCR).*~isnan(RmStar);
    
    psi = optimalPsiMat(CR_i, :);
    
    h= optimalhMat(CR_i, :);
    L = optimalWidthMat(CR_i, :);
    l = optimallMat(CR_i, :);
    
    
    chi_h = optimalhPorosityMat(CR_i, :);
    deltaChi = 1-chi_h;
    
    
    x = RmStar.*hForCR./LForCR;
    x = Rm_arr.*hForCR.^2 ./ LForCR;
    x = (CR./RmStar).^(1/6);
    x = RmStar;
    
    %x = (RmStar/CR).^(1/6);
    
    %x = CR_arr(CR_i)*LForCR;
    
    y = h./LForCR;
    x = sqrt(RmStar);
   % y = h;
   
   x = hForCR./l;
   y = psi*sqrt(Le)/CR;
   
    
    if sum(actualFluxes) > 0
        
        plot(x(actualFluxes==1), y(actualFluxes==1));
        CR_leg{end+1} = ['$\mathcal{C}=', num2str(CR_arr(CR_i)),'$'];
    end
    
end
hold off;
box on;
%xlabel('$ (\mathcal{C}/Rm^*)^{1/6} $'); 
xlabel('$h/l$');
ylabel('$\psi \sqrt{Le} / CR$');


subplot(m, n, 6);

hold on;
CR_leg = {};
for CR_i = 1:length(CR_arr(1:7))
    RmStar = (Rm_arr-Rm_crit(CR_i))./Rm_crit(CR_i);
    %RmStar = (Rm_arr)./Rm_crit(CR_i);
    CR = CR_arr(CR_i);
    
    %RmStar = Rm_arr .* hForCR.^2;
    
    fluxForCR = optimalFluxMat(CR_i, :);
    actualFluxes= ~isnan(fluxForCR).*~isnan(RmStar);
    
    psi = optimalPsiMat(CR_i, :);
    
    h= optimalhMat(CR_i, :);
    L = optimalWidthMat(CR_i, :);
    l = optimallMat(CR_i, :);
    h_channel = optimalChannelHMat(CR_i,:);
    a = optimalChannelWidthMat(CR_i,:);
    
    chi_h = optimalhPorosityMat(CR_i, :);
    deltaChi = 1-chi_h;
    
    
    x = RmStar.*hForCR./LForCR;
    x = Rm_arr.*hForCR.^2 ./ LForCR;
    x = (CR./RmStar).^(1/6);
    x = RmStar;
    
    %x = (RmStar/CR).^(1/6);
    
    %x = CR_arr(CR_i)*LForCR;
    
    y = h./LForCR;
    x = sqrt(RmStar);
   % y = h;
   
   x = CR*hForCR.*l.*Rm_arr*St;
   %x = l.^2.*Rm_arr*St./LForCR;
   x = sqrt(RmStar/CR);
   y = hForCR./LForCR;
   
   x = h_channel./LForCR;
   y = psi;
   
    
    if sum(actualFluxes) > 0
        
        plot(x(actualFluxes==1), y(actualFluxes==1));
        CR_leg{end+1} = ['$\mathcal{C}=', num2str(CR_arr(CR_i)),'$'];
    end
    
end
hold off;
box on;
%xlabel('$ (\mathcal{C}/Rm^*)^{1/6} $'); 
%xlabel('$\mathcal{C} l h Rm St$');
%xlabel('$l^2 Rm St / L$');
xlabel('$h_{channel}/L$');
ylabel('$\psi $');

% Do some transformations


% Now do plotting

doOptimalPlot = false;

if doOptimalPlot
    
    hOptimal = figure();
    set(hOptimal, 'Position', [100 100 1400 900]);
    m = 2; % Num rows
    n = 2; % Num cols
    
    % In the first panel, plot F vs Ra for each CR
    subplot(m, n, 1);
    hold on;
    CR_leg = {};
    for CR_i = 1:length(CR_arr)
        fluxForCR = optimalFluxMat(CR_i, :);
        actualFluxes= ~isnan(fluxForCR);
        if sum(actualFluxes) > 0
            
            plot(Rm_arr(actualFluxes), fluxForCR(actualFluxes));
            CR_leg{end+1} = ['$\mathcal{C}=', num2str(CR_arr(CR_i)),'$'];
        end
    end
    hold off;
    box on;
    xlabel('$Rm$'); ylabel('$F_0$');
    %legend(CR_leg);
    
    % In the second panel, plot width vs Ra for each CR
    subplot(m, n, 2);
    hold on;
    
    for CR_i = 1:length(CR_arr)
        widthForRa = optimalWidthMat(CR_i, :);
        actualWidths= ~isnan(widthForRa);
        if sum(actualWidths) > 0
            plot(Rm_arr(actualWidths), widthForRa(actualWidths));
        end
        
    end
    hold off;
    xlabel('$Rm$'); ylabel('$L$');
    box on;
    legend(CR_leg, 'Location', 'eastoutside');
    
    % In the third panel, plot F vs CR for each Ra
    subplot(m, n, 3);
    hold on;
    leg = {};
    for CR_i = 1:length(Rm_arr)
        fluxForCR = optimalFluxMat(:, CR_i);
        actualFluxes= ~isnan(fluxForCR);
        if sum(actualFluxes) > 0
            
            plot(CR_arr(actualFluxes), fluxForCR(actualFluxes));
            leg{end+1} = ['$Rm=', num2str(Rm_arr(CR_i)),'$'];
        end
    end
    hold off;
    box on;
    xlabel('$C$'); ylabel('$F_0$');
    %legend(Ra_leg);
    
    % In the fourth panel, plot width vs CR for each Ra
    subplot(m, n, 4);
    hold on;
    leg = {};
    for CR_i = 1:length(Rm_arr)
        widthForRa = optimalWidthMat(:, CR_i) ;
        actualWidths= ~isnan(widthForRa);
        if sum(actualWidths) > 0
            
            plot(CR_arr(actualWidths), widthForRa(actualWidths));
            leg{end+1} = ['$Rm=', num2str(Rm_arr(CR_i)),'$'];
        end
    end
    hold off;
    box on;
    xlabel('$C$'); ylabel('$L$');
    legend(leg, 'Location', 'eastoutside');
    
    
    % set(hOptimal,'Units','Inches');
    % pos = get(hOptimal,'Position');
    % set(hOptimal,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    % print(gcf, '-dpdf', [figure_output_dir, optimalFluxes.pdf'], printQuality);
    %
    
    
end




% Make plots of terms in heat, salt and momentum equations

if makeDetailedPlots
    
    if scaleDetailedPlotsWithRmCrit
        Ralabel = '$$(Ra-Ra_{crit})/Ra_{crit}$$';
        Ra_reference = 0*Rm_crit;
    else
        Ralabel = '$$Rm$$';
        Ra_reference = Rm_crit;
    end
    
    
    
    
    
    
    
    hEquations = figure();
    set(hEquations, 'Position', [100 100 1600 1000]);
    m = 2;
    n = 2;
    
    subplot(m, n, 1);
    hold on;
    for CR_i = 1:length(Rm_arr)
        
        fluxForCR = optimalFluxMat(:, CR_i);
        
        advection = optimalHeatAdvection(:, CR_i);
        baroclinic = optimalHeatDiffusion(:, CR_i);
        solidFrame = optimalLatentHeat(:, CR_i);
        liquidSaltFrame = optimalTempFrame(:, CR_i);
        
        actualFluxes= ~isnan(fluxForCR);
        if sum(actualFluxes) > 0
            plot1 = plot(CR_arr(actualFluxes), advection(actualFluxes), 'x-b');
            plot2 = plot(CR_arr(actualFluxes), baroclinic(actualFluxes), 'x-r');
            plot3 = plot(CR_arr(actualFluxes), solidFrame(actualFluxes), 'x-m');
            TframeLeg  = plot(CR_arr(actualFluxes), liquidSaltFrame(actualFluxes), 'x-k');
        end
        
    end
    hold off;
    box on;
    xlabel('$C$'); ylabel('Terms');
    title('Heat  equation');
    legend([plot1, plot2, plot3, TframeLeg], ...
        {'Advection', 'Diffusion', 'Latent heat', 'Frame advection'});
    
    
    subplot(m, n, 2);
    hold on;
    for CR_i = 1:length(Rm_arr)
        
        fluxForCR = optimalFluxMat(:, CR_i);
        
        advection = optimalSaltAdvection(:, CR_i);
        baroclinic = optimalSaltDiffusion(:, CR_i);
        solidFrame = optimalSolidSaltFrame(:, CR_i);
        liquidSaltFrame = optimalLiquidSaltFrame(:, CR_i);
        
        actualFluxes= ~isnan(fluxForCR);
        if sum(actualFluxes) > 0
            plot1 = plot(CR_arr(actualFluxes), advection(actualFluxes), 'x-b');
            plot2 = plot(CR_arr(actualFluxes), baroclinic(actualFluxes), 'x-r');
            plot3 = plot(CR_arr(actualFluxes), solidFrame(actualFluxes), 'x-m');
            TframeLeg  = plot(CR_arr(actualFluxes), liquidSaltFrame(actualFluxes), 'x-k');
        end
        
    end
    hold off;
    box on;
    xlabel('$C$'); ylabel('Terms');
    title('Salt  equation');
    legend([plot1, plot2, plot3, TframeLeg], ...
        {'Advection', 'Diffusion', 'Solid frame adv', 'Liquid frame adv'});
    
    %ylim([0 50]);
    
    
    
    
    subplot(m, n, 3);
    hold on;
    for CR_i = 1:length(Rm_arr)
        
        fluxForCR = optimalFluxMat(:, CR_i);
        
        advection = optimalVorticityDissipation(:, CR_i);
        baroclinic = optimalBaroclinicTorque(:, CR_i);
        solidFrame = optimalVorticityPermeability(:, CR_i);
        
        actualFluxes= ~isnan(fluxForCR);
        if sum(actualFluxes) > 0
            plot1 = plot(CR_arr(actualFluxes), advection(actualFluxes), 'x-b');
            plot2 = plot(CR_arr(actualFluxes), baroclinic(actualFluxes), 'x-r');
            plot3 = plot(CR_arr(actualFluxes), solidFrame(actualFluxes), 'x-m');
        end
        
    end
    hold off;
    box on;
    xlabel('$C$'); ylabel('Terms');
    title('Momentum  equation');
    legend([plot1, plot2, plot3], ...
        {'Vorticity diffusion', 'Baroclinic', 'VorticityPermeability'});
    
    
    
    
    set(hEquations,'Units','Inches');
    pos = get(hEquations,'Position');
    set(hEquations,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(gcf, '-dpdf', [figure_output_dir, 'equationScaling.pdf'], printQuality);
    
    
    
    hSaltEquation = figure();
    set(hSaltEquation, 'Position', [100 100 1600 1000]);
    m = 2;
    n = 2;
    
    subplot(m, n, 1);
    hold on;
    for CR_i = 1:length(Rm_arr)
        
        fluxForCR = optimalFluxMat(:, CR_i);
        
        L = optimalWidthMat(:, CR_i);
        psi = optimalPsiMat(:, CR_i);
        
        advection = optimalSaltAdvection(:, CR_i);
        %baroclinic = optimalSaltDiffusion(:, CR_i);
        %vortPerm = optimalSolidSaltFrame(:, CR_i);
        %liquidSaltFrame = optimalLiquidSaltFrame(:, CR_i);
        
        actualFluxes= ~isnan(fluxForCR);
        if sum(actualFluxes) > 0
            plot1 = plot(psi(actualFluxes)./(L(actualFluxes).^2), advection(actualFluxes));
            %plot2 = plot(CR_arr(actualFluxes), baroclinic(actualFluxes), 'x-r');
            %plot3 = plot(CR_arr(actualFluxes), vortPerm(actualFluxes), 'x-m');
            %TframeLeg  = plot(CR_arr(actualFluxes), liquidSaltFrame(actualFluxes), 'x-k');
        end
        
    end
    hold off;
    box on;
    xlabel('$\bar{\psi}/L^2$'); ylabel('$\mathbf{U} \cdot \nabla S_l$');
    title('Salt  equation - advection');
    
    
    subplot(m, n, 2);
    hold on;
    Ra_leg = {};
    for Rm_i = 1:length(Rm_arr)
        
        fluxForCR = optimalFluxMat(:, Rm_i);
        
        L = optimalWidthMat(:, Rm_i).';
        psi = optimalPsiMat(:, Rm_i).';
        chi = optimalPorosMat(:, Rm_i).';
        
        phi = 1 - chi;
        
        
        solidFrame = optimalSolidSaltFrame(:, Rm_i);
        %liquidSaltFrame = optimalLiquidSaltFrame(:, CR_i);
        
        actualFluxes= ~isnan(fluxForCR);
        if sum(actualFluxes) > 0
            plot1 = plot(CR_arr(actualFluxes).*phi(actualFluxes)./L(actualFluxes), solidFrame(actualFluxes));
            Ra_leg{end+1} = ['$\mathcal{R} = ', num2str(Rm_arr(Rm_i)), '$'];
        end
        
    end
    hold off;
    box on;
    xlabel('$\mathcal{C} \bar{\phi}/L$'); ylabel('$\mathcal{C} d(\chi)/dz$');
    title('Salt  equation - frame advection');
    legend(Ra_leg, 'Location', 'eastoutside');
    
    
    subplot(m, n, 3);
    hold on;
    for Rm_i = 1:length(Rm_arr)
        
        fluxForCR = optimalFluxMat(:, Rm_i);
        
        L = optimalWidthMat(:, Rm_i).';
        psi = optimalPsiMat(:, Rm_i).';
        chi = optimalPorosMat(:, Rm_i).';
        
        phi = 1 - chi;
        
        solidFrame = optimalSolidSaltFrame(:, Rm_i);
        advection = optimalSaltAdvection(:, Rm_i);
        
        actualFluxes= ~isnan(fluxForCR);
        if sum(actualFluxes) > 0
            plot1 = plot(solidFrame(actualFluxes), advection(actualFluxes));
        end
        
    end
    hold off;
    box on;
    xlabel('$\mathcal{C} d(\chi)/dz$'); ylabel('$\mathbf{U} \cdot \nabla S_l$');
    title('Salt  equation - scaling argument');
    
    
    
    
    
    subplot(m, n, 4);
    hold on;
    for Rm_i = 1:length(Rm_arr)
        
        fluxForCR = optimalFluxMat(:, Rm_i);
        
        L = optimalWidthMat(:, Rm_i).';
        psi = optimalPsiMat(:, Rm_i).';
        chi = optimalPorosMat(:, Rm_i).';
        
        phi = 1 - chi;
        
        solidFrame = optimalSolidSaltFrame(:, Rm_i);
        actualFluxes= ~isnan(fluxForCR);
        if sum(actualFluxes) > 0
            plot1 = plot(CR_arr(actualFluxes).*phi(actualFluxes)./L(actualFluxes), ...
                psi(actualFluxes)./L(actualFluxes).^2);
            
        end
        
    end
    hold off;
    box on;
    xlabel('$\mathcal{C} \bar{\phi} / L$'); ylabel('$\bar{\psi} / L^2$');
    title('Salt  equation - scaling result');
    
    
    set(hSaltEquation,'Units','Inches');
    pos = get(hSaltEquation,'Position');
    set(hSaltEquation,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(gcf, '-dpdf', [figure_output_dir, 'saltEquation.pdf'], printQuality);
    
    
    
    
    hHeatEquation = figure();
    set(hHeatEquation, 'Position', [100 100 1600 1000]);
    m = 2;
    n = 2;
    
    subplot(m, n, 1);
    hold on;
    for Rm_i = 1:length(Rm_arr)
        
        fluxForCR = optimalFluxMat(:, Rm_i);
        
        L = optimalWidthMat(:, Rm_i);
        psi = optimalPsiMat(:, Rm_i);
        phi = 1-optimalPorosMat(:, Rm_i);
        
        latent = optimalLatentHeat(:, Rm_i);
        
        actualFluxes= ~isnan(fluxForCR);
        if sum(actualFluxes) > 0
            plot1 = plot(St*phi(actualFluxes)./(L(actualFluxes)), latent(actualFluxes));
        end
        
    end
    hold off;
    box on;
    xlabel('$\mathcal{S} \bar{\phi}/L$'); ylabel('$\mathcal{S} d \chi / dz$');
    title('Heat  equation - latent heat');
    
    
    subplot(m, n, 2);
    hold on;
    Ra_leg = {};
    for Rm_i = 1:length(Rm_arr)
        
        fluxForCR = optimalFluxMat(:, Rm_i);
        
        L = optimalWidthMat(:, Rm_i).';
        psi = optimalPsiMat(:, Rm_i).';
        chi = optimalPorosMat(:, Rm_i).';
        
        phi = 1 - chi;
        
        
        diffusion = optimalHeatDiffusion(:, Rm_i);
        
        actualFluxes= ~isnan(fluxForCR);
        if sum(actualFluxes) > 0
            plot1 = plot(1./L(actualFluxes).^2, diffusion(actualFluxes));
            Ra_leg{end+1} = ['$\mathcal{R} = ', num2str(Rm_arr(Rm_i)), '$'];
        end
        
    end
    hold off;
    box on;
    xlabel('$1/L^2$'); ylabel('$\nabla^2 \theta$');
    title('Heat  equation - diffusion');
    legend(Ra_leg, 'Location', 'eastoutside');
    
    
    subplot(m, n, 3);
    hold on;
    for Rm_i = 1:length(Rm_arr)
        
        fluxForCR = optimalFluxMat(:, Rm_i);
        
        L = optimalWidthMat(:, Rm_i).';
        psi = optimalPsiMat(:, Rm_i).';
        chi = optimalPorosMat(:, Rm_i).';
        
        phi = 1 - chi;
        
        latent = optimalLatentHeat(:, Rm_i);
        diffusion = optimalHeatDiffusion(:, Rm_i);
        
        actualFluxes= ~isnan(fluxForCR);
        if sum(actualFluxes) > 0
            plot1 = plot(latent(actualFluxes), diffusion(actualFluxes));
        end
        
    end
    hold off;
    box on;
    xlabel('$\mathcal{S} d(\chi)/dz$'); ylabel('$\nabla^2 \theta$');
    title('Heat  equation - scaling argument');
    
    
    
    subplot(m, n, 4);
    hold on;
    for Rm_i = 1:length(Rm_arr)
        
        fluxForCR = optimalFluxMat(:, Rm_i);
        
        L = optimalWidthMat(:, Rm_i).';
        psi = optimalPsiMat(:, Rm_i).';
        chi = optimalPorosMat(:, Rm_i).';
        
        phi = 1 - chi;
        
        solidFrame = optimalSolidSaltFrame(:, Rm_i);
        actualFluxes= ~isnan(fluxForCR);
        if sum(actualFluxes) > 0
            plot1 = plot(St*phi(actualFluxes)./L(actualFluxes), 1./L(actualFluxes).^2);
            
        end
        
    end
    hold off;
    box on;
    xlabel('$\mathcal{S} \bar{\phi}/L$'); ylabel('$1/L^2$');
    title('Heat  equation - scaling result');
    
    
    
    set(hHeatEquation,'Units','Inches');
    pos = get(hHeatEquation,'Position');
    set(hHeatEquation,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(gcf, '-dpdf', [figure_output_dir, 'heatEquation.pdf'], printQuality);
    
    
    
    
    
    hMomEquation = figure();
    set(hMomEquation, 'Position', [100 100 1600 1000]);
    m = 2; % Num rows
    n = 3; % Num cols
    
    subplot(m, n, 1);
    hold on;
    for Rm_i = 1:length(Rm_arr)
        
        fluxForCR = optimalFluxMat(:, Rm_i);
        
        L = optimalWidthMat(:, Rm_i);
        psi = optimalPsiMat(:, Rm_i);
        perm = optimalPermMat(:, Rm_i);
        Ra = Rm_arr(Rm_i);
        
        
        torque = optimalBaroclinicTorque(:, Rm_i);
        
        actualFluxes= ~isnan(fluxForCR);
        if sum(actualFluxes) > 0
            plot1 = plot(Ra*perm(actualFluxes)./(L(actualFluxes)), torque(actualFluxes));
        end
        
    end
    hold off;
    box on;
    xlabel('$\mathcal{R} \bar{\Pi}/L$'); ylabel('$\mathcal{R} \Pi d\theta/dx$');
    title('Baroclinic torque');
    
    subplot(m, n, 2);
    hold on;
    for Rm_i = 1:length(Rm_arr)
        
        fluxForCR = optimalFluxMat(:, Rm_i);
        
        L = optimalWidthMat(:, Rm_i);
        psi = optimalPsiMat(:, Rm_i);
        perm = optimalPermMat(:, Rm_i);
        Ra = Rm_arr(Rm_i);
        
        
        diss = optimalVorticityDissipation(:, Rm_i);
        
        actualFluxes= ~isnan(fluxForCR);
        if sum(actualFluxes) > 0
            plot1 = plot(psi(actualFluxes)./(L(actualFluxes).^2), diss(actualFluxes));
        end
        
    end
    hold off;
    box on;
    xlabel('$\bar{\psi}/L^2$'); ylabel('$\nabla^2 \psi$');
    title('Vorticity dissipation');
    
    
    
    subplot(m, n, 3);
    hold on;
    Ra_leg = {};
    for Rm_i = 1:length(Rm_arr)
        
        fluxForCR = optimalFluxMat(:, Rm_i);
        
        L = optimalWidthMat(:, Rm_i).';
        psi = optimalPsiMat(:, Rm_i).';
        chi = optimalPorosMat(:, Rm_i).';
        Pi = optimalPermMat(:, Rm_i).';
        
        phi = 1 - chi;
        
        
        vortPerm = optimalVorticityPermeability(:, Rm_i);
        
        actualFluxes= ~isnan(fluxForCR);
        if sum(actualFluxes) > 0
            plot1 = plot(psi(actualFluxes)./(L(actualFluxes).^2.*Pi(actualFluxes)), ...
                vortPerm(actualFluxes));
            Ra_leg{end+1} = ['$\mathcal{R} = ', num2str(Rm_arr(Rm_i)), '$'];
        end
        
    end
    hold off;
    box on;
    xlabel('$\bar{\psi} / (L^2 \bar{\Pi})$'); ylabel('$(\nabla \psi \cdot \nabla \Pi) / \Pi$');
    title('Vorticity-permeability term');
    
    
    subplot(m, n, 4);
    hold on;
    for Rm_i = 1:length(Rm_arr)
        
        fluxForCR = optimalFluxMat(:, Rm_i);
        
        L = optimalWidthMat(:, Rm_i).';
        psi = optimalPsiMat(:, Rm_i).';
        chi = optimalPorosMat(:, Rm_i).';
        
        phi = 1 - chi;
        
        baroclinic = optimalBaroclinicTorque(:, Rm_i);
        vortPerm = optimalVorticityPermeability(:, Rm_i);
        
        actualFluxes= ~isnan(fluxForCR);
        if sum(actualFluxes) > 0
            plot1 = plot(baroclinic(actualFluxes), vortPerm(actualFluxes));
        end
        
    end
    hold off;
    box on;
    xlabel('$\mathcal{R} \Pi d\theta/dx$'); ylabel('$(\nabla \psi \cdot \nabla \Pi) / \Pi$');
    title('Scaling argument');
    
    
    
    subplot(m, n, [5, 6]);
    hold on;
    for Rm_i = 1:length(Rm_arr)
        
        fluxForCR = optimalFluxMat(:, Rm_i);
        
        L = optimalWidthMat(:, Rm_i).';
        psi = optimalPsiMat(:, Rm_i).';
        chi = optimalPorosMat(:, Rm_i).';
        Pi = optimalPermMat(:, Rm_i).';
        
        phi = 1 - chi;
        
        solidFrame = optimalSolidSaltFrame(:, Rm_i);
        actualFluxes= ~isnan(fluxForCR);
        if sum(actualFluxes) > 0
            plot1 = plot(psi(actualFluxes)./(Pi(actualFluxes).*L(actualFluxes).^2), ...
                Rm_arr(Rm_i)*Pi(actualFluxes)./L(actualFluxes));
            
        end
        
    end
    hold off;
    box on;
    xlabel('$\bar{\psi}/(\bar{\Pi} L^2)$'); ylabel('$\mathcal{R} \bar{\Pi}/L$');
    title('Scaling result');
    legend(Ra_leg, 'Location', 'eastoutside');
    
    
    
    set(hMomEquation,'Units','Inches');
    pos = get(hMomEquation,'Position');
    set(hMomEquation,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(gcf, '-dpdf', [figure_output_dir, 'momEquation.pdf'], printQuality);
    
    
    
    
end % end detailed plots




hFluxWidth = figure();
set(hFluxWidth, 'Position', [100 100 1600 1000]);
m = 2;
n = 3;

subplot(m, n, 1);

CR_mat = repmat(CR_arr, length(Rm_arr), 1).';
plotField = optimalFluxMat./(CR_mat.^0);
makeOptimalPlot(Rm_arr, CR_arr, plotField, optimalFluxMat, 0, ...
    optimalExtrapolationFlag, '$Rm$', '$F_O$', false);
%xlim([0 250]);

subplot(m, n, 2);
plotField = optimalFluxMat./(CR_mat.^1);
makeOptimalPlot(Rm_arr, CR_arr, plotField, optimalFluxMat, 0, ...
    optimalExtrapolationFlag, '$Rm$', '$F_O/(\mathcal{C}^{0})$', false);


subplot(m, n, 3);
makeOptimalPlot(Rm_arr, CR_arr, optimalWidthMat, optimalFluxMat, 0, ...
    optimalExtrapolationFlag, '$Rm$','$L_O$', true);

subplot(m, n, 4);
makeOptimalPlot(Rm_arr, CR_arr, optimalFluxMat, optimalFluxMat, 1, ...
    optimalExtrapolationFlag, '$\mathcal{C}$','$F_O$', false);

subplot(m, n, 5);
plotField = log10(optimalFluxMat);
makeOptimalPlot(Rm_arr, CR_arr, plotField, optimalFluxMat, 1, ...
    optimalExtrapolationFlag, '$log10(\mathcal{C})$','$log10(F_O)$', false, [], 1);
ylim([-4 0]);

subplot(m, n, 6);
makeOptimalPlot(Rm_arr, CR_arr, optimalWidthMat, optimalFluxMat, 1, ...
    optimalExtrapolationFlag, '$\mathcal{C}$','$L_O$', true, [], 1);

set(hFluxWidth,'Units','Inches');
pos = get(hFluxWidth,'Position');
set(hFluxWidth,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf, '-dpdf', [figure_output_dir, 'fluxesWidths.pdf'], printQuality);


hFluxWidth2 = figure();
set(hFluxWidth2, 'Position', [100 100 1600 1000]);
m = 2;
n = 3;

title('Scaling for low C');

subplot(m, n, 1);

% Scale fluxes with C^(5/3), and plot against Ra-Ra_crit
SmallC = CR_arr(1:8);
% 
% CR_mat = repmat(CR_arr(1:8), length(Rm_arr), 1).';
% 
% plotField = (optimalFluxMat(1:8, :)./(CR_mat.^(1))).^(1);
% 
% makeOptimalPlot(Rm_arr, SmallC, plotField, optimalFluxMat, 0, ...
%     optimalExtrapolationFlag, '$(Ra-Ra_{crit})/Ra_{crit}$', '$(F_O/\mathcal{C})$', ...
%     false, Rm_crit, 0);

fluxPower = 1;

hold on;
for CR_i = 1:length(SmallC)
        
        flux = optimalFluxMat(CR_i,:);
        CR = SmallC(CR_i);
        L = optimallMat(CR_i,:);
        h = optimalhMat(CR_i,:);
        psi = optimalPsiMat(CR_i,:);
        
        actualFluxes= ~isnan(flux);
        
        % Need to use this form
        RmStar = (Rm_arr-Rm_crit(CR_i))/Rm_crit(CR_i);
        
        %L = L(actualFluxes);
        %RmStar = Rm_arr/Rm_crit(CR_i);
        
        
        x =(St/Le) * CR^(1)*RmStar;
        
        x =(1/(St*Le^(3/2))) * CR^(1)*RmStar;
        
        %x = sqrt(CR)*sqrt(RmStar)*5;
        
        %x = 5*psi.*h./L;
        
        y = flux.^(fluxPower);
        
        if sum(actualFluxes) > 0
            plot1 = plot(x(actualFluxes) , y(actualFluxes));
            CR_leg{end+1} = ['$\mathcal{C}=', num2str(CR), '$'];
        end
    end
    hold off;
    box on;
    
    xlabel(['$ (1 / St  Le^{3/2}) \mathcal{C} (Rm^*)^{',num2str(fluxPower),'}$']); 
   % xlabel('$ \sqrt{\mathcal{C} Rm^*} St $'); 
    %xlabel('$\psi_a St h / L$');
    ylabel('$F_O$');
   
    

LRmPower = 2;
    
subplot(m, n, 2);

hold on;
for Ra_i = 1:length(Rm_arr)
        
        flux = optimalFluxMat(:, Ra_i);
        Rm = Rm_arr(:, Ra_i);
        L = optimalWidthMat(:, Ra_i);
        
        actualFluxes= ~isnan(flux);
        
        RmStar = Rm./Rm_crit;
       % RmStar = (Rm_arr-Rm_crit(CR_i))./Rm_crit(CR_i);
        
        L = L(actualFluxes);
        
        x = CR_arr.^(1/4)./RmStar;
        x = 1./CR_arr;
        
        if sum(actualFluxes) > 0
            plot1 = plot(x(actualFluxes) , L);
            
        end
    end
    hold off;
    box on;
    xlabel(['$  1/ \mathcal{C} $']); 
    ylabel('$L$');

% hold on;
% for CR_i = 1:length(SmallC)
%         
%         flux = optimalFluxMat(CR_i,:);
%         CR = SmallC(CR_i);
%         L = optimalWidthMat(CR_i,:);
%         
%         actualFluxes= ~isnan(flux);
%         
%         RmStar = Rm_arr/Rm_crit(CR_i);
%        % RmStar = (Rm_arr-Rm_crit(CR_i))./Rm_crit(CR_i);
%         
%         L = L(actualFluxes);
%         
%         x = CR^(1/4)./RmStar;
%         x = CR^(1/3)./(RmStar.^LRmPower);
%         
%         if sum(actualFluxes) > 0
%             plot1 = plot(x(actualFluxes) , L);
%             CR_leg{end+1} = ['$\mathcal{C}=', num2str(CR), '$'];
%         end
%     end
%     hold off;
%     box on;
%     xlabel(['$   \mathcal{C}^{1/3} / (Rm^*)^{',num2str(LRmPower),'} $']); 
%     ylabel('$L$');
   
%     
% subplot(m, n, 3);
% 
% hold on;
% for CR_i = 1:length(SmallC)
%         
%         flux = optimalFluxMat(CR_i,:);
%         CR = SmallC(CR_i);
%         L = optimalWidthMat(CR_i,:);
%         
%         H = optimalHMat(CR_i, :);
%         
%         actualFluxes= ~isnan(flux);
%         
%         RmStar = Rm_arr/Rm_crit(CR_i);
%         
%         L = L(actualFluxes);
%         H = H(actualFluxes);
%         
%         if sum(actualFluxes) > 0
%             plot1 = plot(CR./RmStar(actualFluxes) , H );
%             CR_leg{end+1} = ['$\mathcal{C}=', num2str(CR), '$'];
%         end
%     end
%     hold off;
%     box on;
%     xlabel('$  \mathcal{C} / (Rm^*)^{1/2} $'); ylabel('$H$');
%    
   

hRmPower = LRmPower - 0.5;
   
subplot(m, n, 3);

hold on;
for CR_i = 1:length(SmallC)
        
        flux = optimalFluxMat(CR_i,:);
        CR = SmallC(CR_i);
        L = optimalWidthMat(CR_i,:);
        
        H = optimalHMat(CR_i, :);
        h = optimalhMat(CR_i, :);
        
        actualFluxes= ~isnan(flux);
        
        RmStar = Rm_arr/Rm_crit(CR_i);
        
        L = L(actualFluxes);
        H = H(actualFluxes);
        h = h(actualFluxes);
        
        x = RmStar;
        y = St*h/(CR^(1/3));
        
        if sum(actualFluxes) > 0
            plot1 = plot(x(actualFluxes) , y );
            CR_leg{end+1} = ['$\mathcal{C}=', num2str(CR), '$'];
        end
    end
    hold off;
    box on;
    xlabel(['$ (Rm^*) $']); 
    ylabel('$St \, h/(\mathcal{C}^{1/3})$');
   
    
    
    
    
  lRmPower = -1;
   
subplot(m, n, 4);

hold on;
for CR_i = 1:length(SmallC)
        
        flux = optimalFluxMat(CR_i,:);
        CR = SmallC(CR_i);
        L = optimalWidthMat(CR_i,:);
        l = optimallMat(CR_i,:);
        
        H = optimalHMat(CR_i, :);
        h = optimalhMat(CR_i, :);
        
        actualFluxes= ~isnan(flux);
        
        RmStar = Rm_arr/Rm_crit(CR_i);
        
        L = L(actualFluxes);
        H = H(actualFluxes);
        h = h(actualFluxes);
        l = l(actualFluxes);
        
        x = CR^(1/3)*(St/Le)*(RmStar).^(lRmPower);
        y = l;
        
        if sum(actualFluxes) > 0
            plot1 = plot(x(actualFluxes) , y );
            CR_leg{end+1} = ['$\mathcal{C}=', num2str(CR), '$'];
        end
    end
    hold off;
    box on;
    xlabel(['$ \mathcal{C}^{1/3} (Rm^*)^{',num2str(lRmPower),'} St/Le $']); 
    ylabel('$l$');
   
    
    
    
    
    
    
    
    
    subplot(m, n, 5);

hold on;
for CR_i = 1:length(SmallC)
        
        flux = optimalFluxMat(CR_i,:);
        CR = SmallC(CR_i);
        L = optimalWidthMat(CR_i,:);
        psi = optimalPsiMat(CR_i, :);
        
        H = optimalHMat(CR_i, :);
        
        actualFluxes= ~isnan(flux);
        
        RmStar = Rm_arr/Rm_crit(CR_i);
        
        L = L(actualFluxes);
        H = H(actualFluxes);
        
        x = CR^(1/2).*RmStar.^(1/2)/Le^(1/2);
        
      % x = CR^(4/4).*Rm_arr.^(1/2);
        
        if sum(actualFluxes) > 0
            plot1 = plot(x(actualFluxes) , psi(actualFluxes) );
            CR_leg{end+1} = ['$\mathcal{C}=', num2str(CR), '$'];
        end
    end
    hold off;
    box on;
    xlabel('$  \sqrt{\mathcal{C} Rm^* / Le}$'); 
    ylabel('$\psi$');
   
    
    
    
    
    subplot(m, n, 6);

hold on;
for Ra_i = 1:length(Rm_arr)
        
        flux = optimalFluxMat(:, Ra_i);
        
        L = optimalWidthMat(:, Ra_i);
        
        H = optimalHMat(:, Ra_i);
        h = optimalhMat(:, Ra_i);
        
        actualFluxes= ~isnan(flux);
        
       % RmStar = Rm_arr/Rm_crit(CR_i);
        
        L = L(actualFluxes);
        H = H(actualFluxes);
        h = h(actualFluxes);
        
        x = (1./RmStar).^(hRmPower);
        x = CR_arr.^(1/4);
        y = h;
        
        if sum(actualFluxes) > 0
            plot1 = plot(x(actualFluxes) , y );
            %CR_leg{end+1} = ['$\mathcal{C}=', num2str(CR), '$'];
        end
    end
    hold off;
    box on;
    xlabel(['$ \mathcal{C}^{1/4} $']); 
    ylabel('$h$');
    xlim([0 0.5]);
   
    
%xlim([-0.2 10])



% makeOptimalPlot(Ra_arr, CR_arr, optimalWidthMat, optimalFluxMat, 0, ...
%     optimalExtrapolationFlag, '$L$', true);
%
% subplot(m, n, 3);
% makeOptimalPlot(Ra_arr, CR_arr, optimalFluxMat, optimalFluxMat, 1, ...
%     optimalExtrapolationFlag, '$F_O$', false);
%
% subplot(m, n, 4);
% makeOptimalPlot(Ra_arr, CR_arr, optimalWidthMat, optimalFluxMat, 1, ...
%     optimalExtrapolationFlag, '$F_O$', true);


makeNextPlot = false;
if makeNextPlot

hOptimalCritWidths = figure();
set(hOptimalCritWidths, 'Position', [100 100 1600 1000]);
m=2;
n=2;
subplot(m, n, 1);
%makeOptimalPlot(CR_arr, optimalWidthMat, optimalFluxMat, 1, '$L$', true);
makeOptimalPlot(Rm_arr, CR_arr, optimalWidthMat, optimalFluxMat, 1, ...
    optimalExtrapolationFlag, '$\mathcal{C}$', '$L_O$', false, [], ...
    0 ...% semi log x
    );


subplot(m, n, 2);
makeOptimalPlot(Rm_arr, CR_arr, optimalCritWidthMat, optimalFluxMat, 1, ...
    optimalFlagCritWidth+optimalExtrapolationFlag, '$\mathcal{C}$','$L_{crit}$', true, [], ...
    0 ... % semi log x
    );


subplot(m, n, 3);
makeOptimalPlot(Rm_arr, CR_arr, optimalWidthMat-optimalCritWidthMat, optimalFluxMat, 1, ...
    optimalFlagCritWidth+optimalExtrapolationFlag, '$\mathcal{C}$','$L_O - L_{crit}$', false, [], ...
    0 ... % semi log x
    );
ylim([0, 0.5]);

end

if makeStandardPlots
    
    
    hWidths = figure();
    set(hWidths, 'Position', [100 100 1600 1000]);
    m = 2;
    n = 5; % num columns
    
    subplot(m, n, 1);
    hold on;
    Ra_leg = {};
    for CR_i = 1:length(Rm_arr)
        
        fluxForCR = optimalFluxMat(:, CR_i);
        Ra = Rm_arr(CR_i);
        smallL = optimallMat(:, CR_i);
        
        actualFluxes= ~isnan(fluxForCR);
        if sum(actualFluxes) > 0
            plot1 = plot(CR_arr(actualFluxes), smallL(actualFluxes));
            Ra_leg{end+1} = ['$\mathcal{R} = ', num2str(Ra), '$'];
        end
    end
    hold off;
    box on;
    xlabel('$\mathcal{C}$'); ylabel('$l$');
    title('l');
    %legend(Ra_leg, 'Location', 'eastoutside');
    
    
    subplot(m, n, 2);
    hold on;
    Ra_leg = {};
    for CR_i = 1:length(Rm_arr)
        
        fluxForCR = optimalFluxMat(:, CR_i);
        Ra = Rm_arr(CR_i);
        H = optimalHMat(:, CR_i);
        
        actualFluxes= ~isnan(fluxForCR);
        if sum(actualFluxes) > 0
            plot1 = plot(CR_arr(actualFluxes), H(actualFluxes));
            Ra_leg{end+1} = ['$\mathcal{R} = ', num2str(Ra), '$'];
        end
    end
    hold off;
    box on;
    xlabel('$\mathcal{C}$'); ylabel('$H$');
    title('Mush height');
    %legend(Ra_leg, 'Location', 'eastoutside');
    
    
    
    subplot(m, n, 3);
    hold on;
    Ra_leg = {};
    for CR_i = 1:length(Rm_arr)
        
        fluxForCR = optimalFluxMat(:, CR_i);
        Ra = Rm_arr(CR_i);
        h = optimalhMat(:, CR_i);
        
        actualFluxes= ~isnan(fluxForCR);
        if sum(actualFluxes) > 0
            plot1 = plot(CR_arr(actualFluxes), h(actualFluxes));
            Ra_leg{end+1} = ['$\mathcal{R} = ', num2str(Ra), '$'];
        end
    end
    hold off;
    box on;
    xlabel('$\mathcal{C}$'); ylabel('$h$');
    title('Boundary layer height');
    
    subplot(m, n, [4,5]);
    hold on;
    Ra_leg = {};
    for CR_i = 1:length(Rm_arr)
        
        fluxForCR = optimalFluxMat(:, CR_i);
        Ra = Rm_arr(CR_i);
        L = optimalWidthMat(:, CR_i);
        
        actualFluxes= ~isnan(fluxForCR);
        if sum(actualFluxes) > 0
            plot1 = plot(CR_arr(actualFluxes), L(actualFluxes));
            Ra_leg{end+1} = ['$\mathcal{R} = ', num2str(Ra), '$'];
        end
    end
    hold off;
    box on;
    xlabel('$\mathcal{C}$'); ylabel('$L$');
    title('Mush width');
    legend(Ra_leg, 'Location', 'eastoutside');
    
    
    
    
    
    subplot(m, n, 6);
    hold on; CR_leg = {};
    for CR_i = 1:length(CR_arr)
        
        fluxForCR = optimalFluxMat(CR_i,:);
        CR = CR_arr(CR_i);
        L = optimallMat(CR_i,:);
        
        actualFluxes= ~isnan(fluxForCR);
        
        L = L(actualFluxes);
        
        if sum(actualFluxes) > 0
            plot1 = plot(Rm_arr(actualFluxes) , L);
            CR_leg{end+1} = ['$\mathcal{C}=', num2str(CR), '$'];
        end
    end
    hold off;
    box on;
    xlabel('$\mathcal{R}$'); ylabel('$l$');
    title('Optimal small L');
    %legend(CR_leg, 'Location', 'eastoutside');
    
    
    
    
    subplot(m, n, 7);
    hold on; CR_leg = {};
    for CR_i = 1:length(CR_arr)
        
        fluxForCR = optimalFluxMat(CR_i,:);
        CR = CR_arr(CR_i);
        L = optimalHMat(CR_i,:);
        
        actualFluxes= ~isnan(fluxForCR);
        
        L = L(actualFluxes);
        
        if sum(actualFluxes) > 0
            plot1 = plot(Rm_arr(actualFluxes) , L);
            CR_leg{end+1} = ['$\mathcal{C}=', num2str(CR), '$'];
        end
    end
    hold off;
    box on;
    xlabel('$\mathcal{R}$'); ylabel('$H$');
    title('Mush height');
    %legend(CR_leg, 'Location', 'eastoutside');
    
    
    
    
    subplot(m, n, 8);
    hold on; CR_leg = {};
    for CR_i = 1:length(CR_arr)
        
        fluxForCR = optimalFluxMat(CR_i,:);
        CR = CR_arr(CR_i);
        H = optimalhMat(CR_i,:);
        
        actualFluxes= ~isnan(fluxForCR);
        
        H= H(actualFluxes);
        
        if sum(actualFluxes) > 0
            plot1 = plot(Rm_arr(actualFluxes) , H);
            CR_leg{end+1} = ['$\mathcal{C}=', num2str(CR), '$'];
        end
    end
    hold off;
    box on;
    xlabel('$\mathcal{R}$'); ylabel('$h$');
    title('Boundary layer height');
    
    
    subplot(m, n, [9,10]);
    hold on; CR_leg = {};
    for CR_i = 1:length(CR_arr)
        
        fluxForCR = optimalFluxMat(CR_i,:);
        CR = CR_arr(CR_i);
        H = optimalWidthMat(CR_i,:);
        
        actualFluxes= ~isnan(fluxForCR);
        
        H= H(actualFluxes);
        
        if sum(actualFluxes) > 0
            plot1 = plot(Rm_arr(actualFluxes) , H);
            CR_leg{end+1} = ['$\mathcal{C}=', num2str(CR), '$'];
        end
    end
    hold off;
    box on;
    xlabel('$\mathcal{R}$'); ylabel('$L$');
    title('Mush width');
    legend(CR_leg, 'Location', 'eastoutside');
    
    
    set(hWidths,'Units','Inches');
    pos = get(hWidths,'Position');
    set(hWidths,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(gcf, '-dpdf', [figure_output_dir, 'lengthScales.pdf'], printQuality);
    
    
    
    hLengthScalings = figure();
    set(hLengthScalings, 'Position', [100 100 1600 1000]);
    m = 2;
    n = 2;
    
    
    subplot(m, n, 1);
    hold on; CR_leg = {};
    for CR_i = 1:length(CR_arr)
        
        fluxForCR = optimalFluxMat(CR_i,:);
        CR = CR_arr(CR_i);
        h = optimalhMat(CR_i,:);
        
        
        phi = 1-optimalPorosMat(CR_i, :);
        
        actualFluxes= ~isnan(fluxForCR);
        
        h= h(actualFluxes);
        
        if sum(actualFluxes) > 0
            plot1 = plot(CR^(1/2)./Rm_arr(actualFluxes).^(1/3) , h);
            CR_leg{end+1} = ['$\mathcal{C}=', num2str(CR), '$'];
        end
    end
    hold off;
    box on;
    xlabel('$\mathcal{C}^{1/2} / Ra^{1/3}$'); ylabel('$h$');
    title('Boundary layer height');
    
    subplot(m, n, 2);
    hold on; CR_leg = {};
    for CR_i = 1:length(CR_arr)
        
        fluxForCR = optimalFluxMat(CR_i,:);
        CR = CR_arr(CR_i);
        H = optimalHMat(CR_i,:);
        h = optimalhMat(CR_i,:);
        L = optimalWidthMat(CR_i, :);
        
        phi = 1-optimalPorosMat(CR_i, :);
        
        actualFluxes= ~isnan(fluxForCR);
        
        H= H(actualFluxes);
        h= h(actualFluxes);
        L = L(actualFluxes);
        
        if sum(actualFluxes) > 0
            plot1 = plot(h , H);
            CR_leg{end+1} = ['$\mathcal{C}=', num2str(CR), '$'];
        end
    end
    hold off;
    box on;
    xlabel('$h$'); ylabel('$H$');
    title('Boundary layer height');
    legend(CR_leg, 'Location', 'eastoutside');
    
    
    
    subplot(m, n, 3);
    hold on; CR_leg = {};
    for CR_i = 1:length(CR_arr)
        
        fluxForCR = optimalFluxMat(CR_i,:);
        CR = CR_arr(CR_i);
        h = optimalhMat(CR_i,:);
        
        
        phi = 1-optimalPorosMat(CR_i, :);
        perm = optimalPermMat(CR_i, :);
        
        actualFluxes= ~isnan(fluxForCR);
        
        h= h(actualFluxes);
        perm = perm(actualFluxes);
        
        if sum(actualFluxes) > 0
            plot1 = plot(perm.^0.5 , h);
            CR_leg{end+1} = ['$\mathcal{C}=', num2str(CR), '$'];
        end
    end
    hold off;
    box on;
    xlabel('$\Pi$'); ylabel('$h$');
    title('Boundary layer height');
    
    
    
    
    
    
    set(hLengthScalings,'Units','Inches');
    pos = get(hLengthScalings,'Position');
    set(hLengthScalings,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(gcf, '-dpdf', [figure_output_dir, 'lengthScaling.pdf'], printQuality);
    
    
    
    
    
    hFields = figure();
    set(hFields, 'Position', [100 100 1600 1000]);
    m = 2;
    n = 2;
    
    subplot(m, n, 1);
    hold on;
    minChi = 1; maxChi = 0;
    for Rm_i = 1:length(Rm_arr)
        
        fluxForCR = optimalFluxMat(:, Rm_i);
        Ra = Rm_arr(Rm_i);
        chi = averagePorosMat(:, Rm_i);
        chi = optimalPorosMat(:, Rm_i);
        
        actualFluxes= ~isnan(fluxForCR);
        if sum(actualFluxes) > 0
            plot1 = plot((CR_arr(actualFluxes)/Ra).^(1/6), chi(actualFluxes));
            
            minChi = min(minChi, min(chi(actualFluxes)));
            maxChi = max(maxChi, max(chi(actualFluxes)));
        end
    end
    hold off;
    box on;
    xlabel('$(\mathcal{C}/\mathcal{R})^{1/6}$'); ylabel('$\chi$');
    title('Median porosity');
    if isnan(minChi) || isinf(minChi)
        minChi = 0;
    end
    
    if isnan(maxChi) || isinf(maxChi)
        maxChi = 1;
    end
    
    if maxChi < minChi
        maxChi = minChi + 0.1;
    end
    
    ylim([minChi*0.9 maxChi*1.1]);
    
    
    subplot(m, n, 2);
    hold on;
    for CR_i = 1:length(CR_arr)
        
        fluxForCR = optimalFluxMat(CR_i, :);
        %Ra = Ra_arr(Ra_i);
        Pi = averagePermMat(CR_i, :);
        Pi = optimalPermMat(CR_i, :);
        
        CR = CR_arr(CR_i);
        
        RmCrit = Rm_crit(CR_i);
        
        Rm_star = (Rm_arr-RmCrit)/RmCrit;
        
        actualFluxes= ~isnan(fluxForCR);
        if sum(actualFluxes) > 0
            %plot1 = plot((CR./Rm_star).^(1/2), Pi(actualFluxes));
            plot1 = plot((CR./Rm_star(actualFluxes)), Pi(actualFluxes));
        end
    end
    hold off;
    box on;
    %xlabel('$(\mathcal{C}/\mathcal{R})^{1/2}$'); ylabel('$\Pi$');
    xlabel('$(\mathcal{C}/\mathcal{Ra}^*)$'); ylabel('$\Pi$');
    title('Median permeability');
    ylim([0 1]);
    
    
    subplot(m, n, 3);
    hold on;
    for Rm_i = 1:length(Rm_arr)
        
        fluxForCR = optimalFluxMat(:, Rm_i);
        Ra = Rm_arr(Rm_i);
        psi = optimalPsiMat(:, Rm_i).';
        L = optimalWidthMat(:, Rm_i);
        phi = 1 - optimalPorosMat(:, Rm_i);
        H = optimalHMat(:, Rm_i);
        
        actualFluxes= ~isnan(fluxForCR);
        
        L = L(actualFluxes).';
        phi = phi(actualFluxes).';
        H = H(actualFluxes).';
        
        if sum(actualFluxes) > 0
            % plot(CR_arr(actualFluxes).*phi(actualFluxes)./L(actualFluxes), ...
            %   psi(actualFluxes)./L(actualFluxes).^2);
            %plot1 = plot(CR_arr(actualFluxes).*L.*phi , psi(actualFluxes));
            
            %plot1 = plot(phi.*CR_arr(actualFluxes).*L, psi(actualFluxes));
            plot1 = plot(phi.*CR_arr(actualFluxes)./(L.*H), psi(actualFluxes)./L.^2);
        end
    end
    hold off;
    box on;
    xlabel('$\bar{\phi} \mathcal{C}/L$');
    %xlabel('$\mathcal{C} \mathcal{R}^{1/2}$');
    ylabel('$\bar{\psi}/L^2$');
    title('Streamfunction');
    
    subplot(m, n, 4);
    hold on;
    Ra_leg = {};
    for Rm_i = 1:length(Rm_arr)
        
        fluxForCR = optimalFluxMat(:, Rm_i);
        Ra = Rm_arr(Rm_i);
        psi = optimalPsiMat(:, Rm_i).';
        L = optimalWidthMat(:, Rm_i);
        phi = 1 - optimalPorosMat(:, Rm_i);
        
        actualFluxes= ~isnan(fluxForCR);
        
        L = L(actualFluxes).';
        phi = phi(actualFluxes).';
        
        if sum(actualFluxes) > 0
            %plot1 = plot(CR_arr(actualFluxes).*L.*phi , psi(actualFluxes));
            %plot1 = plot(CR_arr(actualFluxes)*Ra^(1/2)/St, psi(actualFluxes)); %
            
            plot1 = plot((L.*CR_arr(actualFluxes).^2*Ra/St).*(1/3), psi(actualFluxes));
            
            Ra_leg{end+1} = ['$\mathcal{R}=', num2str(Ra), '$'];
        end
    end
    hold off;
    box on;
    %xlabel('$\mathcal{C} \mathcal{R}^{1/2}/\mathcal{S}$');
    %ylabel('$\bar{\psi}$');
    xlabel('$\mathcal{C}^{2/3} (L Ra/St)^{1/3} $');
    ylabel('$\bar{\psi}$');
    
    title('Streamfunction');
    legend(Ra_leg, 'Location', 'eastoutside');
    
    
    
    set(hFields,'Units','Inches');
    pos = get(hFields,'Position');
    set(hFields,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(gcf, '-dpdf', [figure_output_dir, 'fieldsScaling.pdf'], printQuality);
    
    
end % end if make standard plots




% Special plots for e.g. presentations
if makeOceanPlots
    
    hSaltFlux = figure();
    set(hSaltFlux, 'Position', [300 300 800 500]);
    hold on;
    salinity_legend = {};
    
    for CR_i = 1:length(CR_arr)
        
        fluxForCR = optimalFluxMat(CR_i,:);
        CR = CR_arr(CR_i);
        psi = optimalPsiMat(CR_i, :);
        L = optimalWidthMat(CR_i,:);
        Pi = optimalPermMat(CR_i,:);
        
        actualFluxes= ~isnan(fluxForCR);
        
        L = L(actualFluxes);
        Pi = Pi(actualFluxes);
        
        Se = 230;
        S0 = CR*Se/(CR+1);
        
        label = sprintf('$S_0 = %2.0f$', S0);
        
        if sum(actualFluxes) > 0
            salinity_legend{end+1} = label;
            plot1 = plot(Rm_arr(actualFluxes) , fluxForCR(actualFluxes)/CR);
            
        end
    end
    hold off;
    box on;
    xlabel('Ra'); ylabel('$F/S_0$');
    title('Ice-ocean salt flux');
    
    legend(salinity_legend, 'location', 'eastoutside');
    
    set(hSaltFlux,'Units','Inches');
    pos = get(hSaltFlux,'Position');
    set(hSaltFlux,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(gcf, '-dpdf', [figure_output_dir, 'saltFlux.pdf'], printQuality);
    
    
    
    hFluxCR = figure();
    set(hFluxCR, 'Position', [300 300 800 500]);
    hold on;
    Ra_leg = {};
    for CR_i = 1:length(Rm_arr)
        
        fluxForCR = optimalFluxMat(:, CR_i);
        Ra = Rm_arr(CR_i);
        
        if ~(Ra == 100 || Ra == 150 || Ra == 200)
            continue
        end
        
        chi = averagePorosMat(:, CR_i);
        
        Se = 230;
        S0_arr = CR_arr*Se./(CR_arr+1);
        
        
        actualFluxes= ~isnan(fluxForCR);
        if sum(actualFluxes) > 0
            Ra_leg{end+1} = ['Ra = ', num2str(Ra)];
            plot1 = plot(S0_arr(actualFluxes), fluxForCR(actualFluxes));
        end
    end
    hold off;
    box on;
    xlabel('$S_0$'); ylabel('$F_0$');
    title('Ice-ocean salt flux');
    
    xlim([0 70]);
    
    legend(Ra_leg, 'Location', 'eastoutside');
    
    set(hFluxCR,'Units','Inches');
    pos = get(hFluxCR,'Position');
    set(hFluxCR,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(gcf, '-dpdf', [figure_output_dir, 'saltFluxCRZoom.pdf'], printQuality);
    
end




if makeFluxDecompositionPlots
    
    
    hFluxDecomp = figure();
    set(hFluxDecomp, 'Position', [100 100 1600 1000]);
    m=2;
    n=2;
    
    denom = CR_mat.*optimalWidthMat;
    denom = optimalWidthMat;
    denom = 1;
    
    denom_str = '\mathcal{C}';
    denom_str = 'L';
    denom_str = '1';
    
    subplot(m, n, 1);
    plotField = optimalL2DiffusionFluxMat./(denom);
    makeOptimalPlot(Rm_arr, CR_arr, plotField, optimalFluxMat, 0, ...
        optimalExtrapolationFlag, '$Rm^*$', ['$F_{diffusion}/ (',denom_str,')$'], false, Rm_crit);
    %ylim([0 20]);
    
    subplot(m, n, 2);
    plotField = optimalL2FrameFluxMat./(denom);
    makeOptimalPlot(Rm_arr, CR_arr, plotField, optimalFluxMat, 0, ...
        optimalExtrapolationFlag, '$Rm^*$',['$F_{frame}/(',denom_str,')$'], false, Rm_crit);
    % ylim([0 100]);
    
    subplot(m, n, 3);
    plotField = optimalL2FluidFluxMat./(denom);
    makeOptimalPlot(Rm_arr, CR_arr, plotField, optimalFluxMat, 0, ...
        optimalExtrapolationFlag, '$Rm^*$',['$F_{fluid}/(',denom_str,')$'], true, Rm_crit);
    %ylim([0 50]);
    
end