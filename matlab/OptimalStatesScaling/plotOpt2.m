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

figure_output_dir = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/matlab/Scaling/';

% load data:
load([data_dir, 'optimalVals.mat']);

printQuality = '-r100';
makeStandardPlots = true;
makeOceanPlots  = false;
makeDetailedPlots = true;
makeStabilityPlots = false;
makeFluxDecompositionPlots = false;

scaleDetailedPlotsWithRaCrit = true;


RA_MAX = 450;
CR_MAX = 2.0;


% Get Rayleigh numbers and concentration ratios

% Get all concentration ratios and rayleigh numbers
k = keys(optimalVals);
CR_arr = [];
Ra_arr = [];
for CR_i = 1:length(k)
    keyPair = k(CR_i);
    keyPair = keyPair{1};
    temp = keyPair(1);  thisC  = temp{1};
    temp = keyPair(2);  thisRa = temp{1};
    
    
    if  ~ismember(thisC, CR_arr)
        CR_arr(end+1) = thisC;
    end
    
    if ~ismember(thisRa, Ra_arr)
        Ra_arr(end+1) = thisRa;
    end
    
end

Ra_arr = sort(Ra_arr); CR_arr = sort(CR_arr);

%Stefan number
St = 5;



% Now, for each Ra, C pair, get matrices of optimal values

optimalFluxMat = NaN*ones(length(CR_arr), length(Ra_arr));
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
    
    Ra_i = find(Ra_arr == Ra);
    CR_i = find(CR_arr == CR);
    
    % Let's skip some things for now
    if Ra > RA_MAX || CR > CR_MAX
        continue
    end
    
    optimalFluxMat(CR_i, Ra_i) = state.flux;
    optimalWidthMat(CR_i, Ra_i) = state.fullWidth;
    
    if isfield(state, 'criticalWidth')  
    optimalCritWidthMat(CR_i, Ra_i) = state.criticalWidth;
    else
         optimalCritWidthMat(CR_i, Ra_i) = NaN;
    end
    
    if isfield(state, 'flagCritWidth')
        optimalFlagCritWidth(CR_i, Ra_i) = state.flagCritWidth;
    else
        optimalFlagCritWidth(CR_i, Ra_i) = false;
    end
    
    optimallMat(CR_i, Ra_i) = state.smallLScale;
    optimalHMat(CR_i, Ra_i) = state.H;
    optimalhMat(CR_i, Ra_i) = state.h;
    optimalChannelHMat(CR_i, Ra_i) = state.channelHeight;
    
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
    
    optimalPorosMat(CR_i, Ra_i) = state.maxPorosityInterior;
    optimalPermMat(CR_i, Ra_i) = state.maxPermeabilityInterior;
    optimalPsiMat(CR_i, Ra_i) = state.maxStreamfunctionMush;
    averagePorosMat(CR_i, Ra_i) = state.avPorosityInterior;
    averagePermMat(CR_i, Ra_i) = state.avPermeabilityInterior;
    
    if isfield(state, 'L2FsVertDiffusion')
    optimalL2DiffusionFluxMat(CR_i, Ra_i) = state.L2FsVertDiffusion;
    optimalL2FrameFluxMat(CR_i, Ra_i) = state.L2FsVertFrame;
    optimalL2FluidFluxMat(CR_i, Ra_i) = state.L2FsVertFluid;
    
    optimalL0DiffusionFluxMat(CR_i, Ra_i) = state.L0FsVertDiffusion;
    optimalL0FrameFluxMat(CR_i, Ra_i) = state.L0FsVertFrame;
    optimalL0FluidFluxMat(CR_i, Ra_i) = state.L0FsVertFluid;
    else
        optimalL2DiffusionFluxMat(CR_i, Ra_i) = NaN;
    optimalL2FrameFluxMat(CR_i, Ra_i) = NaN;
    optimalL2FluidFluxMat(CR_i, Ra_i) = NaN;
    
    optimalL0DiffusionFluxMat(CR_i, Ra_i) = NaN;
    optimalL0FrameFluxMat(CR_i, Ra_i) = NaN;
    optimalL0FluidFluxMat(CR_i, Ra_i) = NaN;
    end
    
    % Mush:
    
        optimalSaltAdvection(CR_i, Ra_i) = state.saltAdvectionMush;
        optimalSaltDiffusion(CR_i, Ra_i) = state.saltDiffusionMush;
        optimalSolidSaltFrame(CR_i, Ra_i) = state.solidSalinityFrameMush;
        optimalLiquidSaltFrame(CR_i, Ra_i) = state.liquidSalinityFrameMush;
    
        optimalVorticityDissipation(CR_i, Ra_i) = state.vorticityDiffusionMush;
        optimalBaroclinicTorque(CR_i, Ra_i) = state.baroclinicTorqueMush;
        optimalVorticityPermeability(CR_i, Ra_i) = state.vorticityPermeabilityMush;
    
        optimalHeatAdvection(CR_i, Ra_i) = state.heatAdvectionMush;
        optimalHeatDiffusion(CR_i, Ra_i) = state.heatDiffusionMush;
        optimalLatentHeat(CR_i, Ra_i) = state.latentHeatMush;
        optimalTempFrame(CR_i, Ra_i) = state.TFrameAdvectionMush;
        
        if isfield(state, 'extrapolatedMax') && ~isnan(state.extrapolatedMax)
            optimalExtrapolationFlag(CR_i, Ra_i) = state.extrapolatedMax;
        else
            optimalExtrapolationFlag(CR_i, Ra_i) = false;
        end
    
end

% Let's calculate critical rayleigh numbers properly, by extrapolating
% fluxes to zero

Ra_crit = CR_arr;
LOpt_crit = CR_arr;

for CR_i = 1:length(CR_arr)
    CR = CR_arr(CR_i);
    
    fluxForCR = optimalFluxMat(CR_i, :);
    LForRa = optimalWidthMat(CR_i, :);
    actualFluxes= ~isnan(fluxForCR);
    
    
    % Take smallest two fluxes and fit a straight line
    fluxesToFit = [];
    RaToFit = []; LToFit = [];
    for i=1:length(fluxForCR)
        if ~isnan(fluxForCR(i)) && fluxForCR(i)/CR > 1e-3
            fluxesToFit(end+1) = fluxForCR(i);
            RaToFit(end+1) = Ra_arr(i);
            LToFit(end+1) = LForRa(i);
        end
    end
    
    if length(RaToFit) > 1
        Ra_crit(CR_i) = RaToFit(1) - fluxesToFit(1)*(RaToFit(2)-RaToFit(1))/(fluxesToFit(2)-fluxesToFit(1));
        LOpt_crit(CR_i) = LToFit(1) - fluxesToFit(1)*(LToFit(2)-LToFit(1))/(fluxesToFit(2)-fluxesToFit(1));
    elseif ~isempty(RaToFit)
        Ra_crit(CR_i) = RaToFit(1);
        LOpt_crit(CR_i) = LToFit(1);
    else
        Ra_crit(CR_i) = NaN;
        LOpt_crit(CR_i) = NaN;
    end
    
    fprintf('CR = %1.3f, Ra_crit = %1.3f, L_crit = %1.3f\n', CR, Ra_crit(CR_i), LOpt_crit(CR_i));
    
    %Ra_crit(CR_i) = 30 + 2/CR;
    %Ra_crit(CR_i) = 0;
end

%figure(); semilogx(CR_arr, 1./L_crit); xlabel('CR'); ylabel('$a_c$');

%T = table(CR_arr, Ra_crit)

%Ra_crit


h = figure();
set(h, 'Position', [100 100 1300 1000]);

%m = 2; n = 2;

%subplot(m, n, 1);
hold on;
CR_leg = {};
for CR_i = 1:length(CR_arr)
    
    fluxForCR = optimalFluxMat(CR_i, :);
   % Ra = Ra_arr(Ra_i);
    CR = CR_arr(CR_i);
    psi = optimalPsiMat(CR_i, :);
    L = optimalWidthMat(CR_i, :);
    
    chi = optimalPorosMat(CR_i, :);
    
    %chi = averagePorosMat(CR_i, :);
    
    
    H = optimalHMat(CR_i, :);
     H = optimalhMat(CR_i, :);
    
    actualFluxes= ~isnan(fluxForCR);
    
    L = L(actualFluxes);
    chi = chi(actualFluxes);
     H = H(actualFluxes);
     psi = psi(actualFluxes);
    
    phi = 1 - chi;
    Ra_plot = Ra_arr;
    Ra_star = (Ra_arr-Ra_crit(CR_i))/Ra_crit(CR_i);
    
    Ra_star = Ra_star(actualFluxes);
    Ra_plot = Ra_plot(actualFluxes);
    
    F = fluxForCR(actualFluxes);
    
    if sum(actualFluxes) > 0
        

        %p1 = plot(Ra_star.^0.5, psi/CR^0.5);
        %p2 = plot(log(Ra_star), log(F/CR));
        power = 1/4; %1/4 or maybe 1/3
        %p3 = plot(Ra_star.^(-power), L/(CR^(power)));
        %p3 = plot(log10(Ra_star), log10( L/CR^(1/4)) );
        p5 = plot(Ra_star.^(-0.5), H/CR^0.25);
        
        
        %p4 = plot(log10(1./Ra_plot), log10(chi/(CR^0.5)));
         
        
        CR_leg{end+1} = ['$\mathcal{C} = ', num2str(CR),'$'];
    end
end
hold off;
box on;
%xlabel('$\bar{\phi} \mathcal{C}/L$');
%xlabel('$\mathcal{C} \mathcal{R}^{1/2}$');
%ylabel('$\bar{\psi}/L^2$');
%title('Streamfunction');
legend(CR_leg, 'Location', 'eastoutside');

%ylim([0 0.3]);

%ylim([0 10]);