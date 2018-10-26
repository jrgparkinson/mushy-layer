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

figure_output_dir = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/matlab/ConfirmationOfStatus/';

% load data:
load([data_dir, 'optimalVals.mat']);

printQuality = '-r100';

if strcmp(output_folder,'optimalStates-highRes-new/')
    RM_MAX = 800;
    CR_MAX = 20.0;
    CR_MIN = 0.005;
    Da = 1.0;
else
    RM_MAX = 1e5;
    CR_MAX = 40.0;
    Da = 5e-3;
end

Le = 200;
St = 5;


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
        Rm_arr(end+1) = thisRa*Da;
    end
    
end

Rm_arr = sort(Rm_arr); CR_arr = sort(CR_arr);

%Stefan number
St = 5;



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
optimalChannelHMat = optimalFluxMat;
optimalaMat = optimalFluxMat;
optimalhChiMat = optimalFluxMat;
optimalhPsiMat = optimalFluxMat;

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
    
    
if strcmp(output_folder,'optimalStates-Brinkman/')
    Rm = Da*Ra;
else
    Rm = Ra;
end
    
    Rm_i = find(Rm_arr == Rm);
    CR_i = find(CR_arr == CR);
    
    % Let's skip some things for now
    if Rm > RM_MAX || CR > CR_MAX || CR < CR_MIN
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
    optimalaMat(CR_i, Rm_i) = state.channelWidth;
    optimalhChiMat(CR_i, Rm_i) = NaN;
    optimalhPsiMat(CR_i, Rm_i) = NaN;
    
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
    optimalPsiMat(CR_i, Rm_i) = state.maxStreamfunctionInterior;
    optimalPsiMat(CR_i, Rm_i) = state.maxStreamfunctionMush;
    averagePorosMat(CR_i, Rm_i) = state.avPorosityMush;
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
psi_crit = CR_arr;

for CR_i = 1:length(CR_arr)
    CR = CR_arr(CR_i);
    
    fluxForCR = optimalFluxMat(CR_i, :);
    LForRm = optimalWidthMat(CR_i, :);
    actualFluxes= ~isnan(fluxForCR);
    %PiForRm = (optimalhPorosityMat(CR_i, :)).^3; % TODO check this is right
    PiForRm = (averagePermMat(CR_i, :)); % TODO check this is right
    hForRm = optimalhMat(CR_i, :); % TODO check this is right
    psi = optimalPsiMat(CR_i, :); % TODO check this is right
    
    %hForRm = optimalChannelHMat(CR_i, :);
    
    % Take smallest two fluxes and fit a straight line
    fluxesToFit = [];
    RmToFit = []; LToFit = []; PiToFit = []; hToFit = []; psiToFit = [];
    for i=1:length(fluxForCR)
        if ~isnan(fluxForCR(i)) && fluxForCR(i)/CR > 1e-3
            fluxesToFit(end+1) = fluxForCR(i);
            RmToFit(end+1) = Rm_arr(i);
            LToFit(end+1) = LForRm(i);
            PiToFit(end+1) = PiForRm(i);
            hToFit(end+1) = hForRm(i);
            psiToFit(end+1) = psi(i);
        end
    end
    
    if length(RmToFit) > 1
        Rm_crit(CR_i) = RmToFit(1) - fluxesToFit(1)*(RmToFit(2)-RmToFit(1))/(fluxesToFit(2)-fluxesToFit(1));
        LOpt_crit(CR_i) = LToFit(1) - fluxesToFit(1)*(LToFit(2)-LToFit(1))/(fluxesToFit(2)-fluxesToFit(1));
        Pi_crit(CR_i) = PiToFit(1) - fluxesToFit(1)*(PiToFit(2)-PiToFit(1))/(fluxesToFit(2)-fluxesToFit(1));
        h_crit(CR_i) = hToFit(1) - fluxesToFit(1)*(hToFit(2)-hToFit(1))/(fluxesToFit(2)-fluxesToFit(1));
        psi_crit(CR_i) = psiToFit(1) - fluxesToFit(1)*(psiToFit(2)-psiToFit(1))/(fluxesToFit(2)-fluxesToFit(1));
    elseif ~isempty(RmToFit)
        Rm_crit(CR_i) = RmToFit(1);
        LOpt_crit(CR_i) = LToFit(1);
        Pi_crit(CR_i) = PiToFit(1);
        h_crit(CR_i) = hToFit(1);
        psi_crit(CR_i) = psiToFit(1);
    else
        Rm_crit(CR_i) = NaN;
        LOpt_crit(CR_i) = NaN;
        Pi_crit(CR_i) = NaN;
        h_crit(CR_i) = NaN;
        psi_crit(CR_i) = NaN;
    end
    
    if psi_crit(CR_i)  < 0
        psi_crit(CR_i)  = NaN;
    end
    
    fprintf('CR = %1.3f, Rm_crit = %1.3f, L_crit = %1.3f\n', CR, Rm_crit(CR_i), LOpt_crit(CR_i));
    
    
end

noDodgy = 0*optimalExtrapolationFlag;



hStability = figure();
set(hStability, 'Position', [100 100 1600 900]);
m = 2; n = 3;
subplot(m, n, 1);
hold on
plot(log10(CR_arr), log10(1./Rm_crit));
plot([-2 1], [-2.5 -0.5], '--');
plot([-2 0], [-2.75 -0.75], '--');
plot([-2 -1], [-3 -1], '--');

hold off
box on;
xlabel('log$_{10}(\mathcal{C})$');
ylabel('log$_{10}(1/\hat{Rm})$');
legend({'Data',  '$\mathcal{C}^{2/3}$', '$\mathcal{C}$', '$\mathcal{C}^2$'}, 'Location', 'northwest');



subplot(m, n, 2);
hold on
plot(log10(CR_arr), log10(h_crit));
quarter = plot([-2 0], [-1.8 -1.2], '--');
half = plot([-2 0], [-1.8 -0.8], '--');
hold off;
box on;

legend([quarter half], {'$\mathcal{C}^{1/4}$', '$\mathcal{C}^{1/2}$'}, 'Location', 'southeast');

hold off;
xlabel('log$_{10}(\mathcal{C})$');
ylabel('log$_{10}(\hat{h})$');


subplot(m, n, 3);
hold on
plot(log10(CR_arr), log10(LOpt_crit));
quarter = plot([-2 0], [-0.5 0], '--');
sixth = plot([-2 1], [-0.5 0], '--');
hold off;
box on;

legend([quarter sixth], {'$\mathcal{C}^{1/4}$', '$\mathcal{C}^{1/6}$'}, 'Location', 'southeast');

hold off;
xlabel('log$_{10}(\mathcal{C})$');
ylabel('log$_{10}(\hat{L})$');


subplot(m, n, 4);
hold on
plot(log10(CR_arr), log10(Pi_crit));
threeQuarters = plot([-3 0], [-2 0], '--');
half = plot([-2 1], [-1.8 -0.3], '--');
hold off;
box on;

legend([threeQuarters half], {'$\mathcal{C}^{2/3}$', '$\mathcal{C}^{1/2}$'}, 'Location', 'southeast');

hold off;
xlabel('log$_{10}(\mathcal{C})$');
ylabel('log$_{10}(\hat{\Pi})$');


subplot(m, n, 5);
hold on
plot(log10(CR_arr(~isnan(psi_crit))), log10(psi_crit(~isnan(psi_crit))));
half = plot([-2 2], [-2.2 -0.2], '--');
quarter = plot([-2 2], [-2.1 -1.1], '--');
hold off;
box on;

legend([half quarter], {'$\mathcal{C}^{1/2}$', '$\mathcal{C}^{1/4}$'}, 'Location', 'southeast');

hold off;
xlabel('log$_{10}(\mathcal{C})$');
ylabel('log$_{10}(\hat{\psi})$');

subplot(m, n, 6);
yyaxis left;
plot(log10(CR_arr), h_crit.*Rm_crit.*Pi_crit./LOpt_crit);
ylabel('$\hat{Rm}  \hat{h} \hat{\Pi}/\hat{L}$');

yyaxis right;
hold on
plot(log10(CR_arr), St*Rm_crit.*h_crit.^2);
%plot(log10(CR_arr), (1/St)*psi_crit./h_crit.^2);
hold off;
ylabel('$St \hat{Rm}  \hat{h}^2$');

xlabel('log$_{10}(\mathcal{C})$');


 set(hStability,'Units','Inches');
 pos = get(hStability,'Position');
 set(hStability,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
 print(gcf, '-dpdf', [figure_output_dir, 'stabilityPlots.pdf'], printQuality);


 
 
 
 
 
