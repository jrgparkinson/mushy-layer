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
    CR_MIN = 0.01;
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
    optimalhChiMat(CR_i, Rm_i) = state.hChi;
    optimalhPsiMat(CR_i, Rm_i) = state.hPsi;
    
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




hScalingArgsLog = figure();
set(hScalingArgsLog, 'Position', [100 100 1500 1000]);
m = 3;
n = 3;

axH = 0.18; axW = 0.18; offsetZ = 0.1; offsetX = 0.1;
axPos(1, :) = [offsetX             offsetZ+0.65 axW axH];
axPos(2, :) = [offsetX+axW*1.55    offsetZ+0.65 axW axH];
axPos(3, :) = [offsetX+2*axW*1.55  offsetZ+0.65 axW axH];
axPos(4, :) = [offsetX             offsetZ+0.3 axW axH];
axPos(5, :) = [offsetX+axW*1.55    offsetZ+0.3 axW axH];
axPos(6, :) = [offsetX+2*axW*1.55  offsetZ+0.3 axW axH];
axPos(7, :) = [offsetX             offsetZ axW axH];
axPos(8, :) = [offsetX+axW*1.55    offsetZ axW axH];
axPos(9, :) = [offsetX+2*axW*1.55  offsetZ axW axH];

%axSize = [0.32 0.8];
%axY = 0.12;


subplot(m, n, 1);

fluxPower = 1;

hold on;
CR_leg = {};
for CR_i = 1:length(CR_arr)
    
    flux = optimalFluxMat(CR_i,:);
    CR = CR_arr(CR_i);
    L = optimallMat(CR_i,:);
    h = optimalhMat(CR_i,:);
    psi = optimalPsiMat(CR_i,:);
    
    actualFluxes= ~isnan(flux).*(flux>1e-4);
    
    % Need to use this form
    RmStar = (Rm_arr-Rm_crit(CR_i))/Rm_crit(CR_i);
    
    %L = L(actualFluxes);
    %RmStar = Rm_arr/Rm_crit(CR_i);
    if (CR > 0.5) || CR <0.02 || isnan(CR)
        continue
    end
    
    
    x =(St/Le) * CR^(1)*RmStar;
    
    x =(CR^(1))*RmStar;
    
    x = CR^(5/3)*(Rm_arr-Rm_crit(CR_i));
    
    x = (Rm_arr-Rm_crit(CR_i));
    %x = CR^(5/3)*(Rm_arr);
    
    %x = Rm_arr;
    
    %x = sqrt(CR)*sqrt(RmStar)*5;
    
    %x = 5*psi.*h./L;
    
    y = flux.^(fluxPower);
    
    if sum(actualFluxes) > 0
        plot1 = plot(log10(x(actualFluxes==1) ),log10( y(actualFluxes==1)));
        CR_leg{end+1} = ['$\mathcal{C}=', num2str(CR), '$'];
    end
end

%scale_F_rm = plot([-1 2], [-2.9 0.1], '--');
scale_F_rm = plot([1 3], [-4.5 -2.5], '--');
scale_F_rm2 = plot([1 3], [-4.5 -3], '--');

hold off;
box on;

legend([scale_F_rm,scale_F_rm2], {'$Rm^{1}$', '$Rm^{3/4}$'}, 'Location', 'northwest');


xlabel(['$ log10 \mathcal{C} (Rm-\hat{Rm})/\hat{Rm}$']);
% xlabel('$ \sqrt{\mathcal{C} Rm^*} St $');
%xlabel('$\psi_a St h / L$');
xlabel('$log_{10}(Rm-\hat{Rm})$');
%xlabel('$log_{10}(Rm)$');
ylabel('$log_{10} F_O$');
set(gca, 'Position', axPos(1, :));





subplot(m, n, 2);


hold on;
%CR_leg = {};
for CR_i = 1:length(CR_arr)
    
    flux = optimalFluxMat(CR_i,:);
    CR = CR_arr(CR_i);
    L = optimallMat(CR_i,:);
    h = optimalhMat(CR_i,:);
    psi = optimalPsiMat(CR_i,:);
    
    actualFluxes= ~isnan(flux).*(flux>1e-5).*(psi>1e-2);
    
    
    % Need to use this form
    %RmStar = (Rm_arr-Rm_crit(CR_i))/Rm_crit(CR_i);
    RmStar = (Rm_arr)/Rm_crit(CR_i);
    
    if (CR > 0.5) || CR <0.02 || isnan(CR)
        continue
    end
    
    x =sqrt(CR/Le)*sqrt(RmStar);
    x = CR^2*(Rm_arr);
    x = Rm_arr; % - Rm_crit(CR_i);
    
    y = psi;
    
    if sum(actualFluxes) > 0
        plot1 = plot(log10(x(actualFluxes==1)) , log10(y(actualFluxes==1)));
        % CR_leg{end+1} = ['$\mathcal{C}=', num2str(CR), '$'];
    end
end

scale_psi_rm = plot([1.5 3], [-1 -0.25], '--');

hold off;
box on;

legend(scale_psi_rm, '$Rm^{1/2}$');
%xlabel(['$ log10 [ \sqrt{(\mathcal{C}/Le) Rm/\hat{Rm}}]$']);
% xlabel('$ \sqrt{\mathcal{C} Rm^*} St $');
%xlabel('$\psi_a St h / L$');
xlabel('$log_{10} [ \mathcal{C}^2 Rm ]$');
xlabel('$log_{10} [ Rm ]$');
ylabel('log10 max($\psi$)');

set(gca, 'Position', axPos(2,:));



subplot(m, n, 3);

hold on;
%CR_leg = {};
for CR_i = 1:length(CR_arr)
    
    flux = optimalFluxMat(CR_i,:);
    CR = CR_arr(CR_i);
    L = optimalWidthMat(CR_i,:);
    h = optimalhMat(CR_i,:);
    psi = optimalPsiMat(CR_i,:);
    
    actualFluxes= ~isnan(flux).*(flux>1e-5).*(psi>1e-2);
    
    
    % Need to use this form
    RmStar = (Rm_arr-Rm_crit(CR_i))/Rm_crit(CR_i);
    %RmStar = (Rm_arr)/Rm_crit(CR_i);
    
    if (CR > 0.5) || CR <0.02 || isnan(CR)
        continue
    end
    
    %(CR^(1/6)) also works fairly well
    x =1./(St*(RmStar));
    x = 1./(Rm_arr-Rm_crit(CR_i));
    
    y = L;
    
    if sum(actualFluxes) > 0
        plot(log10(x(actualFluxes==1)) ,log10( y(actualFluxes==1)));
        % CR_leg{end+1} = ['$\mathcal{C}=', num2str(CR), '$'];
    end
end
scale_L_Rm = plot([-3 0], [-0.75 0], '--');

hold off;
box on;




%xlabel(['$log_{10}  1 /  [ St ( Rm -\hat{Rm} )/ \hat{Rm} ]$']);
xlabel(['$log_{10}  [ 1/  (Rm-\hat{Rm}) ]$']);
% xlabel('$ \sqrt{\mathcal{C} Rm^*} St $');
%xlabel('$\psi_a St h / L$');
ylabel('$log_{10}(L)$');


legend(CR_leg, 'Location', 'eastoutside');
set(gca, 'Position', axPos(3,:));

newAxes = axes('Position', get(gca, 'Position'), 'visible', 'off');
legend(newAxes, [scale_L_Rm], {'$Rm^{1/4}$'}, 'Location', 'southeast');

set(gca, 'Position', axPos(3,:));



subplot(m, n, 4);

hold on;
%CR_leg = {};
for CR_i = 1:length(CR_arr)
    
    flux = optimalFluxMat(CR_i,:);
    CR = CR_arr(CR_i);
    L = optimalWidthMat(CR_i,:);
    h = optimalhMat(CR_i,:);
    %h = optimalChannelHMat(CR_i, :);
    psi = optimalPsiMat(CR_i,:);
    
    actualFluxes= ~isnan(flux).*(flux>1e-5).*(psi>1e-2);
    
    
    % Need to use this form
    RmStar = (Rm_arr-Rm_crit(CR_i))/Rm_crit(CR_i);
    RmStar = (Rm_arr)/Rm_crit(CR_i);
    
    if (CR > 0.5) || CR <0.02 || isnan(CR)
        continue
    end
    
    x = 1./Rm_arr;
    
    y = h;
    
    if sum(actualFluxes) > 0
        plot(log10(x(actualFluxes==1)) , log10(y(actualFluxes==1)));
        % CR_leg{end+1} = ['$\mathcal{C}=', num2str(CR), '$'];
    end
end

scale_h_Rm1 = plot([-3 -2], [-1.8 -1.55], '--');
scale_h_Rm2 = plot([-2.9 -2], [-1.7 -1.55], '--');

hold off;
box on;

legend([scale_h_Rm1,scale_h_Rm2], {'$Rm^{-1/4}$', '$Rm^{-1/6}$'}, 'Location', 'southeast');

xlabel(['$ log10 ( \mathcal{C}/St^2 )   [ Rm / \hat{Rm} ]$']);
% xlabel('$ \sqrt{\mathcal{C} Rm^*} St $');
xlabel('$log_{10}(1/Rm)$');
ylabel('$log_{10} (h)$');

set(gca, 'Position', axPos(4,:));

% subplot(m, n, 5);
% 
% hold on;
% %CR_leg = {};
% for CR_i = 1:length(CR_arr)
%     
%     flux = optimalFluxMat(CR_i,:);
%     CR = CR_arr(CR_i);
%     L = optimalWidthMat(CR_i,:);
%     %h = optimalhMat(CR_i,:);
%     h = optimalChannelHMat(CR_i, :);
%     psi = optimalPsiMat(CR_i,:);
%     H = optimalHMat(CR_i,:);
%     
%     actualFluxes= ~isnan(flux).*(flux>1e-5).*(psi>1e-2);
%     
%     
%     % Need to use this form
%     RmStar = (Rm_arr-Rm_crit(CR_i))/Rm_crit(CR_i);
%     RmStar = (Rm_arr)/Rm_crit(CR_i);
%     
%     if (CR > 0.5) || CR <0.02 || isnan(CR)
%         continue
%     end
%     
%     %x =CR./(St*(RmStar));
%     x = 1./Rm_arr;
%     y = H;
%     
%     if sum(actualFluxes) > 0
%         plot(log10(x(actualFluxes==1)) , log10(y(actualFluxes==1)));
%         % CR_leg{end+1} = ['$\mathcal{C}=', num2str(CR), '$'];
%     end
% end
% 
% scale_H_Rm = plot([-3 -1.5], [-0.88 -0.88], '--');
% hold off;
% box on;
% 
% legend(scale_H_Rm, '$Rm^0$');
% 
% %xlim([0 0.005]);
% %ylim([0 0.01]);
% 
% %xlabel(['$log_{10}  \mathcal{C} (1/St) 1/[  Rm/ \hat{Rm} ]$']);
% xlabel(['$log_{10} (1 /  Rm) $']);
% ylabel('$log_{10} H$');
% 
% set(gca, 'Position', axPos(5,:));
% 
subplot(m, n, 5);

hold on;
%CR_leg = {};
for CR_i = 1:length(CR_arr)
    
    flux = optimalFluxMat(CR_i,:);
    CR = CR_arr(CR_i);
    L = optimalWidthMat(CR_i,:);
    %h = optimalhMat(CR_i,:);
    h = optimalChannelHMat(CR_i, :);
    hChi = optimalhChiMat(CR_i,:);
    psi = optimalPsiMat(CR_i,:);
    H = optimalHMat(CR_i,:);
    
    actualFluxes= ~isnan(flux).*(flux>1e-5).*(psi>1e-2);
    
    
    % Need to use this form
    RmStar = (Rm_arr-Rm_crit(CR_i))/Rm_crit(CR_i);
    RmStar = (Rm_arr)/Rm_crit(CR_i);
    
    if (CR > 0.5) || CR <0.02 || isnan(CR)
        continue
    end
    
    %x =CR./(St*(RmStar));
    x = 1./Rm_arr;
    y = hChi;
    
    if sum(actualFluxes) > 0
        plot(log10(x(actualFluxes==1)) , log10(y(actualFluxes==1)));
        % CR_leg{end+1} = ['$\mathcal{C}=', num2str(CR), '$'];
    end
end

scale_H_Rm = plot([-3 -1], [-2.5 -2], '--');
hold off;
box on;

legend(scale_H_Rm, '$Rm^0$');

%xlim([0 0.005]);
%ylim([0 0.01]);

%xlabel(['$log_{10}  \mathcal{C} (1/St) 1/[  Rm/ \hat{Rm} ]$']);
xlabel(['$log_{10}  1/Rm $']);
ylabel('$log_{10} h_\chi$');

set(gca, 'Position', axPos(5,:));




subplot(m, n, 6);

hold on;
%CR_leg = {};
for CR_i = 1:length(CR_arr)
    
    flux = optimalFluxMat(CR_i,:);
    CR = CR_arr(CR_i);
    L = optimalWidthMat(CR_i,:);
    %h = optimalhMat(CR_i,:);
    h = optimalChannelHMat(CR_i, :);
    a = optimalaMat(CR_i, :);
    psi = optimalPsiMat(CR_i,:);
    Pi = averagePermMat(CR_i,:); %
    % Pi = optimalPermMat(CR_i,:); %pi at edge of confinement
    
    actualFluxes= ~isnan(flux).*(flux>1e-5).*(psi>1e-2);
    
    
    % Need to use this form
    RmStar = (Rm_arr-Rm_crit(CR_i))/Rm_crit(CR_i);
    RmStar = (Rm_arr)/Rm_crit(CR_i);
    
    if (CR > 0.5) || CR <0.02 || isnan(CR)
        continue
    end
    
    % For edge of confinement
    %x  = (CR.^(2.5)./(RmStar));
    
    x  = (CR./(RmStar));
    x = 1./Rm_arr;
    
    y = Pi;
    
    if sum(actualFluxes) > 0
        plot(log10(x(actualFluxes==1)) , log10(y(actualFluxes==1)));
        % CR_leg{end+1} = ['$\mathcal{C}=', num2str(CR), '$'];
    end
end

scale_Pi_rm1 = plot([-2.9 -1.8], [-0.7 -0.5], '--');
scale_Pi_rm2 = plot([-3 -2], [-1.75 -1.25], '--');

hold off;
box on;


legend([scale_Pi_rm1,scale_Pi_rm2], {'$Rm^{1/6}$','$Rm^{1/2}$'}, 'Location', 'southeast');

xlabel(['$log_{10}(\mathcal{C}/(Rm/\hat{Rm}))$']);
xlabel(['$log_{10}(1/Rm)$']);

ylabel('$log_{10} \bar{\Pi}$');
set(gca, 'Position', axPos(6,:));




subplot(m, n, 7);

hold on;
%CR_leg = {};
for CR_i = 1:length(CR_arr)
    
    flux = optimalFluxMat(CR_i,:);
    CR = CR_arr(CR_i);
    L = optimalWidthMat(CR_i,:);
    %h = optimalhMat(CR_i,:);
    h_channel = optimalChannelHMat(CR_i, :);
    a = optimalaMat(CR_i, :);
    psi = optimalPsiMat(CR_i,:);
    Pi = averagePermMat(CR_i,:); %
    l = optimallMat(CR_i, :);
    % Pi = optimalPermMat(CR_i,:); %pi at edge of confinement
    
    actualFluxes= ~isnan(flux).*(flux>1e-5).*(psi>1e-2);
    
    
    % Need to use this form
    RmStar = (Rm_arr-Rm_crit(CR_i))/Rm_crit(CR_i);
    RmStar = (Rm_arr)/Rm_crit(CR_i);
    
    if (CR > 0.5) || CR <0.02 || isnan(CR)
        continue
    end
    
    % For edge of confinement
    %x  = (CR.^(2.5)./(RmStar));
    
    x  = (CR.^(1/2)./(RmStar));
    x = (Rm_arr);
    
    y = h_channel;
    
    if sum(actualFluxes) > 0
        plot(log10(x(actualFluxes==1)) , log10(y(actualFluxes==1)));
        % CR_leg{end+1} = ['$\mathcal{C}=', num2str(CR), '$'];
    end
end


scale_h_Rm1 = plot([1.8 3], [-2.4 -2.1], '--');
scale_h_Rm2 = plot([1.8 3], [-2.5 -2.5], '--');

hold off;
box on;

legend([scale_h_Rm2, scale_h_Rm1], {'$Rm^{0}$', '$Rm^{1/4}$'});


%xlabel(['$log_{10} ( \mathcal{C}^{1/2}/(Rm/\hat{Rm})^{1})$']);
xlabel(['$log_{10} [ Rm ]$']);
ylabel('$log_{10} h_{channel}$');
set(gca, 'Position', axPos(7,:));



subplot(m, n, 8);

hold on;
%CR_leg = {};
for CR_i = 1:length(CR_arr)
    
    flux = optimalFluxMat(CR_i,:);
    CR = CR_arr(CR_i);
    L = optimalWidthMat(CR_i,:);
    %h = optimalhMat(CR_i,:);
    h = optimalChannelHMat(CR_i, :);
    a = optimalaMat(CR_i, :);
    psi = optimalPsiMat(CR_i,:);
    Pi = averagePermMat(CR_i,:); %
    l = optimallMat(CR_i, :);
    % Pi = optimalPermMat(CR_i,:); %pi at edge of confinement
    
    actualFluxes= ~isnan(flux).*(flux>1e-5).*(psi>1e-2);
    
    
    % Need to use this form
    RmStar = (Rm_arr-Rm_crit(CR_i))/Rm_crit(CR_i);
    RmStar = (Rm_arr)/Rm_crit(CR_i);
    
    if (CR > 0.5) || CR <0.02 || isnan(CR)
        continue
    end
    
    % For edge of confinement
    %x  = (CR.^(2.5)./(RmStar));
    
    x  = (CR.^(1/2)./(RmStar));
    x = 1./(Rm_arr);
    
    y = a;
    
    if sum(actualFluxes) > 0
        plot(log10(x(actualFluxes==1)) , log10(y(actualFluxes==1)));
        % CR_leg{end+1} = ['$\mathcal{C}=', num2str(CR), '$'];
    end
end

scale_a_Rm = plot([-2.75 -1.75], [-2.75 -2], '--');

hold off;
box on;

legend([scale_a_Rm], {'$Rm^{3/4}$'}, 'Location', 'northwest');

%xlabel(['$log_{10} ( \mathcal{C}^{1/2}/(Rm/\hat{Rm})^{1})$']);
xlabel(['$log_{10} [ 1/(Rm)]$']);
ylabel('$log_{10} a$');
set(gca, 'Position', axPos(8,:));



subplot(m, n, 9);

hold on;
%CR_leg = {};
for CR_i = 1:length(CR_arr)
    
    flux = optimalFluxMat(CR_i,:);
    CR = CR_arr(CR_i);
    L = optimalWidthMat(CR_i,:);
    %h = optimalhMat(CR_i,:);
    h = optimalChannelHMat(CR_i, :);
    hChi = optimalhChiMat(CR_i,:);
    hPsi = optimalhPsiMat(CR_i,:);
    psi = optimalPsiMat(CR_i,:);
    H = optimalHMat(CR_i,:);
    
    actualFluxes= ~isnan(flux).*(flux>1e-5).*(psi>1e-2);
    
    
    % Need to use this form
    RmStar = (Rm_arr-Rm_crit(CR_i))/Rm_crit(CR_i);
    RmStar = (Rm_arr)/Rm_crit(CR_i);
    
    if (CR > 0.5) || CR <0.02 || isnan(CR)
        continue
    end
    
    %x =CR./(St*(RmStar));
    x = 1./Rm_arr;
    y = hPsi;
    
    if sum(actualFluxes) > 0
        plot(log10(x(actualFluxes==1)) , log10(y(actualFluxes==1)));
        % CR_leg{end+1} = ['$\mathcal{C}=', num2str(CR), '$'];
    end
end

scale_H_Rm = plot([-3 -1], [-2.2 -1.7], '--');
hold off;
box on;

legend(scale_H_Rm, '$Rm^{1/4}$');

%xlim([0 0.005]);
%ylim([0 0.01]);

%xlabel(['$log_{10}  \mathcal{C} (1/St) 1/[  Rm/ \hat{Rm} ]$']);
xlabel(['$log_{10}  1/Rm $']);
ylabel('$log_{10} h_\psi$');

set(gca, 'Position', axPos(9,:));





%TODO: scaling in terms of C

hScalingArgsLogC = figure();
set(hScalingArgsLogC, 'Position', [100 100 1500 1000]);
m = 3;
n = 3;

axH = 0.18; axW = 0.18; offsetZ = 0.1; offsetX = 0.1;
axPos(1, :) = [offsetX             offsetZ+0.65 axW axH];
axPos(2, :) = [offsetX+axW*1.55    offsetZ+0.65 axW axH];
axPos(3, :) = [offsetX+2*axW*1.55  offsetZ+0.65 axW axH];
axPos(4, :) = [offsetX             offsetZ+0.3 axW axH];
axPos(5, :) = [offsetX+axW*1.55    offsetZ+0.3 axW axH];
axPos(6, :) = [offsetX+2*axW*1.55  offsetZ+0.3 axW axH];
axPos(7, :) = [offsetX             offsetZ axW axH];
axPos(8, :) = [offsetX+axW*1.55    offsetZ axW axH];
axPos(9, :) = [offsetX+2*axW*1.55  offsetZ axW axH];

%axSize = [0.32 0.8];
%axY = 0.12;


subplot(m, n, 1);

hold on
for Rm_i = 1:length(Rm_arr)
    
    flux = optimalFluxMat(:, Rm_i);
    
    L = optimallMat(:, Rm_i);
    h = optimalhMat(:, Rm_i);
    psi = optimalPsiMat(:, Rm_i);
    
    actualFluxes= ~isnan(flux).*(flux>1e-4).*(CR_arr.'<0.25).*(CR_arr.'>0.01);
    
    % Need to use this form
    Rm = Rm_arr(Rm_i);
    RmStar = (Rm-Rm_crit)./Rm_crit;
    
    %x = CR_arr.*RmStar;
    %x = CR_arr.*(Rm-Rm_crit);
    x = CR_arr;
    y = flux./(Rm-Rm_crit.');
    
    if sum(actualFluxes) > 0
        plot(log10(x(actualFluxes==1) ),log10( y(actualFluxes==1)));
        
    end
end

%scale_F_C = plot([0.5 2], [-4 -1.5], '--');
%scale_F_C = plot([0.5 2.5], [-4 -1], '--');
scale_F_C = plot([-1.8 -0.8], [-5.3 -4.3], '--');
scale_F_C2 = plot([-1.5 -0.5], [-5.2 -3.7], '--');
scale_F_C3 = plot([-1.5 -0.5], [-5.2 -3.2], '--');

hold off;
box on;

legend([scale_F_C, scale_F_C2, scale_F_C3], {'$\mathcal{C}$','$\mathcal{C}^{3/2}$','$\mathcal{C}^{2}$'}, 'Location','northwest');

%%xlabel(['$ log10 (\mathcal{C} (Rm-\hat{Rm})/\hat{Rm})$']);
%xlabel(['$ log10 (\mathcal{C} (Rm-\hat{Rm}))$']);
xlabel(['$ log10 [\mathcal{C} ]$']);
ylabel('$log_{10} (F_O/(Rm-\hat{Rm}))$');
set(gca, 'Position', axPos(1, :));


subplot(m, n, 2);


hold on
for Rm_i = 1:length(Rm_arr)
    
    flux = optimalFluxMat(:, Rm_i);
    
    L = optimallMat(:, Rm_i);
    h = optimalhMat(:, Rm_i);
    psi = optimalPsiMat(:, Rm_i);
    
    actualFluxes= ~isnan(flux).*(flux>1e-4).*(CR_arr.'<0.5).*(CR_arr.'>0.02);
    
    % Need to use this form
    Rm = Rm_arr(Rm_i);
    %RmStar = (Rm-Rm_crit)./Rm_crit;
    RmStar = (Rm)./Rm_crit;


    %x = (CR_arr./(RmStar.^(2)));
    %x = (CR_arr.^(1-4/3)/Rm.^2.5);
    x = (CR_arr);
    
    y = L.*(Rm-Rm_crit.').^(1/4);
    
    if sum(actualFluxes) > 0
        plot(log10(x(actualFluxes==1) ),log10( y(actualFluxes==1)));
        
    end
end

scale_L_C = plot([-1.5 -0.5], [ -1.58 -1.58], '--');

hold off;
box on;

legend([scale_L_C], {'$\mathcal{C}^{0}$'}, 'Location','northwest');

xlabel(['$ log10 [ \mathcal{C}]$']);
ylabel('$log_{10} [L  (Rm-\hat{Rm})^{1/4} ]$');
set(gca, 'Position', axPos(2, :));




subplot(m, n, 3);
Rm_leg = {};

hold on
for Rm_i = 1:length(Rm_arr)
    
    flux = optimalFluxMat(:, Rm_i);
    
    L = optimallMat(:, Rm_i);
    h = optimalhMat(:, Rm_i);
    psi = optimalPsiMat(:, Rm_i);
    
    actualFluxes= ~isnan(flux).*(flux>1e-4).*(CR_arr.'<0.5).*(CR_arr.'>0.02);
    
    % Need to use this form
    Rm = Rm_arr(Rm_i);
    %RmStar = (Rm-Rm_crit)./Rm_crit;
    RmStar = (Rm)./Rm_crit;


    
   % x = CR_arr.*(RmStar);
    x = CR_arr;
    
    y = psi./sqrt(Rm); %
     %
    if sum(actualFluxes) > 0
        plot(log10(x(actualFluxes==1) ),log10( y(actualFluxes==1)));
       Rm_leg{end+1} = ['$Rm=',num2str(Rm), '$']; 
    end
end

%plot([-1.5 0.5], [-1.8 -0.8], '--');
scale_psi_C = plot([-1.5 -0.5], [-3 -2], '--');
scale_psi_C2 = plot([-1.5 -0.5], [-3 -2.125], '--');

hold off;
box on;

legend(Rm_leg, 'Location', 'eastoutside');
set(gca, 'Position', axPos(3, :));



%xlabel(['$ log10 (\mathcal{C} \sqrt{(Rm)/\hat{Rm}})$']);

xlabel(['$ log10 (\mathcal{C})$']);
ylabel('$log_{10} (\psi/\sqrt{Rm})$');

newAxes = axes('Position', get(gca, 'Position'), 'visible', 'off');
legend(newAxes, [scale_psi_C,scale_psi_C2], {'$\mathcal{C}^{1}$','$\mathcal{C}^{7/8}$'}, 'Location','northwest');

set(gca, 'Position', axPos(3, :));





subplot(m, n, 4);

hold on
for Rm_i = 1:length(Rm_arr)
    
    flux = optimalFluxMat(:, Rm_i);
    
    L = optimallMat(:, Rm_i);
    h = optimalhMat(:, Rm_i);
    psi = optimalPsiMat(:, Rm_i);
    Pi = averagePermMat(:, Rm_i);
    
    actualFluxes= ~isnan(flux).*(flux>1e-4).*(CR_arr.'<0.5).*(CR_arr.'>0.02);
    
    % Need to use this form
    Rm = Rm_arr(Rm_i);
    %RmStar = (Rm-Rm_crit)./Rm_crit;
    RmStar = (Rm)./Rm_crit;


    %x = (CR_arr./(RmStar.^(2)));
    %x = (CR_arr.^(1-4/3)/Rm.^2.5);
    x = (CR_arr);
    
    y = Pi*Rm^(1/2);
    
    if sum(actualFluxes) > 0
        plot(log10(x(actualFluxes==1) ),log10( y(actualFluxes==1)));
        
    end
end

scale_Pi_C = plot([-1.5 -0.5], [0 0.5] , '--');

hold off;
box on;
legend([scale_Pi_C], {'$\mathcal{C}^{1/2}$'}, 'Location','northwest');

xlabel(['$ log10 ( \mathcal{C})$']);
ylabel('$log_{10} (\Pi Rm^{1/2})$');
set(gca, 'Position', axPos(4, :));





subplot(m, n, 5);

hold on
for Rm_i = 1:length(Rm_arr)
    
    flux = optimalFluxMat(:, Rm_i);
    
    L = optimallMat(:, Rm_i);
    h = optimalhMat(:, Rm_i);
   % h = optimalChannelHMat(:, Rm_i);
    psi = optimalPsiMat(:, Rm_i);
    Pi = averagePermMat(:, Rm_i);
    
    actualFluxes= ~isnan(flux).*(flux>1e-4).*(CR_arr.'<0.5).*(CR_arr.'>0.02);
    
    % Need to use this form
    Rm = Rm_arr(Rm_i);
    %RmStar = (Rm-Rm_crit)./Rm_crit;
    RmStar = (Rm)./Rm_crit;


    %x = (CR_arr./(RmStar.^(2)));
    %x = (CR_arr.^(1-4/3)/Rm.^2.5);
    x = (CR_arr);
    
    y = h.*(Rm.^(1/4));
    
    if sum(actualFluxes) > 0
        plot(log10(x(actualFluxes==1) ),log10( y(actualFluxes==1)));
        
    end
end

scale_h_C = plot([-1.4 -0.5 ], [-1 -0.85] , '--');
%scale_h_C2 = plot([-2 -1], [-1.55 -1.3] , '--');

hold off;
box on

legend([scale_h_C], {'$\mathcal{C}^{1/6}$'}, 'Location','northwest');

xlabel(['$ log10 ( \mathcal{C}  )$']);
ylabel('$log_{10} (h Rm^{1/4})$');
set(gca, 'Position', axPos(5, :));



subplot(m, n, 6);

hold on
for Rm_i = 1:length(Rm_arr)
    
    flux = optimalFluxMat(:, Rm_i);
    
    L = optimallMat(:, Rm_i);
    %h = optimalhMat(:, Rm_i);
    h_channel = optimalChannelHMat(:, Rm_i);
    psi = optimalPsiMat(:, Rm_i);
    Pi = averagePermMat(:, Rm_i);
    
    actualFluxes= ~isnan(flux).*(flux>1e-4).*(CR_arr.'<0.5).*(CR_arr.'>0.02);
    
    % Need to use this form
    Rm = Rm_arr(Rm_i);
    %RmStar = (Rm-Rm_crit)./Rm_crit;
    RmStar = (Rm)./Rm_crit;


    %x = (CR_arr./(RmStar.^(2)));
    %x = (CR_arr.^(1-4/3)/Rm.^2.5);
    x = (CR_arr);
    
    y = h_channel./(Rm.^(1/4));
    
    if sum(actualFluxes) > 0
        plot(log10(x(actualFluxes==1) ),log10( y(actualFluxes==1)));
        
    end
end

scale_h_C = plot([-1.5 -0.5 ], [-2.4 -1.4] , '--');
scale_h_C2 = plot([-1.5 -0.5], [-2.4 -1.525] , '--');

hold off;
box on

legend([scale_h_C,scale_h_C2], {'$\mathcal{C}$','$\mathcal{C}^{7/8}$'}, 'Location','northwest');

xlabel(['$ log10 ( \mathcal{C}  )$']);
ylabel('$log_{10} (h_{channel} / Rm^{1/4})$');
set(gca, 'Position', axPos(6, :));





subplot(m, n, 7);

hold on
for Rm_i = 1:length(Rm_arr)
    
    flux = optimalFluxMat(:, Rm_i);
    
    L = optimallMat(:, Rm_i);
    %h = optimalhMat(:, Rm_i);
    h_channel = optimalChannelHMat(:, Rm_i);
    psi = optimalPsiMat(:, Rm_i);
    Pi = averagePermMat(:, Rm_i);
    a = optimalaMat(:, Rm_i);
    
    actualFluxes= ~isnan(flux).*(flux>1e-4).*(CR_arr.'<0.5).*(CR_arr.'>0.02);
    
    % Need to use this form
    Rm = Rm_arr(Rm_i);
    %RmStar = (Rm-Rm_crit)./Rm_crit;
    RmStar = (Rm)./Rm_crit;


    %x = (CR_arr./(RmStar.^(2)));
    %x = (CR_arr.^(1-4/3)/Rm.^2.5);
    x = (CR_arr);
    
    y = a./(Rm.^(3/4));
    
    if sum(actualFluxes) > 0
        plot(log10(x(actualFluxes==1) ),log10( y(actualFluxes==1)));
        
    end
end

%scale_h_C = plot([-1.5 -0.5 ], [-2.4 -1.4] , '--');
%scale_h_C2 = plot([-1.5 -0.5], [-2.4 -1.525] , '--');

hold off;
box on

%legend([scale_h_C,scale_h_C2], {'$\mathcal{C}$','$\mathcal{C}^{7/8}$'}, 'Location','northwest');

xlabel(['$ log10 ( \mathcal{C}  )$']);
ylabel('$log_{10} (a / Rm^{3/4})$');
set(gca, 'Position', axPos(7, :));



subplot(m, n, 8);

hold on
for Rm_i = 1:length(Rm_arr)
    
    flux = optimalFluxMat(:, Rm_i);
    
    L = optimallMat(:, Rm_i);
    %h = optimalhMat(:, Rm_i);
    h_channel = optimalChannelHMat(:, Rm_i);
    hPsi = optimalhPsiMat(:, Rm_i);
    hChi = optimalhChiMat(:, Rm_i);
    psi = optimalPsiMat(:, Rm_i);
    Pi = averagePermMat(:, Rm_i);
    a = optimalaMat(:, Rm_i);
    
    actualFluxes= ~isnan(flux).*(flux>1e-4).*(CR_arr.'<0.5).*(CR_arr.'>0.02);
    
    % Need to use this form
    Rm = Rm_arr(Rm_i);
    %RmStar = (Rm-Rm_crit)./Rm_crit;
    RmStar = (Rm)./Rm_crit;


    %x = (CR_arr./(RmStar.^(2)));
    %x = (CR_arr.^(1-4/3)/Rm.^2.5);
    x = (CR_arr);
    
    y = hChi./(Rm.^(1/4));
    
    if sum(actualFluxes) > 0
        plot(log10(x(actualFluxes==1) ),log10( y(actualFluxes==1)));
        
    end
end

scale_h_C = plot([-1.4 -0.6 ], [-3.2 -2.4] , '--');
%scale_h_C2 = plot([-1.5 -0.5], [-2.4 -1.525] , '--');

hold off;
box on

legend([scale_h_C], {'$\mathcal{C}$'}, 'Location','northwest');

xlabel(['$ log10 ( \mathcal{C}  )$']);
ylabel('$log_{10} (h_\chi / Rm^{1/4})$');
set(gca, 'Position', axPos(8, :));



subplot(m, n, 9);

hold on
for Rm_i = 1:length(Rm_arr)
    
    flux = optimalFluxMat(:, Rm_i);
    
    L = optimallMat(:, Rm_i);
    %h = optimalhMat(:, Rm_i);
    h_channel = optimalChannelHMat(:, Rm_i);
    hPsi = optimalhPsiMat(:, Rm_i);
    hChi = optimalhChiMat(:, Rm_i);
    psi = optimalPsiMat(:, Rm_i);
    Pi = averagePermMat(:, Rm_i);
    a = optimalaMat(:, Rm_i);
    
    actualFluxes= ~isnan(flux).*(flux>1e-4).*(CR_arr.'<0.5).*(CR_arr.'>0.02);
    
    % Need to use this form
    Rm = Rm_arr(Rm_i);
    %RmStar = (Rm-Rm_crit)./Rm_crit;
    RmStar = (Rm)./Rm_crit;


    %x = (CR_arr./(RmStar.^(2)));
    %x = (CR_arr.^(1-4/3)/Rm.^2.5);
    x = (CR_arr);
    
    y = hPsi./(Rm.^(1/4));
    
    if sum(actualFluxes) > 0
        plot(log10(x(actualFluxes==1) ),log10( y(actualFluxes==1)));
        
    end
end

scale_h_C = plot([-1.4 -0.6 ], [-3 -2.2] , '--');
scale_h_C2 = plot([-1.5 -0.5], [-2.85 -2.1] , '--');
scale_h_C3 = plot([-1.5 -0.5], [-2.7 -2.2] , '--');

hold off;
box on

legend([scale_h_C, scale_h_C2,scale_h_C3], {'$\mathcal{C}$','$\mathcal{C}^{3/4}$','$\mathcal{C}^{1/2}$'}, 'Location','northwest');

xlabel(['$ log10 ( \mathcal{C}  )$']);
ylabel('$log_{10} (h_\psi / Rm^{1/4})$');
set(gca, 'Position', axPos(9, :));







% Plot relationships between scales
plotRel = false;

if plotRel

hRelationships = figure();
set(hRelationships, 'Position', [100 100 1500 1000]);
m = 3;
n = 3;

axH = 0.18; axW = 0.18; offsetZ = 0.1; offsetX = 0.1;
axPos(1, :) = [offsetX             offsetZ+0.65 axW axH];
axPos(2, :) = [offsetX+axW*1.55    offsetZ+0.65 axW axH];
axPos(3, :) = [offsetX+2*axW*1.55  offsetZ+0.65 axW axH];
axPos(4, :) = [offsetX             offsetZ+0.3 axW axH];
axPos(5, :) = [offsetX+axW*1.55    offsetZ+0.3 axW axH];
axPos(6, :) = [offsetX+2*axW*1.55  offsetZ+0.3 axW axH];
axPos(7, :) = [offsetX             offsetZ axW axH];
axPos(8, :) = [offsetX+axW*1.55    offsetZ axW axH];
axPos(9, :) = [offsetX+2*axW*1.55  offsetZ axW axH];

%axSize = [0.32 0.8];
%axY = 0.12;


subplot(m, n, 1);

hold on
for CR_i = 1:length(CR_arr)
    
    flux = optimalFluxMat(CR_i, :);
    
    L = optimallMat(CR_i, :);
    h = optimalhMat(CR_i, :);
    h_channel = optimalChannelHMat(CR_i, :);
    psi = optimalPsiMat(CR_i, :);
    CR = CR_arr(CR_i);
    a = optimalaMat(CR_i, :);
    
    actualFluxes= ~isnan(flux).*(flux>1e-4);
    
    % Need to use this form
    %Rm = Rm_arr(Rm_i);
    %RmStar = (Rm-Rm_crit)./Rm_crit;
    
   
    x = sqrt(1./(Rm_arr.*h_channel));
    y = a;
    
    if sum(actualFluxes) > 0
        plot(log10(x(actualFluxes==1) ),log10( y(actualFluxes==1)));
        
    end
end


hold off;
box on;

xlabel(['$H h_\chi / (Le Rm h_c)$']);
ylabel('$a^2$');
set(gca, 'Position', axPos(1, :));

end