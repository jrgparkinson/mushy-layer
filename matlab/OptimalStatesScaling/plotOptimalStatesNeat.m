% Plot optimal values

function plotOptimalStatesNeat(output_folder, savePlots)

close all;

if nargin < 2
    savePlots = true;
end

if nargin < 1
    output_folder = 'optimalStates-Brinkman-lowRes/';
    output_folder = 'optimalStates-highRes-new/';
end



%set(groot, 'DefaultAxesColorOrder', [1 0 0; 0 1 0; 0 0 1], ...
%    'DefaultAxesLineStyleOrder', 'x-|x--|x:');
set(groot, 'DefaultAxesLineStyleOrder', 'x-|x--|x:');
%set(groot, 'DefaultAxesLineStyleOrder',{'--',':'})

data_dir = getDataDir(output_folder);

% figure_output_dir = data_dir; %'/home/parkinsonjl/convection-in-sea-ice/MushyLayer/matlab/ConfirmationOfStatus/';
figure_output_dir = '/home/parkinsonjl/convection-in-sea-ice/figures/optimalStates/';
figure_output_dir = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/matlab/OptimalStatesScaling/';

% load data:
load([data_dir, 'optimalVals.mat']);

printQuality = '-r100';

% Whether to plot how all the different metrics scale
doScalingArgPlots= true;

% Different possible scalings
scaling = getExponents(0); % old scaling
scaling = getExponents(1); % new dodgy scaling
scaling = getExponents(2); % new better scaling

if strcmp(output_folder,'optimalStates-highRes-new/')
    RM_MAX = 800;
    CR_MAX = 20.0;
    CR_MIN = 0.005;
    Da = 1.0;
else
    RM_MAX = 1e5;
    CR_MAX = 40.0;
    Da = 5e-3;
    CR_MIN = 0.1;
end

Le = 200;
St = 5;
Pr = 10;
Vf = 1;

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
    
    
    if strcmp(output_folder,'optimalStates-Brinkman-lowRes/')
        Rm = Da*Ra*Pr*Pr;
    else
        Rm = Ra;
    end
    
    Rm_i = find(Rm_arr == Rm);
    CR_i = find(CR_arr == CR);
    
    % Let's skip some things for now
    if Rm > RM_MAX || CR > CR_MAX || CR < CR_MIN
        continue
    end
    
    if CR==0.25 && Rm > 300 || CR==0.08 && Rm > 450
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
    %hForRm = optimalhMat(CR_i, :); % TODO check this is right
    hForRm = optimalChannelHMat(CR_i, :); % TODO check this is right
    
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
    
    if  CR_arr(CR_i) < 0.01
        Rm_crit(CR_i) = NaN;
        LOpt_crit(CR_i) = NaN;
        Pi_crit(CR_i) = NaN;
        h_crit(CR_i) = NaN;
        psi_crit(CR_i) = NaN;
    end
    fprintf('CR = %1.3f, Rm_crit = %1.3f, L_crit = %1.3f\n', CR, Rm_crit(CR_i), LOpt_crit(CR_i));
    
    
end

noDodgy = 0*optimalExtrapolationFlag;


unscaledPlot = false;
if unscaledPlot
    
    hFluxRmCR = figure();
    set(hFluxRmCR, 'Position', [100 100 1000 800]);
    m = 2;
    n = 2;
    
    axSize = [0.3 0.35];
    axPos(1, :) = [0.1 0.6 axSize(1) axSize(2)];
    axPos(2, :) = [0.5 0.6 axSize(1) axSize(2)];
    axPos(3, :) = [0.1 0.1 axSize(1) axSize(2)];
    axPos(4, :) = [0.5 0.1 axSize(1) axSize(2)];
    
    %axY = 0.12;
    
    
    subplot(m, n, 1);
    
    CR_mat = repmat(CR_arr, length(Rm_arr), 1).';
    plotField = optimalFluxMat./(CR_mat.^0);
    %plotField = log10(optimalFluxMat);
    plotField = optimalFluxMat;
    
    plotField([find(CR_arr==20):length(CR_arr)], :) = NaN;
    CR_arr(CR_arr==20) = NaN;
    
    
    CR_leg = {};
    hold on;
    for CR_i = 1:length(CR_arr)
        pf = squeeze(plotField(CR_i, :));
        valid = ~isnan(pf);
        
        CR = CR_arr(CR_i);
        if CR == 10.0
            valid(Rm_arr==60) = 0; % these are dodgy
            valid(Rm_arr==70) = 0;
        end
        
        plot(Rm_arr(valid), pf(valid));
        CR_leg{end+1} = ['$\mathcal{C}=',num2str(CR),'$'];
    end
    hold off;
    xlabel('$Rm_S$'); ylabel('$F_O$');
    title('(a)');
    box on;
    
    ax_a = gca;
    
    xlim([0 250]);
    ylim([0 1.5]);
    
    set(gca, 'Position', axPos(1, :));
    
    subplot(m, n, 2);
    % Only plot small C values
    plotField = optimalFluxMat;
    plotField([find(CR_arr==1.0):length(CR_arr)], :) = NaN;
    
    
    hold on;
    for CR_i = 1:length(CR_arr)
        pf = squeeze(plotField(CR_i, :));
        valid = ~isnan(pf);
        
        CR = CR_arr(CR_i);
        %if CR == 10.0
        %    valid(Rm_arr==60) = 0; % these are dodgy
        %    valid(Rm_arr==70) = 0;
        %end
        
        plot(Rm_arr(valid), pf(valid));
        
    end
    hold off;
    xlabel('$Rm_S$'); ylabel('$F_O$');
    title('(b)');
    box on;
    
    ylim([0 0.05]);
    
    % Want legend on right hand plot but for both sets of data
    lgd = legend(ax_a, CR_leg, 'location', 'eastoutside');
    thisAxPos = axPos(2, :);
    lgdPos(1) = thisAxPos(1) + thisAxPos(3)*1.2;
    lgdPos(2) = thisAxPos(2);
    lgdPos(3) = 0.1;
    lgdPos(4) = thisAxPos(4);
    set(lgd, 'position', lgdPos);
    
    set(gca, 'Position', axPos(2, :));
    
    
    subplot(m, n, 3);
    %plotField = log10(optimalFluxMat);
    plotField = optimalFluxMat;
    ReducedRm = Rm_arr;
    
    %RemoveRm = [450, 500, 45, 60, 125,170,70,63,40,175,90,120,50,80,140];
    RemoveRm = [500,450,  170,140, 125, 120, 100, 63, 90, 70,  65, 40,35,  30];
    KeepRm = [35, 40, 45, 60, 100, 250 300 350 400, 500];
    for i=1:length(RemoveRm)
        ReducedRm(find(Rm_arr==RemoveRm(i))) = NaN;
    end
    
    % 'log$_{10}$($\mathcal{C}$)'
    
    hold on;
    for Rm_i = 1:length(Rm_arr)
        pf = squeeze(plotField(:, Rm_i));
        valid = ~isnan(pf);
        
        Rm = Rm_arr(Rm_i);
        
        if sum(find(KeepRm==Rm)) > 0
            
            plot(CR_arr(valid), pf(valid));
            
        end
        
    end
    hold off;
    xlabel('$\mathcal{C}$'); ylabel('$F_O$');
    title('(c)');
    box on;
    
    set(gca, 'Position', axPos(3, :));
    
    subplot(m, n, 4);
    
    plotField = optimalFluxMat;
    
    Rm_leg = {};
    hold on;
    for Rm_i = 1:length(Rm_arr)
        pf = squeeze(plotField(:, Rm_i));
        valid = ~isnan(pf);
        
        Rm = Rm_arr(Rm_i);
        
        if sum(find(KeepRm==Rm)) > 0
            
            if Rm == 100
                valid(CR_arr < 0.25) = 0;
            end
            
            plot(CR_arr(valid), pf(valid));
            if sum(valid) > 0
                Rm_leg{end+1} = ['$Rm_s = ',num2str(Rm),'$'];
            end
        end
        
    end
    hold off;
    xlabel('$\mathcal{C}$'); ylabel('$F_O$');
    title('(d)');
    box on;
    
    ylim([0,0.04]);
    xlim([0,0.25]);
    
    lgd = legend(Rm_leg, 'Location', 'eastoutside');
    %columnlegend(2,Rm_leg,'location', 'eastoutside'); %
    
    
    set(gca, 'Position', axPos(4, :));
    
    if savePlots
        
        set(hFluxRmCR,'Units','Inches');
        pos = get(hFluxRmCR,'Position');
        set(hFluxRmCR,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(gcf, '-dpdf', [figure_output_dir, 'neatFluxesRmCR.pdf'], printQuality);
        
    end
    
    
    
end





doScaledPlots = false;
if doScaledPlots
    
    hFluxRmScaled = figure();
    set(hFluxRmScaled, 'Position', [100 100 1500 600]);
    m = 1;
    n = 2;
    
    axSize = [0.32 0.8];
    axY = 0.12;
    
    
    subplot(m, n, 1);
    
    CR_mat = repmat(CR_arr, length(Rm_arr), 1).';
    plotField = optimalFluxMat./(CR_mat.^0);
    %plotField = log10(optimalFluxMat);
    plotField = optimalFluxMat;
    
    plotField([find(CR_arr==20):length(CR_arr)], :) = NaN;
    CR_arr(CR_arr==20) = NaN;
    makeOptimalPlot(Rm_arr, CR_arr, plotField, optimalFluxMat, 0, ...
        noDodgy, '$Rm_S^*$', '$F_O$   ', false, Rm_crit);
    
    set(gca, 'Position', [0.06 axY axSize]);
    xlim([0 5]);
    
    subplot(m, n, 2);
    % Only plot small C values
    plotField = optimalFluxMat;
    plotField([find(CR_arr==0.25):length(CR_arr)], :) = NaN;
    makeOptimalPlot(Rm_arr, CR_arr, plotField, optimalFluxMat, 0, ...
        noDodgy, '$Rm_S^*$', '$F_O$   ', true, Rm_crit);
    %xlim([0 250]);
    xlim([0 2]);
    
    set(gca, 'Position', [0.55 axY axSize]);
    
    
    %if makePlots
    % set(hFluxRmScaled,'Units','Inches');
    % pos = get(hFluxRmScaled,'Position');
    % set(hFluxRmScaled,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    % print(gcf, '-dpdf', [figure_output_dir, 'neatFluxesRmscaled.pdf'], printQuality);
    %end
end



if doScalingArgPlots
    
    hScalingArgs = figure();
    set(hScalingArgs, 'Position', [100 100 1500 900]);
    m = 3;
    n = 3;
    
    axH = 0.18; axW = 0.18; offsetZ = 0.1; offsetX = 0.1;
    axPos(1, :) = [offsetX             offsetZ+0.67 axW axH];
    axPos(2, :) = [offsetX+axW*1.55    offsetZ+0.67 axW axH];
    axPos(3, :) = [offsetX+2*axW*1.55  offsetZ+0.67 axW axH];
    axPos(4, :) = [offsetX             offsetZ+0.34 axW axH];
    axPos(5, :) = [offsetX+axW*1.55    offsetZ+0.34 axW axH];
    axPos(6, :) = [offsetX+2*axW*1.55  offsetZ+0.34 axW axH];
    axPos(7, :) = [offsetX             offsetZ axW axH];
    axPos(8, :) = [offsetX+axW*1.55    offsetZ axW axH];
    axPos(9, :) = [offsetX+2*axW*1.55  offsetZ axW axH];
    
    %axSize = [0.32 0.8];
    %axY = 0.12;
    
    
    subplot(m, n, 1);
    
    plotScaling('$F_O/\mathcal{C}$', scaling.F, optimalFluxMat, Rm_arr, Rm_crit, CR_arr, St, Le, Vf);
    
    title('(a)');
    set(gca, 'Position', axPos(1, :));
    
    
    
    
    
    subplot(m, n, 2);
    
    plotScaling('max($\psi$)', scaling.psi, optimalPsiMat, Rm_arr, [], CR_arr, St, Le, Vf);
    
    
    
    title('(b)');
    set(gca, 'Position', axPos(2,:));
    
    
    
    subplot(m, n, 3);
    
    
    plotScaling('$L$', scaling.L, optimalWidthMat, Rm_arr, Rm_crit, CR_arr, St, Le, Vf);
    %plotScaling('$L$', scaling.L, optimalWidthMat, Rm_arr, [], CR_arr, St, Le, Vf);
    
    title('(c)');
    set(gca, 'Position', axPos(3,:));
    
    
    subplot(m, n, 4);
    
    
    plotScaling('$h_c$', scaling.h_c, optimalChannelHMat, Rm_arr, [], CR_arr, St, Le, Vf);
    title('(d)');
    set(gca, 'Position', axPos(4,:));
    
    subplot(m, n, 5);
    
    
    
    [CR_plots, CR_leg] = plotScaling('$h_\chi$', scaling.h_chi, optimalhPsiMat, Rm_arr, [], CR_arr, St, Le, Vf);
    
    title('(e)');
    set(gca, 'Position', axPos(5,:));
    
    
    
    lgd =  legend(CR_plots, CR_leg, 'Location', 'west');
    lgdPos = axPos(6,:);
    lgdPos(3) = 0.1;
    lgdPos(2) = lgdPos(2) - 0.1;
    set(lgd, 'position',   lgdPos);
    
    
    
    
    subplot(m, n, 7);
    
    
    plotScaling('$a$', scaling.a, optimalaMat, Rm_arr, [], CR_arr, St, Le, Vf);
    
    title('(f)');
    set(gca, 'Position', axPos(7,:));
    
    
    
    
    
    
    
    
    
    
    subplot(m, n, 8);
    
    
    scaling.H.S = 0;
    scaling.H.V = 0;
    scaling.H.R = 1;
    
    plotScaling('$H \, \mathcal{S} \, V_f$', scaling.H, optimalHMat.*(St*Vf), Rm_arr, [], CR_arr, St, Le, Vf);
    
    ylim([0 1]);
    title('(g)');
    set(gca, 'Position', axPos(8,:));
    
    
    
    if savePlots
        set(hScalingArgs,'Units','Inches');
        pos = get(hScalingArgs,'Position');
        set(hScalingArgs,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(gcf, '-dpdf', [figure_output_dir, 'empiricalScaling.pdf'], printQuality);
        print(gcf, '-depsc', [figure_output_dir, 'empiricalScaling.eps'], printQuality);
    end
end



doScalingArgPlotsTalk= false;
if doScalingArgPlotsTalk
    
    hScalingArgsTalk = figure();
    set(hScalingArgsTalk, 'Position', [100 100 1200 360]);
    m = 1;
    n = 3;
    
    axW = 0.19; offsetZ = 0.27; offsetX = 0.08; axH = 1-offsetZ-0.05;
    expon=1.5;
    
    axPos(1, :) = [offsetX                offsetZ axW axH];
    axPos(2, :) = [offsetX+axW*expon    offsetZ axW axH];
    axPos(3, :) = [offsetX+2*axW*expon  offsetZ axW axH];
    
    subplot(m, n, 1);
    
    hold on;
    CR_leg = {};
    for CR_i = 1:length(CR_arr)
        
        flux = optimalFluxMat(CR_i,:);
        CR = CR_arr(CR_i);
        L = optimallMat(CR_i,:);
        h = optimalhMat(CR_i,:);
        psi = optimalPsiMat(CR_i,:);
        
        actualFluxes= ~isnan(flux);
        
        % Need to use this form
        RmStar = (Rm_arr-Rm_crit(CR_i))/Rm_crit(CR_i);
        
        %L = L(actualFluxes);
        %RmStar = Rm_arr/Rm_crit(CR_i);
        if (CR > 0.5) || CR <0.02 || isnan(CR)
            continue
        end
        
        
        x = (Rm_arr-Rm_crit(CR_i)).^(3/4) / (St^(1/2)*sqrt(Le));
        
        y = flux/(CR); %.^(fluxPower);
        
        if sum(actualFluxes) > 0
            plot1 = plot(x(actualFluxes) , y(actualFluxes));
            CR_leg{end+1} = ['$\mathcal{C}=', num2str(CR), '$'];
        end
    end
    hold off;
    box on;
    
    ylim([0 0.5]);
    
    xlabel(['$   \hat{Rm}_S^{3/4} / (St Le)^{1/2}$']);
    ylabel('$F_O/\mathcal{C}$');
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
        
        x =CR*sqrt(Rm_arr)/(St*sqrt(Le));
        
        y = psi;
        
        if sum(actualFluxes) > 0
            plot1 = plot(x(actualFluxes==1) , y(actualFluxes==1));
            % CR_leg{end+1} = ['$\mathcal{C}=', num2str(CR), '$'];
        end
    end
    hold off;
    box on;
    
    xlim([0 1]);
    xlim([0 0.1]);
    % ylim([0 0.1]);
    
    xlabel(['$ \mathcal{C} Rm_S^{1/2}/[St Le^{1/2}]$']);
    % xlabel('$ \sqrt{\mathcal{C} Rm^*} St $');
    %xlabel('$\psi_a St h / L$');
    ylabel('max($\psi_{channel}$)');
    % title('(b)');
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
        x =1./(sqrt(St)*((Rm_arr-Rm_crit(CR_i)).^(1/4)));
        
        y = L;
        
        if sum(actualFluxes) > 0
            plot(x(actualFluxes==1) , y(actualFluxes==1));
            % CR_leg{end+1} = ['$\mathcal{C}=', num2str(CR), '$'];
        end
    end
    hold off;
    box on;
    
    xlim([0.05 0.45]);
    ylim([0.05 0.6]);
    
    xlabel(['$ 1 /  [ (St V_f)^{1/2} \hat{Rm}_S^{1/4}]$']);
    % xlabel('$ \sqrt{\mathcal{C} Rm^*} St $');
    %xlabel('$\psi_a St h / L$');
    ylabel('$L$');
    %title('(c)');
    legend(CR_leg, 'Location', 'eastoutside');
    
    set(gca, 'Position', axPos(3,:));
    
    
    
    if savePlots
        set(hScalingArgsTalk,'Units','Inches');
        pos = get(hScalingArgsTalk,'Position');
        set(hScalingArgsTalk,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(gcf, '-dpdf', [figure_output_dir, 'scalingResultTalk.pdf'], printQuality);
        print(gcf, '-depsc', [figure_output_dir, 'scalingResultTalk.eps'], printQuality);
        
        
        
        print(gcf, '-dsvg', [figure_output_dir, 'scalingResultTalk.svg'], printQuality);
        fprintf('Saved file to %s \n', [figure_output_dir, 'scalingResultTalk.svg']);
    end
end




% Plot L(C)
plotL = false;

if plotL
    hLPlot = figure();
    set(hLPlot, 'Position', [300 300 450 300]);
    
    hold on;
    %CR_leg = {};
    Rm_leg = {};
    for Rm_i = 1:length(Rm_arr)
        
        flux = optimalFluxMat(:, Rm_i);
        L = optimalWidthMat(:, Rm_i);
        h = optimalChannelHMat(:, Rm_i);
        a = optimalaMat(:, Rm_i);
        psi = optimalPsiMat(:, Rm_i);
        Pi = averagePermMat(:, Rm_i); %
        H = optimalHMat(:, Rm_i);
        
        actualFluxes= ~isnan(flux).*(flux>1e-5).*(psi>1e-2);
        
        
        % Need to use this form
        RmStar = (Rm_arr-Rm_crit(CR_i))/Rm_crit(CR_i);
        RmStar = (Rm_arr)/Rm_crit(CR_i);
        
        RmPlot = [];
        
        Rm = Rm_arr(Rm_i);
        if Rm < 200 || Rm > 600 || Rm == 350 || Rm == 450 || Rm == 230 || Rm == 240
            continue
        end
        
        if Rm == 300
            actualFluxes = actualFluxes.*(CR_arr.'<0.5);
        end
        % if (CR > 0.5) || CR <0.02 || isnan(CR)
        %      continue
        %  end
        
        
        
        x  = CR_arr;
        
        y = L;
        
        if sum(actualFluxes) > 0
            %plot(log10(x(actualFluxes==1)) , log10(y(actualFluxes==1)));
            plot((x(actualFluxes==1)) , (y(actualFluxes==1)));
            % CR_leg{end+1} = ['$\mathcal{C}=', num2str(CR), '$'];
            Rm_leg{end+1} = ['$Rm_S=',num2str(Rm_arr(Rm_i)),'$'];
        end
    end
    hold off;
    box on;
    
    xlim([0 0.5]);
    
    xlabel('$\mathcal{C}$');
    
    ylabel('$L$');
    
    legend(Rm_leg, 'Location', 'eastoutside');
    
    % if makePlots
    %     set(hLPlot,'Units','Inches');
    %     pos = get(hLPlot,'Position');
    %     set(hLPlot,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    %     print(gcf, '-dpdf', [figure_output_dir, 'OptimalWidthLowC.pdf'], printQuality);
    % end
    
end



doStabilityPlots = false;

if doStabilityPlots
    
    
    
    
    
    
    
    hStability = figure();
    set(hStability, 'Position', [100 100 800 600]);
    m = 2; n = 2;
    
    
    h_valid =(h_crit > 0).*(h_crit < 0.4);
    
    
    subplot(m, n, 1);
    hold on
    plot(log10(CR_arr(h_valid==1)), log10(h_crit(h_valid==1)));
    third = plot([-2.4 0], [-2.4 -1.6], '--');
    half = plot([-2.5 -0.1], [-2.3 -1.1], '--');
    hold off;
    box on;
    
    legend([third half], {'$\mathcal{C}^{1/3}$', '$\mathcal{C}^{1/2}$'}, 'Location', 'southeast');
    
    hold off;
    xlabel('log$_{10}(\mathcal{C})$');
    ylabel('log$_{10}(\hat{h}_c)$');
    
    title('(a)');
    
    
    subplot(m, n, 2);
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
    title('(b)');
    
    
    subplot(m, n, 3);
    hold on
    %plot(log10(CR_arr), log10(1./Rm_crit));
    %twoThirds = plot([-2 1], [-2.5 -0.5], '--');
    %one = plot([-2 0], [-2.75 -0.75], '--');
    %two = plot([-2 -1], [-3 -1], '--');
    
    plot(log10(CR_arr), log10(Rm_crit));
    twoThirds = plot([-2 1], [2.5 0.5], '--');
    %one = plot([-2 0], [2.75 0.75], '--');
    %two = plot([-2 -1], [3 1], '--');
    
    hold off
    box on;
    xlabel('log$_{10}(\mathcal{C})$');
    ylabel('log$_{10}(\hat{Rm})$');
    %legend([twoThirds, one, two], {'$\mathcal{C}^{2/3}$', '$\mathcal{C}$', '$\mathcal{C}^2$'}, 'Location', 'northwest');
    legend([twoThirds,], ...
        {'$\mathcal{C}^{-2/3}$'}, ...
        'Location', 'southwest');
    title('(c)');
    
    
    
    subplot(m, n, 4);
    
    
    min_i = find(CR_arr==0.01);
    
    
    hold on;
    y = Rm_crit.*Pi_crit;
    Pi = plot(log10(CR_arr), y/nanmean(y(1:5)));
    %Pi = plot(log10(CR_arr), y/y(min_i));
    
    
    y = Rm_crit.*h_crit.^2;
    hSquared = plot(log10(CR_arr(h_valid==1)), y(h_valid==1)/nanmean(y(1:5)));
    %hSquared = plot(log10(CR_arr(h_valid==1)), y(h_valid==1)/y(min_i));
    y = Rm_crit.*h_crit.*Pi_crit;
    hPi= plot(log10(CR_arr(h_valid==1)), y(h_valid==1)/nanmean(y(1:5)));
    %hPi= plot(log10(CR_arr(h_valid==1)), y(h_valid==1)/y(min_i));
    
    y = Rm_crit.*h_crit;
    h= plot(log10(CR_arr(h_valid==1)), y(h_valid==1)/nanmean(y(1:5)));
    
    hold off;
    
    
    ylim([0.3 6]);
    xlim([-2 0]);
    ylabel('$y/\bar{y}(\mathcal{C} < 0.1)$');
    
    legend([Pi, hSquared, hPi, h],...
        {'$y=\hat{Rm}_S \hat{\Pi}$', '$y=\hat{Rm}_S  \hat{h}_c^2$',...
        '$y=\hat{Rm}_S  \hat{h}_c \hat{\Pi}$',  '$y=\hat{Rm}_S  \hat{h}_c$'}, ...
        'Location', 'northwest');
    
    xlabel('log$_{10}(\mathcal{C})$');
    title('(d)');
    box on;
    
    if savePlots
        set(hStability,'Units','Inches');
        pos = get(hStability,'Position');
        set(hStability,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(gcf, '-dpdf', [figure_output_dir, 'convectiveOnset.pdf'], printQuality);
    end
    
    
    
    
    
    
end



% Poster plot
doPosterScalingPlot = false;
if doPosterScalingPlot
    
    hPosterScaling = figure();
    set(hPosterScaling, 'Position', [200 200 450 450]);
    hold on;
    CR_leg = {};
    for CR_i = 1:length(CR_arr)
        
        flux = optimalFluxMat(CR_i,:);
        CR = CR_arr(CR_i);
        L = optimallMat(CR_i,:);
        h = optimalhMat(CR_i,:);
        psi = optimalPsiMat(CR_i,:);
        
        actualFluxes= ~isnan(flux);
        
        % Need to use this form
        RmStar = (Rm_arr-Rm_crit(CR_i))/Rm_crit(CR_i);
        
        %L = L(actualFluxes);
        %RmStar = Rm_arr/Rm_crit(CR_i);
        if (CR > 0.5) || CR <0.02 || isnan(CR)
            continue
        end
        
        x = (Rm_arr-Rm_crit(CR_i)).^(3/4) / (St^(1/2)*sqrt(Le));
        
        y = flux/(CR.^(1.0));
        
        if sum(actualFluxes) > 0
            plot1 = plot(x(actualFluxes) , y(actualFluxes));
            CR_leg{end+1} = ['$\mathcal{C}=', num2str(CR), '$'];
        end
    end
    hold off;
    box on;
    
    % xlim([0 0.1]);
    ylim([0 0.5]);
    
    
    xlabel(['$   (Rm_S-\hat{Rm}_S)^{3/4} [V_f / (St \, Le ) ]^{1/2}$']);
    % xlabel('$ \sqrt{\mathcal{C} Rm^*} St $');
    %xlabel('$\psi_a St h / L$');
    ylabel('$F_O/\mathcal{C}$');
    
    %legend(CR_leg, 'Location', 'eastoutside');
    legend(CR_leg, 'Location', 'northoutside', 'Orientation','horizontal');
    
    set(gca, 'Position', [0.21 0.2 0.65 0.5]);
    
    set(hPosterScaling,'Units','Inches');
    pos = get(hPosterScaling,'Position');
    set(hPosterScaling,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(gcf, '-depsc', [figure_output_dir, 'fluxScalingPoster.eps'], printQuality);
    print(gcf, '-dsvg', [figure_output_dir, 'fluxScalingPoster.svg'], printQuality);
    
end

end %end plotOptimalStates function


% All fields:
% H - mush height
% a - channel width
% psi_c - channel streamfunction
% h_c - channel height
% L - channel spacing
% psi - streamfunction in mush
% h_chi - height of high porous region
% chi - typical porosity
%
% will scale with
% S - stefan number
% V - frame advection
% R - Rayleigh number
% C - concentration ratio
% L - lewis number
%

% Scaling type can be:
% 0 - empirical scaling
% 1 - new justified scaling

function s = getExponents(scalingType)

baseScaling = struct('S', 0, 'V', 0, 'R', 0, 'C', 0, 'L', 0);

s.H = baseScaling;
s.a = baseScaling;
s.psi_c = baseScaling;
s.h_c = baseScaling;
s.L = baseScaling;
s.psi = baseScaling;
s.h_chi = baseScaling;
s.chi = baseScaling;
s.F = baseScaling;

%H, a, psi_c, h_c always scale the same
s.H.S = -1;
s.H.V = -1;

s.a.R = -3/4;
s.a.L = -1/2;
s.a.S = -1/2;
s.a.V = -1/2;

s.psi_c.R = 1/2;
s.psi_c.C = 1;
s.psi_c.L = -1/2;
s.psi_c.S = -1;

s.h_c.C = 1;
s.h_c.R = 1/4;
s.h_c.S = 3/2;
s.h_c.V = -1/2;

% L, psi, h_chi, chi have different scalings
if scalingType == 0
    
    s.L.R = -1/4;
    s.L.S = -1/2;
    s.L.V = -1/2;
    
    s.psi.C = 1;
    s.psi.R = 1/2;
    s.psi.S = -1;
    
    s.h_chi.C = 1 ;
    s.h_chi.R = -1/4 ;
    s.h_chi.S = -3/2;
    s.h_chi.V = -1/2;
    
    s.F.C = 1;
    s.F.R = 3/4;
    s.F.S = -1/2;
    s.F.L = -1/2;
    s.F.V = 1/2;
    
    %s.chi. = ;
elseif scalingType == 1
    
    s.L.R = -3/40;
    s.L.C = 1/5;
    s.L.V = -11/20;
    
    s.psi.C = 4/5;
    s.psi.R = 13/40;
    s.psi.S = -3/4;
    s.psi.V = 1/20;
    
    s.h_chi.C = 3/5;
    s.h_chi.R = -3/5;
    s.h_chi.S = -1;
    s.h_chi.V = -2/5;
    
    s.chi.C = 1/5;
    s.chi.V = 1/5;
    s.chi.R = -1/5;
    
    s.F.C = 4/5; %4/5;
    s.F.R = 23/40;
    s.F.S = -1;
    s.F.L = -1/2;
    s.F.V = 11/20;
    
elseif scalingType == 2
    
    s.L.R = -1/4;
    s.L.S = -2/3;
    s.L.V = -1/2;
    
    s.psi.C = 1;
    s.psi.R = 1/2;
    s.psi.S = -5/6;
    
    s.h_chi.C = 1 ;
    s.h_chi.R = -1/4 ;
    s.h_chi.S = -7/6;
    s.h_chi.V = -1/2;
    
    s.F.C = 1;
    s.F.R = 3/4;
    s.F.S = -1/3;
    s.F.L = -1/2;
    s.F.V = 1/2;
   
        
end




end



function  [plts, CR_leg] = plotScaling(ylab, scaling, optimalMetric, Rm_arr, Rm_crit, CR_arr, St, Le, Vf)

% Should do something different for ylab=='$F_O/\mathcal{C}$'

hold on;
CR_leg = {};
plts = [];

xlab = '';

scalingString = '';
if length(Rm_crit) > 0
    scalingString = '(Rm_S-\hat{Rm}_S)';
else
    scalingString = 'Rm_S';
end

xlab = printScaling(xlab, scalingString, scaling.R);


%Special cases
if sum(strfind(ylab, 'F_O')) == 0
    xlab = printScaling(xlab, '\mathcal{C}', scaling.C);
end


xlab = printScaling(xlab, 'V_f', scaling.V);
xlab = printScaling(xlab, '\mathcal{S}', scaling.S);
xlab = printScaling(xlab, 'Le', scaling.L);

for CR_i = 1:length(CR_arr)
    
    metrics = optimalMetric(CR_i,:);
    CR = CR_arr(CR_i);
    
    validMetrics = ~isnan(metrics);
    
    % Need to use this form
    if length(Rm_crit) > 0
        % Rm_plot = (Rm_arr-Rm_crit(CR_i))/Rm_crit(CR_i);
        Rm_plot = Rm_arr-Rm_crit(CR_i);
    else
        Rm_plot = Rm_arr;
    end
    
    if (CR > 0.5) || CR <0.02 || isnan(CR)
        continue
    end
    
    x = Rm_plot.^scaling.R ...
        * St^scaling.S ...
        * Le^scaling.L ...
        * Vf^scaling.V;
    
    y = metrics;
    
    if sum(strfind(ylab, 'F_O')) > 0
        y = y/(CR^scaling.C);
    else
        x = x*CR^scaling.C;
    end
    
    
    
    if sum(validMetrics) > 0
        plts(end+1) =  plot(x(validMetrics) , y(validMetrics));
        CR_leg{end+1} = ['$\mathcal{C}=', num2str(CR), '$'];
    end
end
hold off;
box on;


xlabel(['$',xlab,'$']);
ylabel(ylab);
ax = gca;
ax.XRuler.Exponent = 0; % Turn off scientific notation

end

function xlab = printScaling(xlab, text, power)
if power ~= 0
    
    xlab = [xlab, text];
    
    if power == 1
        xlab = [xlab, '\,'];
    else
        xlab = [xlab, '^{',  strtrim(rats(power)), '}'];
    end
    
end



end
