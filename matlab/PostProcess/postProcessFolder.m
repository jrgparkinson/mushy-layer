function postProcessFolder(output_dir, makePlots, redoRuns, katzUnits, requireSteady, simpleAnalysis)
close all;

%data_dir = ['/home/parkinsonjl/mnt/sharedStorage/', ...
%     'optimalStates-highRes-new/'];

if nargin < 6
    simpleAnalysis = true;
end

if nargin < 5
    requireSteady = false; % only process steady folders
end

if nargin < 4
    katzUnits = true;
end

if nargin < 3
    redoRuns = true;
end

if nargin < 2
    makePlots = true;
end

csvFileName = 'postProcess.csv';
if exist([output_dir, csvFileName] ,'file') == 2 && ~redoRuns
    return
end

%set(groot, 'DefaultTextInterpreter', 'latex')
%set(groot, 'DefaultLegendInterpreter', 'latex')

% Options:
scalePlots = false;

ending = '';
if requireSteady
    ending = 'steady';
end

fprintf('====================================\n');
fprintf('Processing folder: %s \n', output_dir);

outPref = '  -';

%folders_regex = ['CR(\d+\.?\d+?)RaC(\d+\.?\d+?)Le(\d+e?\.?\d+?)(\w+)pts(\d+)-',ending]; %Make sure it ends in steady

%full_output_dir = [data_dir, folderName, '/'];
%output_dir = [data_dir, folderName, '/'];


% New method - get these params from the inputs file

try
    inputs = readInputs([output_dir, 'inputs']);
    St = str2double(inputs.stefan);
    frameAdv = str2double(inputs.nonDimVel);
    
    CR = str2double(inputs.compositionRatio);
    Ra = str2double(inputs.rayleighComp);
    Le = str2double(inputs.lewis);
    %PermType = str2num(inputs.permeabilityFunction);
    num_cells = inputs.num_cells;
    p = strsplit(num_cells, ' ');
   % width = str2num(p{1});
   
catch e %e is an MException struct
    fprintf('%s Could not get parameters from inputs file, skipping \n', outPref);
    fprintf(1,'%s The identifier was:\n%s \n',outPref, e.identifier);
    fprintf(1,'%s There was an error! The message was:\n%s \n',outPref, e.message);
    
    return
end


% Find and load the final plot file, and calculate various
% diagnostics

diagnostics = struct();





% Cover this in a try/catch so we keep processing other folders if
% there's an issue
try
    finalPlotFile = getFinalPlotFile(output_dir);
catch e
    fprintf('%s Could not load plot file, skipping \n', outPref);
    fprintf(1,'%s The identifier was:\n%s \n',outPref, e.identifier);
    fprintf(1,'%s There was an error! The message was:\n%s \n',outPref, e.message);
    return
end

% If we didn't open a valid file, skip this
if length(isprop(finalPlotFile, 'levelArray') ) == 0
    return
elseif length(finalPlotFile.levelArray) == 0
    return
end

% First let's calculate the size of each term in the equations
% a) throughout the entire mushy layer
% b) throughout the interior of the mushy layer

T = finalPlotFile.dataForComp(finalPlotFile.components.Temperature).';
chi = finalPlotFile.dataForComp(finalPlotFile.components.Porosity).';
%     H = finalPlotFile.dataForComp(finalPlotFile.components.Enthalpy).';
S = finalPlotFile.dataForComp(finalPlotFile.components.Bulkconcentration).';
Sl = finalPlotFile.dataForComp(finalPlotFile.components.Liquidconcentration).';
%
%     U = finalPlotFile.dataForComp(finalPlotFile.components.xAdvectionvelocity).';
V = finalPlotFile.dataForComp(finalPlotFile.components.yAdvectionvelocity).';
psi = finalPlotFile.getStreamfunction(1e4);


%     %psi2 =  finalPlotFile.getStreamfunction(1e5);
%     %psiDiff = (psi2-psi)./psi;
%     % psiMaxFracDiff = max(max(psiDiff));
%      chiSl = chi.*Sl;
%
perm = chi.^3;

% Do some conversion to Wells Units
CR = CR - 1; % Always do this, as CR comes from the folder name

if katzUnits
    
    % H = H - 1;
    T = T - 1;
    S = S + 1;
    Sl = Sl + 1;
end






% Also get the flux and width
if katzUnits
    
    % See if we've done post processing via cpp code,
    % which should give us the best salt flux calculation
    
    cpp_diagFileLoc = [output_dir, 'cpp_diagnostics.out'];
    %fprintf('Diag file: %s \n', diagFileLoc);
    
    gotDiags = false;
    
    
    
    if exist(cpp_diagFileLoc, 'file') == 2
    
        try

            cpp_diags = getDiagnostics(cpp_diagFileLoc);

            % Take the negative of this.
            diagnostics.flux = -cpp_diags.Fs_average(end);


            diagnostics.L2FsVertDiffusion = cpp_diags.L2FsVertDiffusion(end);
            diagnostics.L2FsVertFluid = cpp_diags.L2FsVertFluid(end);
            diagnostics.L2FsVertFrame = cpp_diags.L2FsVertFrame(end);

            diagnostics.L1FsVertDiffusion = cpp_diags.L1FsVertDiffusion(end);
            diagnostics.L1FsVertFluid = cpp_diags.L1FsVertFluid(end);
            diagnostics.L1FsVertFrame = cpp_diags.L1FsVertFrame(end);

            diagnostics.L0FsVertDiffusion = cpp_diags.L0FsVertDiffusion(end);
            diagnostics.L0FsVertFluid = cpp_diags.L0FsVertFluid(end);
            diagnostics.L0FsVertFrame = cpp_diags.L0FsVertFrame(end);
            
            gotDiags = true;

        catch e %e is an MException struct
            fprintf('%s Could not get cpp diagnostic file, skipping \n', outPref);
            fprintf(1,'%s The identifier was:\n%s \n',outPref, e.identifier);
            fprintf(1,'%s There was an error! The message was:\n%s \n',outPref, e.message);
            gotDiags = false;


        end
    
    else
        
        % Get diags for latest diagnostics.out
      try
          
          diagFileLoc = [output_dir, 'diagnosticsLatest.out'];
%fprintf('Diag file: %s \n', diagFileLoc);

try
    diagFile = getDiagnostics(diagFileLoc);
catch e %e is an MException struct
    fprintf('%s Could not get diagnostic file, skipping \n', outPref);
    fprintf(1,'%s The identifier was:\n%s \n',outPref, e.identifier);
    fprintf(1,'%s There was an error! The message was:\n%s \n',outPref, e.message);
    
    return
end

           

            % Take the negative of this.
            diagnostics.flux = -diagFile.Fs_vertical_av(end);


            diagnostics.L2FsVertDiffusion = diagFile.L2FsVertDiffusion(end);
            diagnostics.L2FsVertFluid = diagFile.L2FsVertFluid(end);
            diagnostics.L2FsVertFrame = diagFile.L2FsVertFrame(end);

            diagnostics.L1FsVertDiffusion = diagFile.L1FsVertDiffusion(end);
            diagnostics.L1FsVertFluid = diagFile.L1FsVertFluid(end);
            diagnostics.L1FsVertFrame = diagFile.L1FsVertFrame(end);

            diagnostics.L0FsVertDiffusion = diagFile.L0FsVertDiffusion(end);
            diagnostics.L0FsVertFluid = diagFile.L0FsVertFluid(end);
            diagnostics.L0FsVertFrame = diagFile.L0FsVertFrame(end);
            
            gotDiags = true;

        catch e %e is an MException struct
            fprintf('%s Could not get diagnostics from diag file file, skipping \n', outPref);
            fprintf(1,'%s The identifier was:\n%s \n',outPref, e.identifier);
            fprintf(1,'%s There was an error! The message was:\n%s \n',outPref, e.message);
            gotDiags = false;


        end
        
    end
    
    if ~gotDiags
        diagnostics.flux = computeSaltFlux(S, V, Sl, chi, dx, frameAdv, Le);
    end
    
    fprintf('%s Salt flux = %1.7f \n', outPref, diagnostics.flux);
    
else
    
    if isfield(diagFile, 'Fs_vertical_av')
        diagnostics.flux = diagFile.Fs_vertical_av(end);
    elseif isfield(diagFile, 'Fs_bottom')
        diagnostics.flux = diagFile.Fs_bottom(end);
    end
    
    % Convert flux (unless our code has already done it)
    if diagnostics.flux < 0
        
        if katzUnits
            % Katz Units
            % Shouldn't end up here anymore, and this conversion isn't
            % correct
            diagnostics.flux = -1 - diagnostics.flux;
        else
            % Wells units
            diagnostics.flux = -diagnostics.flux;
        end
    end
    
end


% Get all the fields and do more complicated analysis

if ~simpleAnalysis
    
try
    
[heatAdvection, heatDiffusion, latentHeat, TFrameAdvection, ...
    saltAdvection, saltDiffusion, liquidSalinityFrame, solidSalinityFrame, ...
    vorticityDiffusion, baroclinicTorque, vorticityPermeability] = ...
    computeFields(finalPlotFile, frameAdv, St, CR, Ra, Le);

[X, Y] = finalPlotFile.grid();
dx = X(1, 2) - X(1, 1);

% Get estimate for eutectic porosity
chiEutectic = CR/(CR + 1);

% Now just take the mushy region
minPorosity = chiEutectic*1.2; % Make sure we avoid the eutectic
minPorosity = 0.005;
maxPorosity = 0.995;

% Also want to avoid the sharp temp gradients at the eutectic
% Get approx eutectic position and remove cells near this

mushyLayer = (chi > minPorosity).*(chi < maxPorosity).*(T > -0.95);

porosityMush = chi; porosityMush(mushyLayer ~= 1) = NaN;

streamfunctionMush = psi; streamfunctionMush(mushyLayer ~= 1) = NaN;
permeabilityMush = perm; permeabilityMush(mushyLayer ~= 1) = NaN;

heatAdvectionMush = heatAdvection; heatAdvectionMush(mushyLayer ~= 1) = NaN;
heatDiffusionMush = heatDiffusion; heatDiffusionMush(mushyLayer ~= 1) = NaN;
latentHeatMush = latentHeat; latentHeatMush(mushyLayer ~= 1) = NaN;
TFrameAdvectionMush = TFrameAdvection; TFrameAdvectionMush(mushyLayer ~= 1) = NaN;

saltAdvectionMush = saltAdvection; saltAdvectionMush(mushyLayer ~= 1) = NaN;
saltDiffusionMush = saltDiffusion; saltDiffusionMush(mushyLayer ~= 1) = NaN;
liquidSalinityFrameMush = liquidSalinityFrame; liquidSalinityFrameMush(mushyLayer ~= 1) = NaN;
solidSalinityFrameMush = solidSalinityFrame; solidSalinityFrameMush(mushyLayer ~= 1) = NaN;

vorticityDiffusionMush = vorticityDiffusion; vorticityDiffusionMush(mushyLayer ~= 1) = NaN;
baroclinicTorqueMush = baroclinicTorque; baroclinicTorqueMush(mushyLayer ~= 1) = NaN;
vorticityPermeabilityMush = vorticityPermeability; vorticityPermeabilityMush(mushyLayer ~= 1) = NaN;



%Calculate some geometric properties
%1) Channel height
[diagnostics.channelWidth, diagnostics.channelHeight, ~, ~] = finalPlotFile.channelGeometry();

% Width of active region around chimney
% Determined by distance from psi maximum to the edge
%Find position of psi maximum
maxPsi = max(streamfunctionMush(:));

[y_i, x_i] = find(streamfunctionMush==maxPsi);

maxXi = length(X(1, :));

if length(x_i) ~= 1
    % If no unique maximum, return NaN (happens when, e.g. there is
    % no flow)
    diagnostics.largeLScale = NaN;
    diagnostics.smallLScale = NaN;
else
    
    distanceFromLeft = (x_i)*dx;
    distanceFromRight = (maxXi-x_i)*dx;
    
    diagnostics.largeLScale = max(distanceFromLeft, distanceFromRight);
    diagnostics.smallLScale = min(distanceFromLeft, distanceFromRight);
end


% Find penetration depth
psiHorizMax = max(streamfunctionMush, [], 2);
YPsi_i = find(psiHorizMax < maxPsi*0.25);

Y_psi_penetration = min(YPsi_i);

chiHorizMin = min(porosityMush, [], 2);

Y_chi_confinement = min(find(chiHorizMin < 0.75));



% Height of mush-liquid boundary layer
averageMushyChi = nanmean(nanmean(porosityMush));

Y_i_avPoros = [];
Y_i_mushLiquid = [];
Y_i_solid = [];

x = X(1, :); y = Y(:, 1);

for x_i = 1:length(x)
    porosityColumn = squeeze(chi(:, x_i));
    
    
    y_avPorosity_i = find(porosityColumn < averageMushyChi);
    if length(y_avPorosity_i) > 0
        Y_i_avPoros(end+1) = min(y_avPorosity_i); % Should this be max?
    end
    
    y_maxPorosity_i = find(porosityColumn < maxPorosity);
    if length(y_maxPorosity_i) > 0
        Y_i_mushLiquid(end+1) = min(y_maxPorosity_i);
    end
    
    
    y_solid_i = find(porosityColumn < minPorosity);
    if length(y_solid_i) > 0
        Y_i_solid(end+1) = min(y_solid_i);
    end
    
end




median_avPorosity_Yi = median(Y_i_avPoros);
median_mushyLiquid_Yi = median(Y_i_mushLiquid);
median_solid_Yi = median(Y_i_solid);

min_solid_Yi = min(Y_i_solid);

diagnostics.h = abs(median_avPorosity_Yi-median_mushyLiquid_Yi)*dx;
diagnostics.H = abs(median_solid_Yi-median_mushyLiquid_Yi)*dx;

diagnostics.hPsi = abs(Y_psi_penetration - median_mushyLiquid_Yi)*dx;
diagnostics.hChi = abs(Y_chi_confinement - median_mushyLiquid_Yi)*dx;

%[Y_i, X_i] = find(~isnan(porosityMush));

% Also find the average porosity at position h and H
diagnostics.hporosity = nanmean(squeeze(chi(round(median_avPorosity_Yi),:)));
diagnostics.Hporosity = nanmean(squeeze(chi(min_solid_Yi-2, :)));



if makePlots
    
    h = figure(1);
    set(h, 'Position', [50, 50, 1600, 1000]);
    
    m = 3; n = 5;
    
    
    maxT = max([heatAdvectionMush(:); heatDiffusionMush(:); latentHeatMush(:); TFrameAdvectionMush(:)]);
    minT = min([heatAdvectionMush(:); heatDiffusionMush(:); latentHeatMush(:); TFrameAdvectionMush(:)]);
    
    largeT = max(abs(maxT), abs(minT));
    
    subplot(m, n, 1);
    h = pcolor(heatAdvectionMush);
    set(h, 'EdgeColor', 'none');
    if scalePlots
        ax = gca; ax.CLim = [-largeT largeT];
    end
    colormap(gca, bluewhitered);
    colorbar(); set(gca,'YDir','normal')
    title('$U \cdot \nabla T$');
    
    
    subplot(m, n, 2);
    h = pcolor(heatDiffusionMush);
    set(h, 'EdgeColor', 'none');
    if scalePlots
        ax = gca; ax.CLim = [-largeT largeT];
    end
    colormap(gca, bluewhitered);
    colorbar(); set(gca,'YDir','normal')
    title('$\nabla^2 T$');
    
    
    subplot(m, n, 3);
    h = pcolor(latentHeatMush);
    set(h, 'EdgeColor', 'none');
    if scalePlots
        ax = gca; ax.CLim = [-largeT largeT];
    end
    colormap(gca, bluewhitered);
    colorbar(); set(gca,'YDir','normal')
    title('$V St d \chi / dz$');
    
    
    subplot(m, n, 4);
    h = pcolor(TFrameAdvectionMush);
    set(h, 'EdgeColor', 'none');
    if scalePlots
        ax = gca; ax.CLim = [-largeT largeT];
    end
    colormap(gca, bluewhitered);
    colorbar(); set(gca,'YDir','normal')
    title('$V dT/dz$');
    
    
    subplot(m, n, 5);
    h = pcolor(heatAdvectionMush + latentHeatMush + TFrameAdvectionMush - heatDiffusionMush);
    set(h, 'EdgeColor', 'none');
    if scalePlots
        ax = gca; ax.CLim = [-largeT largeT];
    end
    colormap(gca, bluewhitered);
    colorbar(); set(gca,'YDir','normal')
    title('Heat residual');
    
    
    % Salt equation plots
    maxS = max([saltAdvectionMush(:); saltDiffusionMush(:); liquidSalinityFrameMush(:); solidSalinityFrameMush(:)]);
    minS = min([saltAdvectionMush(:); saltDiffusionMush(:); liquidSalinityFrameMush(:); solidSalinityFrameMush(:)]);
    
    largeS = max(abs(maxS), abs(minS));
    
    
    subplot(m, n, 6);
    h = pcolor(saltAdvectionMush);
    set(h, 'EdgeColor', 'none');
    if scalePlots
        ax = gca; ax.CLim = [-largeS largeS];
    end
    colormap(gca, bluewhitered);
    colorbar(); set(gca,'YDir','normal')
    title('$U \cdot \nabla S_l$');
    
    subplot(m, n, 7);
    h = pcolor(saltDiffusionMush);
    set(h, 'EdgeColor', 'none');
    if scalePlots
        ax = gca; ax.CLim = [-largeS largeS];
    end
    colormap(gca, bluewhitered);
    colorbar(); set(gca,'YDir','normal')
    title('$\nabla \cdot \chi \nabla S_l $');
    
    subplot(m, n, 8);
    h = pcolor(liquidSalinityFrameMush);
    set(h, 'EdgeColor', 'none');
    if scalePlots
        ax = gca; ax.CLim = [-largeS largeS];
    end
    colormap(gca, bluewhitered);
    colorbar(); set(gca,'YDir','normal')
    title('$V  d (\chi S_l) / dz$');
    
    subplot(m, n, 9);
    h = pcolor(solidSalinityFrameMush);
    set(h, 'EdgeColor', 'none');
    if scalePlots
        ax = gca; ax.CLim = [-largeS largeS];
    end
    colormap(gca, bluewhitered);
    colorbar(); set(gca,'YDir','normal')
    title('$V C d (1-\chi) /dz$');
    
    subplot(m, n, 10);
    h = pcolor(saltAdvectionMush + liquidSalinityFrameMush - solidSalinityFrameMush - saltDiffusionMush);
    set(h, 'EdgeColor', 'none');
    if scalePlots
        ax = gca; ax.CLim = [-largeS largeS];
    end
    colormap(gca, bluewhitered);
    colorbar(); set(gca,'YDir','normal')
    title('Salt residual');
    
    
    % Momentum equation plots
    maxMom = max([vorticityDiffusionMush(:); baroclinicTorqueMush(:); vorticityPermeabilityMush(:)]);
    minMom = min([vorticityDiffusionMush(:); baroclinicTorqueMush(:); vorticityPermeabilityMush(:)]);
    
    largeMom = max(abs(maxMom), abs(minMom));
    
    
    subplot(m, n, 11);
    h = pcolor(vorticityDiffusionMush);
    set(h, 'EdgeColor', 'none');
    if scalePlots
        ax = gca; ax.CLim = [-largeMom largeMom];
    end
    colormap(gca, bluewhitered);
    colorbar(); set(gca,'YDir','normal')
    title('$\nabla^2 \psi$');
    
    subplot(m, n, 12);
    h = pcolor(baroclinicTorqueMush);
    set(h, 'EdgeColor', 'none');
    if scalePlots
        ax = gca; ax.CLim = [-largeMom largeMom];
    end
    colormap(gca, bluewhitered);
    colorbar(); set(gca,'YDir','normal')
    title('$Rm \Pi dT/dx$');
    
    subplot(m, n, 13);
    h = pcolor(vorticityPermeabilityMush);
    set(h, 'EdgeColor', 'none');
    if scalePlots
        ax = gca; ax.CLim = [-largeMom largeMom];
    end
    colormap(gca, bluewhitered);
    colorbar(); set(gca,'YDir','normal')
    title('$(\nabla \psi \cdot \nabla \Pi)/\Pi$');
    
    subplot(m, n, 14);
    h = pcolor(vorticityDiffusionMush + baroclinicTorqueMush - vorticityPermeabilityMush);
    caxis([-3 3]);
    set(h, 'EdgeColor', 'none');
    if scalePlots
        ax = gca; ax.CLim = [-largeMom largeMom];
    end
    colormap(gca, bluewhitered);
    colorbar(); set(gca,'YDir','normal')
    title('Momentum residual');
    
    
    
    %         subplot(m, n, 15);
    %         h = pcolor(heatAdvectionMush + latentHeatMush + TFrameAdvectionMush - heatDiffusionMush);
    %         set(h, 'EdgeColor', 'none');
    %         colormap(gca, bluewhitered);
    %         colorbar(); set(gca,'YDir','normal')
    %         title('Heat equation residual');
    
    
end % end if make plots


diagnostics.heatAdvectionMush = max(max(abs(heatAdvectionMush)));
diagnostics.latentHeatMush = max(max(abs(latentHeatMush)));
diagnostics.TFrameAdvectionMush = max(max(abs(TFrameAdvectionMush)));
diagnostics.heatDiffusionMush = max(max(abs(heatDiffusionMush)));


diagnostics.saltAdvectionMush = max(max(abs(saltAdvectionMush)));
diagnostics.liquidSalinityFrameMush = max(max(abs(liquidSalinityFrameMush)));
diagnostics.solidSalinityFrameMush = max(max(abs(solidSalinityFrameMush)));
diagnostics.saltDiffusionMush = max(max(abs(saltDiffusionMush)));

diagnostics.baroclinicTorqueMush = max(max(abs(baroclinicTorqueMush)));
diagnostics.vorticityPermeabilityMush = max(max(abs(vorticityPermeabilityMush)));
diagnostics.vorticityDiffusionMush = max(max(abs(vorticityDiffusionMush)));


diagnostics.maxPorosityMush = max(max(abs(porosityMush)));
diagnostics.maxStreamfunctionMush = max(max(abs(streamfunctionMush)));
diagnostics.maxPermeabilityMush = max(max(abs(permeabilityMush)));

diagnostics.avPorosityMush = nanmean(nanmean(porosityMush));
diagnostics.avStreamfunctionMush = nanmean(nanmean(streamfunctionMush));
diagnostics.avPermeabilityMush = nanmean(nanmean(permeabilityMush));


% Now do the same but for interior points

averageMushPorosity = nanmean(nanmean(chi(mushyLayer==1)));

minPorosity = 0.05;
maxPorosity = 0.6;
maxPorosity = nanmedian(nanmedian(porosityMush));
%maxPorosity = (1+averageMushPorosity)*0.5;

% Also want to avoid the sharp temp gradients at the eutectic
% Get approx eutectic position and remove cells near this

InterioryLayer = (chi > minPorosity).*(chi < maxPorosity).*(T > -0.95);
heatAdvectionInterior = heatAdvection; heatAdvectionInterior(InterioryLayer ~= 1) = NaN;
heatDiffusionInterior = heatDiffusion; heatDiffusionInterior(InterioryLayer ~= 1) = NaN;
latentHeatInterior = latentHeat; latentHeatInterior(InterioryLayer ~= 1) = NaN;
TFrameAdvectionInterior = TFrameAdvection; TFrameAdvectionInterior(InterioryLayer ~= 1) = NaN;

saltAdvectionInterior = saltAdvection; saltAdvectionInterior(InterioryLayer ~= 1) = NaN;
saltDiffusionInterior = saltDiffusion; saltDiffusionInterior(InterioryLayer ~= 1) = NaN;
liquidSalinityFrameInterior = liquidSalinityFrame; liquidSalinityFrameInterior(InterioryLayer ~= 1) = NaN;
solidSalinityFrameInterior = solidSalinityFrame; solidSalinityFrameInterior(InterioryLayer ~= 1) = NaN;

vorticityDiffusionInterior = vorticityDiffusion; vorticityDiffusionInterior(InterioryLayer ~= 1) = NaN;
baroclinicTorqueInterior = baroclinicTorque; baroclinicTorqueInterior(InterioryLayer ~= 1) = NaN;
vorticityPermeabilityInterior = vorticityPermeability; vorticityPermeabilityInterior(InterioryLayer ~= 1) = NaN;


porosityInterior = chi;  porosityInterior(InterioryLayer ~= 1) = NaN;
streamfunctionInterior = psi; streamfunctionInterior(InterioryLayer ~= 1) = NaN;
permeabilityInterior = perm; permeabilityInterior(InterioryLayer ~= 1) = NaN;

if makePlots
    
    h = figure(2);
    set(h, 'Position', [50, 50, 1600, 1000]);
    
    m = 3; n = 5;
    
    
    maxT = max([heatAdvectionInterior(:); heatDiffusionInterior(:);
        latentHeatInterior(:); TFrameAdvectionInterior(:)]);
    minT = min([heatAdvectionInterior(:); heatDiffusionInterior(:);
        latentHeatInterior(:); TFrameAdvectionInterior(:)]);
    
    largeT = max(abs(maxT), abs(minT));
    
    
    subplot(m, n, 1);
    h = pcolor(heatAdvectionInterior);
    set(h, 'EdgeColor', 'none');
    if scalePlots
        ax = gca; ax.CLim = [-largeT largeT];
    end
    colormap(gca, bluewhitered);
    colorbar(); set(gca,'YDir','normal')
    title('$U \cdot \nabla T$');
    
    subplot(m, n, 2);
    h = pcolor(heatDiffusionInterior);
    set(h, 'EdgeColor', 'none');
    if scalePlots
        ax = gca; ax.CLim = [-largeT largeT];
    end
    colormap(gca, bluewhitered);
    colorbar(); set(gca,'YDir','normal')
    title('$\nabla^2 T$');
    
    subplot(m, n, 3);
    h = pcolor(latentHeatInterior);
    set(h, 'EdgeColor', 'none');
    if scalePlots
        ax = gca; ax.CLim = [-largeT largeT];
    end
    colormap(gca, bluewhitered);
    colorbar(); set(gca,'YDir','normal')
    title('$V St d \chi / dz$');
    
    subplot(m, n, 4);
    h = pcolor(TFrameAdvectionInterior);
    set(h, 'EdgeColor', 'none');
    if scalePlots
        ax = gca; ax.CLim = [-largeT largeT];
    end
    colormap(gca, bluewhitered);
    colorbar(); set(gca,'YDir','normal')
    title('$V dT/dz$');
    
    subplot(m, n, 5);
    h = pcolor(heatAdvectionInterior + latentHeatInterior + TFrameAdvectionInterior - heatDiffusionInterior);
    set(h, 'EdgeColor', 'none');
    if scalePlots
        ax = gca; ax.CLim = [-largeT largeT];
    end
    colormap(gca, bluewhitered);
    colorbar(); set(gca,'YDir','normal')
    title('Heat residual');
    
    
    % Salt equation plots
    maxS = max([saltAdvectionInterior(:); saltDiffusionInterior(:);
        liquidSalinityFrameInterior(:); solidSalinityFrameInterior(:)]);
    minS = min([saltAdvectionInterior(:); saltDiffusionInterior(:);
        liquidSalinityFrameInterior(:); solidSalinityFrameInterior(:)]);
    
    largeS = max(abs(maxS), abs(minS));
    
    subplot(m, n, 6);
    h = pcolor(saltAdvectionInterior);
    set(h, 'EdgeColor', 'none');
    if scalePlots
        ax = gca; ax.CLim = [-largeS largeS];
    end
    colormap(gca, bluewhitered);
    colorbar(); set(gca,'YDir','normal')
    title('$U \cdot \nabla S_l$');
    
    subplot(m, n, 7);
    h = pcolor(saltDiffusionInterior);
    set(h, 'EdgeColor', 'none');
    if scalePlots
        ax = gca; ax.CLim = [-largeS largeS];
    end
    colormap(gca, bluewhitered);
    colorbar(); set(gca,'YDir','normal')
    title('$\nabla \cdot \chi \nabla S_l $');
    
    subplot(m, n, 8);
    h = pcolor(liquidSalinityFrameInterior);
    set(h, 'EdgeColor', 'none');
    if scalePlots
        ax = gca; ax.CLim = [-largeS largeS];
    end
    colormap(gca, bluewhitered);
    colorbar(); set(gca,'YDir','normal')
    title('$V  d (\chi S_l) / dz$');
    
    subplot(m, n, 9);
    h = pcolor(solidSalinityFrameInterior);
    set(h, 'EdgeColor', 'none');
    if scalePlots
        ax = gca; ax.CLim = [-largeS largeS];
    end
    colormap(gca, bluewhitered);
    colorbar(); set(gca,'YDir','normal')
    title('$V C d (1-\chi) /dz$');
    
    subplot(m, n, 10);
    h = pcolor(saltAdvectionInterior + liquidSalinityFrameInterior - solidSalinityFrameInterior - saltDiffusionInterior);
    set(h, 'EdgeColor', 'none');
    if scalePlots
        ax = gca; ax.CLim = [-largeS largeS];
    end
    colormap(gca, bluewhitered);
    colorbar(); set(gca,'YDir','normal')
    title('Salt residual');
    
    
    % Momentum equation plots
    maxMom = max([vorticityDiffusionInterior(:); baroclinicTorqueInterior(:); vorticityPermeabilityInterior(:)]);
    minMom = min([vorticityDiffusionInterior(:); baroclinicTorqueInterior(:); vorticityPermeabilityInterior(:)]);
    
    largeMom = max(abs(maxMom), abs(minMom));
    
    subplot(m, n, 11);
    h = pcolor(vorticityDiffusionInterior);
    set(h, 'EdgeColor', 'none');
    if scalePlots
        ax = gca; ax.CLim = [-largeMom largeMom];
    end
    colormap(gca, bluewhitered);
    colorbar(); set(gca,'YDir','normal')
    title('$\nabla^2 \psi$');
    
    subplot(m, n, 12);
    h = pcolor(baroclinicTorqueInterior);
    set(h, 'EdgeColor', 'none');
    if scalePlots
        ax = gca; ax.CLim = [-largeMom largeMom];
    end
    colormap(gca, bluewhitered);
    colorbar(); set(gca,'YDir','normal')
    title('$Rm \Pi dT/dx$');
    
    subplot(m, n, 13);
    h = pcolor(vorticityPermeabilityInterior);
    set(h, 'EdgeColor', 'none');
    if scalePlots
        ax = gca; ax.CLim = [-largeMom largeMom];
    end
    colormap(gca, bluewhitered);
    colorbar(); set(gca,'YDir','normal')
    title('$(\nabla \psi \cdot \nabla \Pi)/\Pi$');
    
    subplot(m, n, 14);
    h = pcolor(vorticityDiffusionInterior + baroclinicTorqueInterior - vorticityPermeabilityInterior);
    caxis([-3 3]);
    set(h, 'EdgeColor', 'none');
    if scalePlots
        ax = gca; ax.CLim = [-largeMom largeMom];
    end
    colormap(gca, bluewhitered);
    colorbar(); set(gca,'YDir','normal')
    title('Momentum residual');
    box on;
    
    
    
    %         subplot(m, n, 15);
    %         h = pcolor(heatAdvectionInterior + latentHeatInterior + TFrameAdvectionInterior - heatDiffusionInterior);
    %         set(h, 'EdgeColor', 'none');
    %         colormap(gca, bluewhitered);
    %         colorbar(); set(gca,'YDir','normal')
    %         title('Heat equation residual');
    
    
    
end % end if make plots


diagnostics.heatAdvectionInterior = max(max(abs(heatAdvectionInterior)));
diagnostics.latentHeatInterior = max(max(abs(latentHeatInterior)));
diagnostics.TFrameAdvectionInterior = max(max(abs(TFrameAdvectionInterior)));
diagnostics.heatDiffusionInterior = max(max(abs(heatDiffusionInterior)));

diagnostics.saltAdvectionInterior = max(max(abs(saltAdvectionInterior)));
diagnostics.liquidSalinityFrameInterior = max(max(abs(liquidSalinityFrameInterior)));
diagnostics.solidSalinityFrameInterior = max(max(abs(solidSalinityFrameInterior)));
diagnostics.saltDiffusionInterior = max(max(abs(saltDiffusionInterior)));

diagnostics.baroclinicTorqueInterior = max(max(abs(baroclinicTorqueInterior)));
diagnostics.vorticityPermeabilityInterior = max(max(abs(vorticityPermeabilityInterior)));
diagnostics.vorticityDiffusionInterior = max(max(abs(vorticityDiffusionInterior)));




% compute max and average values
diagnostics.maxPorosityInterior = max(max(abs(porosityInterior)));
diagnostics.maxStreamfunctionInterior = max(max(abs(streamfunctionInterior)));
diagnostics.maxPermeabilityInterior = max(max(abs(permeabilityInterior)));

diagnostics.avPorosityInterior = nanmean(nanmean(porosityInterior));
diagnostics.avStreamfunctionInterior = nanmean(nanmean(streamfunctionInterior));
diagnostics.avPermeabilityInterior = nanmean(nanmean(permeabilityInterior));



if isfield(inputs, 'domain_width')
    diagnostics.fullWidth = str2double(inputs.domain_width)*2;
elseif  isfield(inputs, 'domain_length')
    diagnostics.fullWidth = str2double(inputs.domain_length)*2;
end

% Compute some vertical profiles and save them to a file

% Compute min and max of chi(x).
min_chi = min(chi, [], 2);
max_chi = max(chi, [], 2);

min_psi = min(psi, [], 2);
max_psi = min(psi, [], 2);

z = Y(:, 1);

save([output_dir, 'verticalProfiles.mat'], 'z', 'min_chi', 'max_chi', ...
    'min_psi', 'max_psi');
if  makePlots
    figure();
    hold on;
    plot(z, min_chi);
    plot(z, max_chi);
    hold off;
    
    legend('min($\chi$)', 'max($\chi$)');
end


if  makePlots
    hSummary = figure();
    set(hSummary, 'Position', [200 200 1000 800]);
    colormap(flipud(parula));
    
    % Compute interface positions
    y_avChi = dx*median_avPorosity_Yi;
    y_ml = dx*median_mushyLiquid_Yi ;
    y_solidMush = dx*median_solid_Yi;
    
    xFullWidth = [X(1,1) X(1,end)];
    
    chanSide = 1;
    %diagnostics.channelWidth, diagnostics.channelHeight
    y_chan_top = y_ml +  diagnostics.channelHeight;
    y_chan = y_ml + 0.5*diagnostics.channelHeight;
    
    y_psi_penetration = Y_psi_penetration*dx;
    y_chi_confinement = Y_chi_confinement*dx;
    
    if chanSide > 0
        xChanWidth = [xFullWidth(end)-diagnostics.channelWidth xFullWidth(end)];
        xChanCentre = xFullWidth(end);
    else
        xChanWidth = [xFullWidth(1) xFullWidth(1)+diagnostics.channelWidth];
        xChanCentre = xFullWidth(1);
    end
    
    xMiddle = mean(xFullWidth);
    xQuarter = xFullWidth(1) + 0.25*(xFullWidth(end) - xFullWidth(1));
    xThreeQuarter = xFullWidth(1) + 0.75*(xFullWidth(end) - xFullWidth(1));
    xTwoThirds = xFullWidth(1) + (2/3)*(xFullWidth(end) - xFullWidth(1));
    
    pcolor(X, Y, chi);
    cbar = colorbar('Location', 'westoutside');
    cbar.Label.String = '\chi';
    
    hold on;
    %mushLiquid = plot(xFullWidth, [y_ml y_ml], 'r-', 'LineWidth', 2);
    % solidMush = plot(xFullWidth, [y_solidMush y_solidMush], 'k-', 'LineWidth', 2);
    avChi = plot(xFullWidth, [y_avChi y_avChi], 'r--', 'LineWidth', 2);
    
    hPlot = plot([xMiddle xMiddle], [y_ml y_avChi], 'k--');
    HPlot = plot([xQuarter xQuarter], [y_ml y_solidMush], 'k:');
    
    PsiPenetrationPlot = plot([xThreeQuarter xThreeQuarter], [y_ml y_psi_penetration], 'blue--');
    ChiConfinementPlot = plot([xTwoThirds xTwoThirds], [y_ml y_chi_confinement], 'cyan--');
    
    
    chanWidth = plot(xChanWidth, [y_chan y_chan], 'r-', 'LineWidth',2);
    
    chanHeight = plot([xChanCentre xChanCentre], [y_ml y_chan_top], 'r:');
    
    hold off;
    
    daspect([1 1 1]);
    
    legend([avChi, chanWidth, chanHeight, hPlot, HPlot, PsiPenetrationPlot, ChiConfinementPlot],...
        {'Median height where $\chi = \bar{\chi}$', ...
        '$a$', '$h_{channel}$', '$h$', '$H$', '$h(\psi = 0.1 \psi_{max})$', '$h(\chi=0.1)$'}, ...
        'Location', 'eastoutside');
    
    
end


catch e
    fprintf('%s Some issue with analysis \n', outPref);
    fprintf(1,'%s The identifier was:\n%s \n',outPref, e.identifier);
    fprintf(1,'%s There was an error! The message was:\n%s \n',outPref, e.message);
    
end

end % end complicated analysis

% Write out diagnostics as a csv file
t = struct2table(diagnostics);
writetable(t, [output_dir, csvFileName])



end