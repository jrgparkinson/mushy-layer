function plotLambda

close all;

%base_dir = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/';
base_dir = '/home/parkinsonjl/mnt/sharedStorage/lambdaConvergenceTest/';

if 7~=exist(base_dir,'dir')
   fprintf('Base directory does not exist - are external drives mounted? \n'); 
   return;
end


foldersVecUpper = {'coupledEta0.0NoBuoyancy',...
    'coupledEta0.95NoBuoyancy'};

%foldersVecUpper = {'coupledEta0.0',...
%    'coupledEta0.95', 'coupledEta0.95AllPorous'}; 
foldersVecUpper = {'coupledEta0.0middle',...
    'coupledEta0.95', 'coupledEta0.95AllPorous'}; 


foldersVecLower = {'coupledEta0.0AllPorousmiddle',...
    'coupledEta0.95MinRegridTime0.3customLambdaTest', 'coupledEta0.95AllPorous'};


% Also try Porosity1.0Eta0.95middleuDelUfixedDt
foldersVecUpper = {'Porosity1.0Eta0.0middlefixedDt', 'coupledEta0.0middlefixedDt',   'coupledEta0.0middle'}; 
foldersVecLower = {'Porosity1.0Eta0.95middlefixedDt', 'coupledEta0.95middlefixedDt', 'coupledEta0.95middle'};


m = 2; n = 3;
resolutions = [32 64 128];
lineTypes = {'-.',':', '-', '--'};



detailedPlots = false;
summaryPlot = true;

if summaryPlot
    hSummary = figure();
hSummary.Position=[200 200 1500 800];

   m=2; n=3;
   
   noEtaFolder = 'coupledEta0.0middlefixedDt';
   etaFolder = 'coupledEta0.95middlefixedDt';
   
   
   
   convergenceMax = 0.0; convergenceSum  = 0.0;
   
   subplot(m, n, 1);
   fullFolder = [base_dir, noEtaFolder, '/Nx', num2str(resolutions(end)), '-0'];
   plotFinalState(noEtaFolder, fullFolder, convergenceMax, convergenceSum);
   title('$\eta=0.0$');
   
   subplot(m, n, 4);
   fullFolder = [base_dir, etaFolder, '/Nx', num2str(resolutions(end)), '-0'];
   plotFinalState(etaFolder, fullFolder, convergenceMax, convergenceSum);
   title('$\eta=0.95$');
   
   % Plot time series
   subplot(m, n, 2);
   steady_max_no_eta = plotMaxLambda(base_dir, noEtaFolder, resolutions);
   
   subplot(m, n, 5);
   steady_max_eta = plotMaxLambda(base_dir, etaFolder, resolutions);
   

   % Now plot convergence
   subplot(m, n, 3);
   plotSimpleSummary(resolutions, steady_max_no_eta, 1);
   
   subplot(m, n, 6);
   plotSimpleSummary(resolutions, steady_max_eta, 2);
   
    
end

if detailedPlots
    
    % Plot timeseries of lambda error
h = figure();
h.Position = [200 400 1500 800];

[steady_max_err1,steady_err_sq1] = plotRow(m, n, 1, resolutions, lineTypes, base_dir, ...
    foldersVecUpper);
[steady_max_err2,steady_err_sq2] = plotRow(m, n, 4, resolutions, lineTypes, base_dir, ...
    foldersVecLower);

% Now plot convergence of lambda error
hSummary = figure();
hSummary.Position=[200 200 1500 800];

left_color = [1 0 0];
right_color = [0 0 1];
set(hSummary,'defaultAxesColorOrder',[left_color; right_color]);

m = 1; n = 2;
plotSummary(m, n, 1, resolutions, lineTypes, steady_max_err1, steady_err_sq1, foldersVecUpper);
plotSummary(m, n, 2, resolutions, lineTypes, steady_max_err2, steady_err_sq2, foldersVecLower);


% Compute convergence rates
convergenceUpperMax = computeConvergence(resolutions, steady_max_err1);
convergenceLowerMax = computeConvergence(resolutions, steady_max_err2);
convergenceUpperSum = computeConvergence(resolutions, steady_err_sq1);
convergenceLowerSum = computeConvergence(resolutions, steady_err_sq2);

convergenceMax = [convergenceUpperMax, convergenceLowerMax];
convergenceSum = [convergenceUpperSum, convergenceLowerSum];

% Plot final states
hFinalStates = figure();
hFinalStates.Position = [0 0 1600 1000];

annotation('textbox',[0 0.8 0.1 0.1],...
    'String','Color: $\lambda-1$. Contours: temperature or porosity',...
    'EdgeColor','none', 'FitBoxToText','off', 'Interpreter', 'latex', 'FontSize', 16);


m = 2; n = 3;
allFolders = [foldersVecUpper, foldersVecLower];
for i = 1:length(allFolders)
   subplot(m, n, i);
   fullFolder = [base_dir, allFolders{i}, '/Nx', num2str(resolutions(end)), '-0'];
   
   
   if 7~=exist(fullFolder,'dir')
       fprintf('Dir does not exist: %s \n', fullFolder); 
       return;
   end
   
   plotFinalState(allFolders{i}, fullFolder, convergenceMax(i), convergenceSum(i));
end



end

end

function plotFinalState(title_str, fullFolder, convMax, convSum)

    output = getFinalPlotFile(fullFolder, 'plt');

    porosity = output.dataForComp(output.components.Porosity).';
    lambda = output.dataForComp(output.components.lambda).';
    T = output.dataForComp(output.components.Temperature).';
    
    fieldToPlot = porosity;
    if max(max(fieldToPlot)) == min(min(fieldToPlot))
        fieldToPlot = T;
    end
       
    [X,Y] = output.grid();
   
   % amrPolyshapes = output.getMeshes();
    
    % a) plot porosity and grids  
%     pcolor(X, Y, porosity);
%     

% 
%     hold on;
%     plot(amrPolyshapes, 'FaceColor', 'magenta', ...
%         'FaceAlpha',0.1, 'EdgeColor', 'magenta', 'LineWidth', 1.5)
%     hold off;
    
    % b) plot porosity contours and lambda
    % Plot rescaled lambda
    
    %lambdaScale = max( max(max(lambda)), -min(min(lambda)));
    %scaledLambda = (lambda/lambdaScale);
    %minL = min(min(scaledLambda)); maxL = max(max(scaledLambda));
    scaledLambda = lambda;
    minL = min(min(scaledLambda)); maxL = max(max(scaledLambda));
    
    for i=1:size(scaledLambda, 1)
        for j=1:size(scaledLambda, 2)
            if scaledLambda(i, j) < 0
                 scaledLambda(i, j) =  scaledLambda(i, j)/abs(minL);
            else
                 scaledLambda(i, j) =  scaledLambda(i, j)/abs(maxL);
            end
        end
    end
    
   % scaledLambda = 0.5*(scaledLambda+1);
    
    newMinL = min(min(scaledLambda)); newMaxL = max(max(scaledLambda));
     
    %pcolor(X, Y, scaledLambda);
    pcolor(X,Y,lambda);
    
    hold on;
   % [C, h] = contour(X, Y, porosity, 5);
    [C, h] = contour(X, Y, fieldToPlot, 7);
    h.LineColor = [0 0 0];
    h.LineWidth = 2.0;
    hold off;

    caxis(gca, [min(min(lambda)) max(max(lambda))]);
    colormap(gca, bluewhitered);
    cbar = colorbar('Location', 'southoutside');
    
    % Take off axis labels
    axPorosity = gca;
    axPorosity.XTick = [];
    axPorosity.YTick = [];

    % Add fine region outline
    title_str = strrep(title_str, 'coupled', '');
    longtitle(title_str, 25);

    if convMax ~= 0 && convSum ~= 0
        text(X(1, 4), Y(10,1), ['Conv: ', num2str(convMax), ' (max), ', num2str(convSum) ,' (sum)' ], 'FontSize', 14);
    end

    drawnow;

end

function plotSummary(m, n, subplot_i, resolutions, lineTypes, steady_max_err, steady_max_sq_err, folders)

subplot(m, n, subplot_i);

l = [];

x = log10(resolutions);
xlabel('log$_{10}$(Nx)'); 

yyaxis left;
ylabel('max$|\lambda-1|$');
hold on;
for i=1:length(folders)
    thisErr = steady_max_err(:, i);
    
    y = log10(thisErr);
    
l(end+1) = plot(x, y,  'LineStyle', lineTypes{i}); %'Color', 'r',

convergence = computeConvergence(resolutions, thisErr);
text(x(end), y(end), ['Rate = ', num2str(convergence)], 'FontSize', 14);

% Draw line with this gradient
conv_plot = y;
for i=2:length(conv_plot)
    conv_plot(i) = conv_plot(i-1)-convergence*(x(i)-x(i-1));
end
%plot(x, conv_plot, ':');

end

yyaxis right;
ylabel('max$|\lambda-1|^2$');
for i=1:length(folders)
    thisErr = steady_max_sq_err(:, i);
    y = log10(thisErr);
plot(x, y, 'LineStyle', lineTypes{i}); % 'Color', 'b',

% Get the gradient of the line
convergence = computeConvergence(resolutions, thisErr);
text(x(end), y(end), ['Rate = ', num2str(convergence)], 'FontSize', 14);

% Draw line with this gradient
conv_plot = y;
for i=2:length(conv_plot)
    conv_plot(i) = conv_plot(i-1)-convergence*(x(i)-x(i-1));
end
%plot(x, conv_plot, ':');


end

firstOrder = max(max(steady_max_err))*ones(1, length(resolutions));
secondOrder = max(max(steady_max_err))*ones(1, length(resolutions));

for i=2:length(resolutions)
   ratio = resolutions(i)/resolutions(i-1);
    
   firstOrder(i) = firstOrder(i-1)/ratio;
   secondOrder(i) = secondOrder(i-1)/ratio^2;
end

plot(x, log10(firstOrder), 'k-');
plot(x, log10(secondOrder), 'k-');

hold off;

 
box on;

xlim([x(1)*0.8 x(end)*1.2]);

leg = folders;
leg{end+1} = '1st order';
leg{end+1} = '2nd order';

legend(l, folders, 'Location', 'northoutside');


end

function [steady_max_err, steady_max_sq_err] = plotRow(m, n, starti, resolutions, ... 
lineTypes, base_dir, foldersVec, ...
doMaxLambda, doSumLambda, doLegend)



scaleLambdaWithVel = false;

leg = {};

subplot(m, n, starti)
hold on;


%resolutions = [32 64 128];

steady_max_err = NaN*ones(length(resolutions), length(foldersVec));
steady_max_sq_err = NaN*ones(length(resolutions), length(foldersVec));

maxMaxLambda = 1e-15;

minFinalLambda = 1e-1;

for f = 1:length(foldersVec)
    ax = gca;
    ax.ColorOrderIndex = 1;
    folder = foldersVec{f};
    lineType = lineTypes{f};
    for i=1:length(resolutions)

        fname = [folder, '/Nx', num2str(resolutions(i)), '-0'];
        
        diags = getDiagnostics([base_dir, fname, '/diagnostics.out']);

        if isfield(diags, 'LambdaMax')
            if (scaleLambdaWithVel)
                plot(diags.time, diags.LambdaMax./diags.maxVel, lineType);
            else
                plot(diags.time, diags.LambdaMax, lineType);
            end
            
            leg{end+1} = fname;
            
            steadyLambdaMax = nanmean(diags.LambdaMax(round(end*0.75):end));
            steady_max_err(i, f) = steadyLambdaMax;
            
            maxMaxLambda = max(maxMaxLambda, max(diags.LambdaMax));
            
            if min(diags.LambdaMax(end)) > 1e-10
                minFinalLambda = min(minFinalLambda, diags.LambdaMax(end));
            end
            
        end
        
        
    end 
    
    ypos = ax.Position(2);
    yheight = ax.Position(4);
    
    ax.Position = [0.08 ypos 0.2 yheight];
end

xlabel('$t$');
if (scaleLambdaWithVel)
    ylabel('max$|\lambda-1|/$max$|\mathbf{U}|$');
else
    ylabel('max$|\lambda-1|$');
end
ylim([minFinalLambda/100 maxMaxLambda*2]);

set(gca,'yscale','log')

box on;
hold off;



minFinalLambda = 1e-1;

subplot(m, n, [starti+1 starti+2]);
hold on;
for f = 1:length(foldersVec)
    ax = gca;
    ax.ColorOrderIndex = 1;
    folder = foldersVec{f};
    lineType = lineTypes{f};
   for i=1:length(resolutions)

        
        diags = getDiagnostics([base_dir, folder, '/Nx', num2str(resolutions(i)), '-0/diagnostics.out']);

        if isfield(diags, 'LambdaSum')
            plot(diags.time, diags.LambdaSum, lineType);
           % leg{end+1} =folders{i};
           if min(diags.LambdaSum(end)) > 1e-10
            minFinalLambda = min(minFinalLambda, diags.LambdaSum(end));
           end
           
           steadyLambdaMax = nanmean(diags.LambdaSum(round(end*0.75):end));
            steady_max_sq_err(i, f) = steadyLambdaMax;
            
        end
   end 
    
   ypos = ax.Position(2);
   yheight = ax.Position(4);
    
   ax.Position = [0.35 ypos 0.2 yheight];
end

xlabel('$t$');
        ylabel('sum$(\lambda-1)^2$');
  ylim([minFinalLambda/10 1e-2]);
set(gca,'yscale','log')
box on;
hold off;


l = legend(leg, 'Location', 'eastoutside');  %'parent',hpanel,'outerposition',[.25 0 .5 1]

l.Position = [(1-l.Position(3)) ypos l.Position(3) yheight];

end

function convergence = computeConvergence(resolutions, err)

numDataSets = size(err, 2);

convergence = NaN*ones(1, numDataSets);

for i = 1:numDataSets
    thisErr = squeeze(err(:, i)).';
    coeffs = polyfit(log10(resolutions), log10(thisErr), 1);
    convergence(i) = -coeffs(1);
end


end



function steady_max_err = plotMaxLambda(base_dir, folder, resolutions)

maxMaxLambda = 1e-15;
minFinalLambda = 1e-1;
scaleLambdaWithVel = false;
lineType = '-';

hold on;

for i=1:length(resolutions)

        fname = [folder, '/Nx', num2str(resolutions(i)), '-0'];
        
        diags = getDiagnostics([base_dir, fname, '/diagnostics.out']);

        if isfield(diags, 'LambdaMax')
            if (scaleLambdaWithVel)
                plot(diags.time, diags.LambdaMax./diags.maxVel, lineType);
            else
                plot(diags.time, diags.LambdaMax, lineType);
            end
            
            %leg{end+1} = fname;
            
            steadyLambdaMax = nanmean(diags.LambdaMax(round(end*0.75):end));
            steady_max_err(i) = steadyLambdaMax;
            
            maxMaxLambda = max(maxMaxLambda, max(diags.LambdaMax));
            
            if min(diags.LambdaMax(end)) > 1e-10
                minFinalLambda = min(minFinalLambda, diags.LambdaMax(end));
            end
            
        end
        
        
end 

xlabel('$t$');
if (scaleLambdaWithVel)
    ylabel('max$|\lambda-1|/$max$|\mathbf{U}|$');
else
    ylabel('max$|\lambda-1|$');
end
ylim([minFinalLambda/100 maxMaxLambda*2]);

set(gca,'yscale','log')

box on;

    hold off;

end


function plotSimpleSummary(resolutions, steady_max, order)

x = 1./resolutions;

ideal_slope = steady_max*0.9;
for i=2:length(ideal_slope)
   ideal_slope(i) = ideal_slope(i-1) * (x(i)/x(i-1))^(order) ;
end

hold on;

plot(x, steady_max, 'rx');
plot(x, ideal_slope, 'k-');

hold off;


lab = '';
if order == 1
    lab = '1st order';
elseif order == 2
    lab = '2nd order';
end
text(x(2)*1.2, ideal_slope(2)*0.8, lab, 'FontSize', 20);


set(gca, 'XScale', 'log', 'YScale', 'log');



box on;

xlabel('$\Delta x$');
ylabel('max$|\lambda(t=0.5)-1|$');


end