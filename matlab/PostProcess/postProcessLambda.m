function  postProcessLambda(  )
close all;
clear all;

Da = {'0.0001', '0.0005', '0.006', '0.004', '0.002', '0.001', '0.01'};

base_dir = '/home/parkinsonjl/mnt/sharedStorage/AMRConvergenceTest/ConvectionDB/';

runs = {'Uniform-convectionDB-64--0', 'VariableMesh2SubcycleRefluxFreestream0.45-convectionDB-64-ref2--0'};
runs_short = {'Uniform', 'Variable Mesh'};

saveLoc = [base_dir, 'maxLambda.mat'];

%1st component is uniform, 2nd AMR
maxLambda = NaN*ones(length(Da), 2);
deltaU = NaN*ones(length(Da), 2);
load(saveLoc);




for i=1:length(Da)
    for j=1:length(runs)
        
        fprintf('%s\n', runs{j});
        
    if isnan(deltaU(i, j) )
    
       folder =  [base_dir, 'chi0.4-Da',Da{i},'-Ra1e4/',runs{j},'/'];

       plt = getFinalPlotFile(folder);

       lambda = plt.dataForComp(plt.components.lambda);

       u_ad = plt.dataForComp(plt.components.xAdvectionvelocity);
       u_darcy = plt.dataForComp(plt.components.xDarcyvelocity);

       
       deltaU(i, j) = abs(max(max(u_ad)) - max(max(u_darcy)));
       maxLambda(i, j) = max(max(abs(lambda)));
   
    
    end
    
    fprintf('Da = %s, max(abs(lambda)) = %1.3f \n', Da{i}, maxLambda(i,j));
  
    
    end
    
    
    
   
    
end

save(saveLoc, 'maxLambda', 'Da', 'deltaU');
fprintf('Saved data to %s \n', saveLoc);

Da_num = NaN*ones(length(Da), 1);

for i = 1:length(Da)
    Da_num(i) = str2num(Da{i});
end

% Post process
h = figure();
set(h, 'position', [200 200 1000 1000]);
m = 2;
n = 2;

subplot(m, n, 1)
plot(Da_num, maxLambda, 'x');
%legend(runs_short, 'Location', 'eastoutside')

xlabel('Da');
ylabel('max(abs($\lambda$))');


subplot(m, n, 2)
plot(log10(Da_num), log10(maxLambda), 'x');
legend(runs_short, 'Location', 'eastoutside')

xlabel('log10(Da)');
ylabel('log10(max(abs($\lambda$)))');


subplot(m, n, 3)
plot(Da_num, deltaU, 'x');
%legend(runs_short, 'Location',  'eastoutside')

xlabel('Da');
ylabel('$max(U_{ad}) - max(U_{cc})$');



subplot(m, n, 4)
plot(log10(Da_num), log10(deltaU), 'x');
legend(runs_short, 'Location',  'eastoutside')

xlabel('log10(Da)');
ylabel('log10($max(U_{ad}) - max(U_{cc})$)');
end

