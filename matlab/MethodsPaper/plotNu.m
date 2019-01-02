close all;

plotSingle = true;

if plotSingle
    
    hSingle = figure();
    hSingle.Position(3) = 1000;
    hSingle.Position(1) = 200;
    
    folder = '/home/parkinsonjl/mushy-layer/execSubcycle/Nu/';
    %folder = ['/home/parkinsonjl/mnt/sharedStorage/TestDiffusiveTimescale/ConvectionDB-cfl0.29', ...
    %    '/chi0.4-Da1.0e-06-Ra1.0e+09/VariableMesh2SubcycleRefluxFreestream0.95-ref2-convectionDB-32--0/'];
    
    folders = {'/home/parkinsonjl/mushy-layer/execSubcycle/Nu/method0/', ...
        '/home/parkinsonjl/mushy-layer/execSubcycle/Nu/method1/'};
    
    UniformPrefix = 'Uniform-convectionDB-';
AMRPrefix = 'VariableMesh2SubcycleRefluxFreestream0.95-ref2-convectionDB-';


    remoteBaseFolder = ['/home/parkinsonjl/mnt/sharedStorage/TestDiffusiveTimescale/ConvectionDB-cfl0.09/'];
    
    folders = {[remoteBaseFolder, 'chi0.4-Da1.0e-06-Ra1.0e+07/', UniformPrefix,'128--0']...
        [remoteBaseFolder, 'chi0.4-Da1.0e-06-Ra1.0e+08/', UniformPrefix, '128--0'], ...
        [remoteBaseFolder, 'chi0.4-Da1.0e-06-Ra1.0e+09/', UniformPrefix, '128--0']};
    
    NuLeBars = 12.9;
    NuLeBars = 3.17;
    NuLeBars = 5.24;
    NuLeBars = 3.07;
    
    NuLeBars = [1.01 1.41 3.17 5.24 1.08 3.07 12.9];
    
    % folder = ['/home/parkinsonjl/mnt/sharedStorage/TestDiffusiveTimescale/', ...
    %    'ConvectionDB-cfl0.26/chi0.4-Da1.0e-06-Ra1.0e+09/Uniform-convectionDB-128--0/'];
    
   
    m = 1; n = 2;
    
    leg = {};
    
    subplot(m, n, 1);
    
    hold on
    for i=1:length(folders)
         d = getDiagnostics(fullfile(folders{i}, 'diagnostics.out'));
        plot(d.time, d.Nusselt);
        
        leg{end+1} = folders{i}(end-10:end);
    end
    
    xlabel('t');
    ylabel('Nu');
    box on;
    
    % Recent average
%     Nu = d.Nusselt;
%     RecentNu = Nu(round(0.8*end):end);
%     avNu = mean(RecentNu);
%     maxNu = max(d.Nusselt);
%     minNu = min(d.Nusselt);
%     %title(['Average Nu = ', num2str(avNu), ', max Nu = ',num2str(maxNu),', min Nu = ', num2str(minNu)]);
%     stdev = std(RecentNu);
%     titleString = sprintf('Average Nu = %1.2f +/- %1.3f', avNu, stdev);
%     %title(['Average Nu = ', num2str(avNu), ' +/- ', num2str(stdev)]);
%     title(titleString);
    
 
    
   % plot(d.time, d.time*0 + avNu, '--');
   % plot(d.time, d.NusseltLeft);
   % plot(d.time, d.NusseltRight);
   xl = xlim;
   for Li = 1:length(NuLeBars)
    plot(xl, xl*0 + NuLeBars(Li), '--');
   end
    hold off;
    
    leg{end+1} = 'Nu Le Bars';
    legend(leg, 'Location', 'southeast');
    
    %ylim([avNu*0.9 avNu*1.05]);
    
%     subplot(m, n, 2);
%     hold on;
%     
%     plot(d.time, d.MaxVHalf);
%     plot(d.time, d.MaxUHalf);
%     
%     hold off;
%     
%     box on;
%     xlabel('time');
%     ylabel('max vel');
%     legend({'Max V', 'Max U'});
%     
%     avV = mean(d.MaxVHalf(round(0.8*end):end));
%     avU = mean(d.MaxUHalf(round(0.8*end):end));
%     
%     title({['Av(Max U) = ', num2str(avU)], ['Av(Max(V)) = ', num2str(avV)]});
%     
    
else
    
    chi = 0.4;
    Da = 1e-6;
    Ra = 1e9;
    
    % subfolder = 'chi0.4-Da1.0e-06-Ra1.0e+09';
    subfolder = sprintf('chi%1.1f-Da%1.1e-Ra%1.1e', chi, Da, Ra);
    res = [64,128];
    leg = {};
    
    %d = getDiagnostics('/home/parkinsonjl/mnt/sharedStorage/TestDiffusiveTimescale/ConvectionDB-cfl0.15/chi0.4-Da1.0e-06-Ra1.0e+09/Uniform-convectionDB-128--0/diagnostics.out');
    
    figure();
    hold on;
    for i=1:length(res)
        d = getDiagnostics(['/home/parkinsonjl/mnt/sharedStorage/TestDiffusiveTimescale/', ...
            'ConvectionDB-cfl0.26/', subfolder, '/Uniform-convectionDB-',num2str(res(i)),'--0/diagnostics.out']);
        plot(d.time, d.Nusselt)
        leg{end+1} = num2str(res(i));
    end
    
    LeBarsNu = getLeBarsNu(chi, Da, Ra);
    plot(d.time, LeBarsNu + d.time*0, '--');
    leg{end+1} = 'Le Bars \& Worster (2006)';
    
    hold off;
    
    title(subfolder);
    xlabel('t');
    ylabel('Nu');
    box on;
    
    legend(leg, 'Location', 'southeast');
end