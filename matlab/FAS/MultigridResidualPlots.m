% script to plot some residuals in order to illustrate how multigrid works
function residualPlots
clear all;
close all;

base_dir = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/';
residuals = {'down8', 'down4', 'down2', ...
    'exactSolve', 'up2', 'up4', 'up8'};

%residuals = {'down8', 'down4'};

m = 1;
n = length(residuals);

%h = figure();
%set(h, 'Position', [100 100 1600 400]);




for i = 1:length(residuals)
    
    
    filename = [base_dir, residuals{i}];
    
    output = ChomboOutput(2,-1,base_dir,residuals{i});
    
    data = output.dataForComp(output.components.component_0);
    dataMin = min(min(data));
    dataMax =  max(max(data));
    
    
    
    %ax(i) = subplot(m, n, i);
    fig = figure();
    set(fig, 'Position', [300 300 600 900]);
    c(i) = output.plot(output.components.component_0);
    
    colormap(bluewhitered);
    c(i).Location = 'southoutside';
    set(c(i),'fontsize',44);
   
    %caxis([dataMin dataMax]);
    
   % if i == 6
   %    colormap(bluewhitered); 
   % end
%     if dataMin == dataMax
%         colormap([1 1 1]);
%     else
%         colormap(bluewhitered);
%     end
    
   
    h = gca;
    set(gca, 'XTickLabel', '');
    set(gca, 'YTickLabel', '');    
    
    print(fig, residuals{i}, '-dpng');
    
  
end


end


