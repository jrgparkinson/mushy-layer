dir = '/home/parkinsonjl/gyre/katz/';
dir = '/home/parkinsonjl/mushylayer/benchmark/bm2_src/output/';
dir = '/home/parkinsonjl/mushylayer/benchmark/bm2_src/output-vertical/';
dir = '/home/parkinsonjl/mushylayer/benchmark/bm3_src/bm3_Ra100/';
dir = './katz/fig9-5/'; cd '/home/parkinsonjl/convection-in-sea-ice/test/';
%dir = '/home/parkinsonjl/mushylayer/benchmark/bm2_src/output-frameAdv/';
close all;

%ml = loadMushyLayerOutput(filename);

figure();
cmin = -4;
cmax = 4;
for i=30:2:1000
    name = sprintf('%04d', i);
    filename = [dir, 'out_', name];
    
    %ml = loadConvectTestOutput(filename);
    ml = loadMushyLayerOutput(filename);
    full_field = log(ml.field.Pi);
    full_field = ml.field.phi;
    plotting_field = full_field(2:end-1, 2:end-1); %remove ghost cells
    %ml.field.phi;
    
    %ml.par.RaT
    
    imagesc(plotting_field); set(gca, 'YDir', 'normal');
    
    colorbar();
    title([name, ' t=', num2str(ml.par.t_new)]);
    
%     if i == 0
%         cmin = min(min(plotting_field));
%         cmax = max(max(plotting_field));
%         caxis([cmin cmax]);
%     end
    
    pause(0.1);
end





