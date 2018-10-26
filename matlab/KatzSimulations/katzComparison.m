%Code to compare an output file from Rich's model with output file from my model

dim = 2;
plot_prefix = 'bm2-Ra10';
plot_prefix = 'bm2-Ra100-32pts-steady';
close all;
dir_parkinson =  '/home/parkinsonjl/convection-in-sea-ice/test/bm2-Ra100/';
dir_katz = '/home/parkinsonjl/mushylayer/benchmark/bm2_src/output-vertical-Ra100/';
dir_katz_frameAdv = '/home/parkinsonjl/mushylayer/benchmark/bm2_src/output-frameAdv/';

%dir_katz = '/home/parkinsonjl/mushylayer/benchmark/bm3_src/bm3_Ra100';
%dir_katz = '/home/parkinsonjl/mushylayer/benchmark/bm2_src/output-vertical-Ra0/';
%filename_parkinson = '/home/parkinsonjl/convection-in-sea-ice/test/bm2-Ra100/bm2-Ra100-32pts-steady.2d.hdf5';

%Ra 50

% Ra 100
filename_katz = '/home/parkinsonjl/mushylayer/benchmark/bm2_src/output/out_0170';
%filename_katz = '/home/parkinsonjl/mushylayer/benchmark/bm2_src/output-frameAdv/out_0389';

% Ra 200

ml_katz = loadConvectTestOutput(filename_katz);
ml_park = MushyLayerOutput(dim, -1, dir_parkinson, plot_prefix);

T_katz = ml_katz.field.T(2:end-1, 2:end-1);
u_katz = ml_katz.field.u(2:end-1, 2:end-1);
v_katz = ml_katz.field.w(2:end-1, 2:end-1);
T_park = ml_park.dataForComp(ml_park.mlComps.theta);
u_park = ml_park.dataForComp(ml_park.mlComps.xFluidVel);
v_park = ml_park.dataForComp(ml_park.mlComps.zFluidVel);

% figure(1);
% imagesc(T_katz); colorbar();  set(gca, 'YDir', 'normal'); title('Katz temperature');
% 
% figure(2);
% imagesc(T_park); colorbar();  set(gca, 'YDir', 'normal'); title('Parkinson temperature');

% figure(3);
% diff = T_katz - T_park;
% imagesc(diff);
% colorbar();
% set(gca, 'YDir', 'normal'); title('Difference');

fig = figure(1);
park = T_park;
set(fig, 'Position', [400 400 1400 400])
crange = [min(min(park)) max(max(park))];
for i=0:100
    name = sprintf('%04d', i);
    filename = [dir_katz, 'out_', name];
    ml = loadConvectTestOutput(filename);
    with_ghost = ml.field.T; 
    katz = ml.field.T(2:end-1, 2:end-1); %remove ghost cells
    
    diff = katz-park;
    
    subplot(1,3,1); 
    imagesc(katz);  set(gca, 'YDir', 'normal'); colorbar(); caxis(crange);
    title(['Katz plot']);
    
    
    subplot(1,3,2);
    imagesc(park);  set(gca, 'YDir', 'normal'); colorbar(); caxis(crange);
    title(['Parkinson plot']);
    
    subplot(1,3,3);
    imagesc(diff);
    set(gca, 'YDir', 'normal');
    colorbar();
    title(['Katz - Park (variable = T), frame = ', name]);
    
    pause(0.01);
end



