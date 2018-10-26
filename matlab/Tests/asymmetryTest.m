close all;

data_dir = '/Users/parkinsonj/Dropbox/DPhil/Data/';
dim = 2;

% ======================================================================
% Choose which file/frame to analyse
% ======================================================================
frame = 0;
rotation = true;

% output = MushyLayerOutput(dim, frame, data_dir, 'oneStep/plot');
% output = MushyLayerOutput(dim, frame, data_dir, 'oneStep/plotCentreDiv');
% output = MushyLayerOutput(dim, frame, data_dir, 'oneStep/plotFixedGrid');
% output = MushyLayerOutput(dim, frame, data_dir, 'oneStep/plotFixedU');
% output = MushyLayerOutput(dim, frame, data_dir, 'oneStep/plotFixedP');

frame = 8;
output = MushyLayerOutput(dim, frame, data_dir, 'bm2/bm2-Ra100-64pts');
output = MushyLayerOutput(dim, frame, data_dir, 'bm2-symm/bm2-Ra100-64pts'); %symmetric perturbation
rotation = false;
%output = MushyLayerOutput(dim, frame, data_dir, 'bm2-fixed/bm2-fixed-Ra100-64pts');

%output.plot(output.mlComps.pressureErr);

% ======================================================================
%Choose which component we want to analyse the symmetry of:
% ======================================================================
component = output.mlComps.divUstarErr;
reflection = false;

component = output.mlComps.pressureErr;
reflection = false;

% component = output.mlComps.xGradPErr;
% reflection = true;

% component = output.mlComps.zFluidVel;
% reflection = true;

% component = output.mlComps.zFluidVelErr;
% reflection = true;


component = output.mlComps.theta;
reflection = true;

% component = output.mlComps.thetaForcingErr;
% reflection = true;


% ======================================================================
% Rotate field by 180 degrees and subtract from original
% ======================================================================
err = output.dataForComp(component);
err_rotated = imrotate(err, 180);

if reflection
    asymmetry = err + err_rotated;
else
    asymmetry = err - err_rotated;
end

asymmetry_rotated = imrotate(asymmetry, 180);

if reflection
    asymmetry_two = asymmetry_rotated - asymmetry;
else
    asymmetry_two = asymmetry_rotated + asymmetry;
end

figure(1);
imagesc(err);
title('Field');
colorbar();
set(gca, 'YDir', 'normal');

figure(3);
imagesc(asymmetry);
title('Asymmetry in field');
colorbar();
set(gca, 'YDir', 'normal');


% figure(5);
% imagesc(asymmetry_two);
% title('Asymmetry in asymmetry field');
% colorbar();
% set(gca, 'YDir', 'normal');
% 
