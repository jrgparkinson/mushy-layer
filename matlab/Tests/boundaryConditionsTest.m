%Load an output file
close all;

data_dir = '/Users/parkinsonj/Dropbox/DPhil/Data/bm2-fixed/';
dim = 2;
frame = 0;

output = MushyLayerOutput(dim, frame, data_dir, 'bm2-fixed-Ra100-64pts'); %, 'bilinear')

%Get theta
theta = output.dataForComp(output.mlComps.theta);
[X,Y] = output.grid();
dx = output.finest_dx();


% contourf(X, Y, theta);

[dtheta_dx, ~] = gradient2order(theta, dx, dx);

figure(1);
imagesc(dtheta_dx);
colorbar();
title('d theta/dx');
