close all;
%clear all;

%cd('../../convection-in-sea-ice/SeaIceModel/execAdvectDiffuse');
filename = 'bm1-512pts-steady.2d.hdf5';
%filename = 'plt000050.2d.hdf5';

plotfile = Plotfile(filename);
%plotfile.resize_data(plotfile.get_finest_dx());

component = 5;
figure(1);
plotfile.plot(component);

