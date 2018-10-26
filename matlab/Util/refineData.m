% Take a chombo output file, refine it for specified variables, and write
% out the results to a csv file
clear all; close all;

%Specify these variables:
dim = 2;
output_dir = '/home/parkinsonjl/gyre/bm2-fixedPorosityField-Ra40/64-1/';
subcycled = false;

old_res = 64;

ml =  MushyLayerOutput(dim, -1, output_dir, 'bm2-fixedPorosityField-Ra40-64pts01687', subcycled,'bicubic'); % 'bicubic'

comps = [ml.mlComps.theta];
comp_filenames = {'theta'};
file_suffix = 'dat';
refinement = 2;
includeGhost = false;

% Don't change anything else

for refinement=[1,2,4, 8]
new_res = old_res*refinement;

lev_i = 1;

   lev = ml.levelArray(lev_i);
   for comp_i = 1:length(comps)
       comp = comps(comp_i);
       
       %dataUnrefined = lev.dataForComp(comp);
       dataRefined = lev.dataForComp(comp, includeGhost, refinement);
       %dataRefined = smoothn(dataRefined, 0.1, 'robust');
       dataRefined=inpaint_nans(dataRefined);
       
       %Fill nan cells with extrapolation from interior
       
       % Write to file
       filename = [output_dir, comp_filenames{comp_i}, '-', num2str(new_res),'.', file_suffix]
       csvwrite(filename, dataRefined);
       
   end
   
   figure();
   imagesc(dataRefined);
   
end





