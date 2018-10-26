 clear all;
 
 output_dir = '/Users/parkinsonj/Dropbox/DPhil/Data/';
 

 cd(output_dir);
 
 % This MUST match exactly the file reference at the 
 % top of out_xxxx.info
 %file = './katz-bm1/out_21034';
 file = './katz/out_0201';
 
 ml = loadMushyLayerOutput(file);
 
 figure();
 mushyLayerDisplay(ml, 'logchi');
 
% mushyLayerDisplayMovie(ml, 'schlogchi');
 