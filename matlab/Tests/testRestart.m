% Test that stopping and restarting gives the same result as running
% a simulation continuously

%output_dir = '/Users/parkinsonj/Dropbox/DPhil/Data/'; %OS X
continuous_dir = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/output/'; %Linux
restart_dir = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/output-restart/'; %Linux
prefix = 'plt';
subcycled = true;
dim = 2;

% Whether test passes or not
success = true;


for frame = 50:60
    mlCont = MushyLayerOutput(dim, frame, continuous_dir, prefix, subcycled);
    mlRestart = MushyLayerOutput(dim, frame, restart_dir, prefix, subcycled);
 
    varToCheck = mlCont.mlComps.yDarcyVel;
    
    chCompare = ChomboCompare(mlCont);
    
    
    %figure();
     %mlCont.plot(varToCheck);
     %pause;
     
     %figure();
     %mlRestart.plot(varToCheck);
     %pause;
     
     chCompare.diff(mlRestart, [varToCheck], [varToCheck]);
     
      %figure();
     %mlRestart.plot(varToCheck);
     %pause;
     
     %figure();
     % chCompare.diff_output.plot(varToCheck);
     
     %close all;
      
      [L1, L2, Max, Sum] = chCompare.err(mlRestart, varToCheck);
      
      fprintf('Sum error for component %d at step %d = %f \n', varToCheck, frame, Sum);
    
      if Sum ~= 0.0
          fprintf('WARNING - difference detected \n ');
          success = false;
      end
      
end

if success
    fprintf('Test passed \n');
else
    fprintf('Test failed \n');
end

