% This post processing script will loop through a timeseries and remove
% failed timesteps, updating the timestep and time for future output files
% Also, if we dropped a checkpoint file then remove it

test_run = false;


time = 0;
step = 0;

true_time = 0;
true_step = 0;

ignored_time = 0;
ignored_steps = 0;

output_dir = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/output/';

plot_prefix = 'plt';
chk_prefix = 'chk';
subcycled = true;
dim = 2;

filesExist = true;

%ml = MushyLayerOutput(dim, step, output_dir, plot_prefix, subcycled);

while (filesExist)
    
fname = ChomboOutput.getFilename(output_dir, plot_prefix, dim, step);
chkname = ChomboOutput.getFilename(output_dir, chk_prefix, dim, step);

if ~exist(fname)
    filesExist = false;
    continue;
end

 fileinfo = hdf5info(fname);

time = h5readatt(fname,'/level_0', 'time');
dt = h5readatt(fname,'/level_0', 'dt');
timestepFailed = h5readatt(fname,'/level_0', 'timestepFailed');

if (timestepFailed)
    % delete file
    
    fprintf('Delete file %s \n', fname);
    ignored_time = ignored_time + dt;
    ignored_steps = ignored_steps + 1;
    
    if test_run ~= true
        delete (fname);
        
        if exist(chkname)
            delete(chkname);
        end
    end
else
    % set step and time correctly
    
    true_time = time - ignored_time;
    true_step = step - ignored_steps;
    
    % Preserve formatting
    str_true_step = num2str(true_step);
    if length(str_true_step) < length(num2str(step))
        str_true_step = ['0', str_true_step];
    end
    
    correct_name = strrep(fname, num2str(step),str_true_step);
    correct_chkname = strrep(chkname, num2str(step),str_true_step);
    
    fprintf('time %f -> %f, filename %s->%s \n', time, true_time, fname, correct_name);
    
    if test_run ~= true
        h5writeatt(fname,'/level_0', 'time', true_time);
        if step ~= true_step
            movefile(fname,correct_name);
           
        end
        
        %Check if there are checkpoint files we should move
            if exist(chkname)
                h5writeatt(chkname,'/level_0', 'time', true_time);
        if step ~= true_step
            movefile(chkname,correct_chkname);
           
        end
            end
    end
    
     
    
end

step = step + 1;

end
            
