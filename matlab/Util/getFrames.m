function [actual_plot_prefix, frames] = getFrames(output_dir, plot_prefix)
frames = [];

search = [fullfile(output_dir, plot_prefix), '*'];
files = dir(search);

actual_plot_prefix = '';
for file_i = 1:length(files)
   % fprintf('   > filename: %s \n', files(file_i).name);
   full_fname = files(file_i).name;
    
    if length(strfind(full_fname, '.2d.hdf5')) > 0 && ...
            length(strfind(full_fname, 'chk')) == 0
        
        
        fname = strrep(full_fname, '.2d.hdf5', '');
        %fprintf('   > fname: %s \n', fname);
        
        fprintf('Doing regex match on filename "%s" \n', fname);
       
        % Get final number - frame
        matchStr = regexp(fname,'\d+$','match');
        
        thisframe = str2num(matchStr{1});
        
        frames(end+1) = thisframe;
        actual_plot_prefix = strrep(fname, matchStr{1}, '');

    
    end
    
end
end