function plotFile = getFinalPlotFile(output_dir, plot_prefix)
if nargin < 2
    plot_prefix = '';
    
    if output_dir(end) ~= '/'
    output_dir = [output_dir, '/'];
    end
end

plotFile = getNPlotFile(output_dir,plot_prefix, -1);

% dim = 2;
% subcycled = true;
% frames = [];
% 
% %plot_prefix = 'AMRConvergenceTest';
% 
% search = [output_dir, plot_prefix, '*'];
% files = dir(search);
% 
% actual_plot_prefix = '';
% for file_i = 1:length(files)
%    % fprintf('   > filename: %s \n', files(file_i).name);
%     
%     if length(strfind(files(file_i).name, '.2d.hdf5')) > 0 && ...
%             length(strfind(files(file_i).name, 'chk')) == 0
%         
%         
%         fname = strrep(files(file_i).name, '.2d.hdf5', '');
%         %fprintf('   > fname: %s \n', fname);
% 
%         parts = strsplit(fname, '-');
% 
%         if length(parts) > 1
%             frames(end+1) = str2num(parts{end});
% 
%            % fprintf('   > frame: %d \n', str2num(parts{end}));
% 
%             actual_plot_prefix = strrep(fname, parts{end}, '');
% 
% 
%            % fprintf('   > plot pref: %s \n', actual_plot_prefix);
%         end
%     
%     end
%     
% end
% 
% 
% frame = max(frames);
% %frame = 0; % testing
% 
% fprintf('    > determined plot prefix = %s and frame = %d \n', actual_plot_prefix, frame);
% 
% if length(frame) > 0 && length(actual_plot_prefix) > 0
%     plotFile = MushyLayerOutput(dim, frame, output_dir, actual_plot_prefix, subcycled);
% else
%     plotFile = []; 
% end

end