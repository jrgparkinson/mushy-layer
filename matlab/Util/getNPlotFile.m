% General function to get either first/last/any other plot file
% Set frame = 0 to get first plot file in folder
% Set frame = -1 to get last plot file in folder
% Use any other to get the specific frame

function plotFile = getNPlotFile(output_dir, plot_prefix, frame)
if nargin < 2
    plot_prefix = '';
    
    
end

if nargin < 3
    frame = -1;
end

if output_dir(end) ~= '/'
    output_dir = [output_dir, '/'];
end

dim = 2;
subcycled = true;


[actual_plot_prefix, frames] = getFrames(output_dir, plot_prefix);

% Default - assume the user has specified the frame they actually want
desiredFrame = frame;

if frame == 0
    % Get the first frame
    desiredFrame = min(frames);
elseif frame == -1
    % Get the last frame
    desiredFrame = max(frames);  
end

%frame = max(frames);
%frame = 0; % testing

fprintf('    > determined plot prefix = %s and frame = %d \n', actual_plot_prefix, desiredFrame);

if length(desiredFrame) > 0 && length(actual_plot_prefix) > 0
    plotFile = MushyLayerOutput(dim, desiredFrame, output_dir, actual_plot_prefix, subcycled);
else
    plotFile = []; 
end

end