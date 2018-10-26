function plotFile = getInitialPlotFile(output_dir, plot_prefix)
if nargin < 2
    plot_prefix = '';
    
    if output_dir(end) ~= '/'
    output_dir = [output_dir, '/'];
    end
end

plotFile = getNPlotFile(output_dir,plot_prefix, 0);

end