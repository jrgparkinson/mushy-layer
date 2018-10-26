function [filenames, frames] = getFiles(d, prefix)
files = dir([d, prefix, '*']);

filenames = {};
frames = [];

for i=1:length(files)
    file = files(i);
    filenames{end+1} = file.name;
    thisFrame = strrep(file.name, 'toPlot', '');
    thisFrame = strrep(thisFrame, '.mat', '');
    frames(end+1) = str2num(thisFrame);
end


%frames = sort(frames);


%frames = files;
end