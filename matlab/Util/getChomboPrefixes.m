function plot_prefixes = getChomboPrefixes(outputDir, includeChkFiles)

if nargin < 2
    includeChkFiles = true;
end

files = dir(outputDir);

plot_prefixes = {};

for i =1:length(files)
    filename = files(i).name;
    if (length(filename) > 5)
        if strcmp( filename(end-4:end) , '.hdf5')
            thisPrefix = regexprep(filename,'\.\dd\.hdf5','');
            thisPrefix = regexprep(thisPrefix,'\d+$','');
            
            % check if we already have this prefix
            a = strfind(plot_prefixes, thisPrefix);
            prefixExists = false;
            for pi = 1:length(a)
                if a{pi} == 1
                    prefixExists = true;
                end
            end
            if ~prefixExists
                if includeChkFiles || ~ strcmp(thisPrefix(1:3), 'chk')
                    plot_prefixes{end+1} = thisPrefix;
                end
            end
            
        end
    end
end

end