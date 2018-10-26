% Get the runs contained in a directory (ignoring similar runs with
% different resolutions)

function runs = getRuns(outputDir)

folders = dir(outputDir);

runs = {};

for i =1:length(folders)
    folderName = folders(i).name;
    if (length(folderName) > 5)
        if strcmp( folderName(end-1:end) , '-0')
            %thisPrefix = regexprep(folderName,'\.\dd\.hdf5','');
            %thisPrefix = regexprep(folderName,'\d+-(ref\d)--0$','$1');
            thisPrefix = regexprep(folderName,'\d+--0$','');
            
            % check if we already have this prefix
            a = strfind(runs, thisPrefix);
            prefixExists = false;
            for pi = 1:length(a)
                if a{pi} == 1
                    prefixExists = true;
                end
            end
            if ~prefixExists
                runs{end+1} = thisPrefix;
            end
            
        end
    end
end

end
