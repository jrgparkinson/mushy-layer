% Be careful with this script - it deletes files.

function deleteOldFiles

base_folder = '/media/parkinsonjl/DATA/optimalStates-Brinkman/';
base_folder = '/home/parkinsonjl/mnt/sharedStorage/optimalStates-Brinkman/';

folders = dir(base_folder);

for i=1:length(folders)
    thisFolder = [base_folder, folders(i).name];
    
    fprintf([folders(i).name, '\n']);
    
    if length(folders(i).name) < 5
        continue
    end
    
    theseFiles = dir(thisFolder);
    
    frames = [];
    
    for j=1:length(theseFiles)
        n = theseFiles(j).name;
        if length(n) > 10
            if n(end-7:end) == '.2d.hdf5'
                
                frames(end+1) = getFrame(n);
                % fprintf([frame, '\n']);
            end
        end
    end
    
    finalFrame = max(frames);
    
    %    fprintf('%d \n', finalFrame);
    %
    %    for j=1:length(frames)
    %    fprintf('  %d \n', frames(j));
    %    end
    %
    %    fprintf('\n');
    
    for j=1:length(theseFiles)
        n = theseFiles(j).name;
        f = getFrame(n);
        
        if f > 0
            
            if sum(f == finalFrame) > 0
                fprintf('Keep %d \n', f);
            else
                fprintf('Delete %d \n', f);
                
                oldFile = [thisFolder, '/', n];
                
                if n(1:2) ~= 'zz'
                
                    %Temporarily, move to a new name

                    newFile = [thisFolder, '/zz', n];

                    fprintf([oldFile, ' to ', newFile, '\n']);
                    movefile(oldFile, newFile);
                
                else
                    delete(oldFile);
                
                end
                
                
            end
        
        end
    end
    
    
end % end loop over folders


end % end function


function f = getFrame(n)

if length(n) < 15
    f = -1;
else
    f = str2num(n(end-13:end-8));
end

end