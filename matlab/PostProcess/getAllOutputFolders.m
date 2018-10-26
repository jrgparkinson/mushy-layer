function [ folders ] = getAllOutputFolders( data_dir )
%GETALLOUTPUTFOLDERS Summary of this function goes here
%   Detailed explanation goes here
folders = {};

% Loop over Da_folders
Da_folders = dir([data_dir, 'Da*']);
for Da_i = 1:length(Da_folders)
    
    thisDafolder = [data_dir, Da_folders(Da_i).name, '/'];
    
    CR_folders = dir([thisDafolder, 'CR*']);
    
    for CR_i = 1:length(CR_folders)
        thisCRfolder = [thisDafolder, CR_folders(CR_i).name, '/'];
        
        RaC_folders = dir([thisCRfolder, 'RaC*']);
        
        for RaC_i = 1:length(RaC_folders)
            thisRaCfolder = [thisCRfolder, RaC_folders(RaC_i).name, '/'];
            
            pts_folders = dir([thisRaCfolder, 'pts*']);
           
            for pts_i = 1:length(pts_folders)
               thisPtsfolder = [thisRaCfolder, pts_folders(pts_i).name, '/'];
               folders{end+1} = thisPtsfolder;
               
               fprintf('%s \n', thisPtsfolder);
            end
            
        end
        
        
    end
    
   
    
end

end

