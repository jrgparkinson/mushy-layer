function ChomboExtend(inputFilename, outputFilename, extrapCells )
%CHOMBOREFINE function to take a single level Chombo checkpoint file,
%periodic in the x-direction, and add some extra extrapolated cells on each
% side, changing the aspect ratio

% % Check inputs
extrapCells = int32(extrapCells);

% % Setup constants
LEVEL0 = '/level_0/';
DX_STR = 'dx';
PROBDOMAIN_STR = 'prob_domain';
BILINEAR = 'bilinear'; BICUBIC = 'bicubic'; NEAREST = 'nearest';
interpolationMethod = BICUBIC; 


%Firstly, make a copy of the original file with the new name
copyfile(inputFilename, outputFilename);

%if numRefinements == 0
%    % Don't need to do anything else
%   return; 
%end

% % Load data
% Now get all the data that we need to edit
dim = h5readatt(outputFilename, '/Chombo_global/', 'SpaceDim');
oldDx = h5readatt(outputFilename, LEVEL0, DX_STR);
oldProbDomain = h5readatt(outputFilename, LEVEL0, PROBDOMAIN_STR);
num_levels = h5readatt(outputFilename, '/', 'num_levels');
%num_comps =  h5readatt(outputFilename, '/', 'num_components');

Nx = oldProbDomain.hi_i + 1 - oldProbDomain.lo_i;
Ny = oldProbDomain.hi_j + 1 - oldProbDomain.lo_j;

if (num_levels > 1)
   error('Currently only able to handle single level files'); 
end



%Then edit data
% First load the data...
% 1) Get a list of data types
info = hdf5info(outputFilename);

comp_names = {};
comp_numComps = {};
objTypes = {};

groups = info.GroupHierarchy.Groups(2).Groups;

for i = 1:length(groups)
   name = groups(i).Name;
   
   %Strip off extra stuff
   name = strrep(name, '/level_0/', '');
   name = strrep(name, '_attributes', '');
   comp_names{end+1} = name;
   
   %thisComp.Name = name;
   
   for at_i = 1:length(groups(i).Attributes)
       att_name = groups(i).Attributes(at_i).Shortname;
       att_val = groups(i).Attributes(at_i).Value;
       
       if strcmp(att_name, 'objectType')
           %type = att_val;
           objTypes{end+1} = att_val.Data;
       elseif strcmp(att_name, 'comps')
            comp_numComps{end+1} = att_val;
        end
   end
   %objectType = groups(i).Attributes(3).Value;
   
  
   
%   components(end+1) = thisComp;
    
end

data = NaN*ones(length(comp_names), 1);

components = struct('Name', comp_names, 'NumComps', comp_numComps, 'Type', objTypes, 'Data', data);


for i = 1:length(comp_names)
    thisComp = components(i);
   
    components(i).Data = h5read(outputFilename,['/level_0/',thisComp.Name,':datatype=0']);
   
end

% Don't think we want to change the processors, but get this anyway
procs = h5read(outputFilename,['/level_0/Processors']);

% % Extrapolation
% We have everything now, so extrapolate all fields
newDx =  oldDx;
newProcs = procs;

newProbDomain = oldProbDomain;
c = class(newProbDomain.hi_i);

% Extend domain in x-direction
newProbDomain.hi_i = cast((newProbDomain.hi_i+2*extrapCells), c);

newNx = Nx + 2*extrapCells;

% Shift and extrapolate data
for i=1:length(components)
    thisComp = components(i);
    oldData = components(i).Data;
    
    if strcmp(thisComp.Type, 'FArrayBox')
         newData = NaN*ones(newNx*Ny, 1);
  
        for j=1:components(i).NumComps
            dataLength = Nx*Ny;
            offset = max(1, 1+dataLength*(j-1));

            thisDirData = oldData(offset:(offset+dataLength-1));
            reshaped = reshape(thisDirData, Nx, Ny);

           % refined = resizem(reshaped, refineFactor, interpolationMethod);
           extended = NaN*ones(newNx, Ny);
           extended(extrapCells+1:end-extrapCells, :) = reshaped;
           for extrap_i=1:(extrapCells)
                extended(extrap_i, :) = reshaped(1, :);
                extended(end+1-extrap_i, :) = reshaped(end, :);
           end
          

            extendedDataLength = size(extended, 1)*size(extended,2);
            extendedOffset = max(1, 1+extendedDataLength*(j-1));

            unshaped = reshape(extended, extendedDataLength, 1);
            newData(extendedOffset:(extendedOffset+extendedDataLength-1)) = unshaped;

        end
       
        components(i).Data = newData;
        
    elseif strcmp(thisComp.Type, 'unknown') && thisComp.NumComps == 1
        %Assume this is a fluxbox.
        
        offset = 1;
        extendedOffset = 1;
        
        for j = 1:dim
            % Flux box has an extra component in this direction
            thisNx = Nx; thisNy = Ny;
            if j == 1
               thisNx = thisNx + 1;
            else
               thisNy = thisNy + 1;
            end
            
            dataLength = thisNx*thisNy;
            extendedDataLength = (thisNx+2*extrapCells)*thisNy;
            
            thisDirData = oldData(offset:(offset+dataLength-1));
            reshaped = reshape(thisDirData, thisNx, thisNy);

            %refined = resizem(reshaped, refineFactor, interpolationMethod);
            extended = NaN*ones(thisNx+2*extrapCells, thisNy);
            extended(extrapCells+1:end-extrapCells, :) = reshaped;
           for extrap_i=1:(extrapCells)
                extended(extrap_i, :) = reshaped(1, :);
                extended(end+1-extrap_i, :) = reshaped(end, :);
           end

            unshaped = reshape(extended, extendedDataLength, 1);
            newData(extendedOffset:(extendedOffset+extendedDataLength-1)) = unshaped;
            
            offset = offset + dataLength; % Get this ready for the next direction
            extendedOffset = extendedOffset + extendedDataLength;
        end
        
        components(i).Data = newData;
       
        
    else
        error('Cannot handle data type');
    end
    
    
end
%Before we can refine data, need to reshape it

% % Replace old data
% Finally, write it back again

% New Dx
h5writeatt(outputFilename, LEVEL0, DX_STR, newDx)

% New problem domain
h5writeattCompound(outputFilename, LEVEL0, PROBDOMAIN_STR, newProbDomain);
h5writeCompound(outputFilename, '/level_0/boxes', newProbDomain);

fid = H5F.open(outputFilename, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
gid  = H5G.open(fid, LEVEL0);

% Refined Data
for i = 1:length(components)
   
    thisComp = components(i);
    dataToWrite = thisComp.Data;
    
    % Delete old data, insert new data
   
    H5L.delete(gid,[thisComp.Name,':datatype=0'], 'H5P_DEFAULT');
    
    h5create(outputFilename,['/level_0/',thisComp.Name,':datatype=0'],length(dataToWrite));
    h5write(outputFilename, ['/level_0/',thisComp.Name,':datatype=0'], dataToWrite);
    
  
    h5write(outputFilename, ['/level_0/',thisComp.Name,':offsets=0'], [0 int64(length(dataToWrite))]);

end

 % New processors
 h5write(outputFilename,['/level_0/Processors'], newProcs);

%H5A.close(attr_id);
H5F.close(fid);

end