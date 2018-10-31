%function ChomboRefine(inputFilename, outputFilename, refineFactor )
function ChomboRefine(inputFilename, outputFilename, refineFactor)
%CHOMBOREFINE function to take a single level Chombo checkpoint file and refine it by the
% specified factor, saving as a new file which we can restart from
initialDir = pwd; 

if nargin == 0
    FilterSpec = {'*.hdf5', 'HDF5 Files (*.hdf5)'};
    DialogTitle = 'Select Chombo checkpoint file';
    [inputFilename, inputPath,FilterIndex] = uigetfile(FilterSpec,DialogTitle);
    
    cd(inputPath);
    
    guessedInputFilename = inputFilename; %TODO: something more clever here

    %outputFilename = input('Enter output filename: ', 's');
    [outputFilename,outputPathName] = uiputfile(FilterSpec,'Save refined file as...',guessedInputFilename);
    refineFactor = input('Enter refinement factor (1,2,4,8...): ');
end

% % Check inputs
numRefinements = log(refineFactor)/log(2); %As log_2(x) = log_y(x)/log_y(2)
if floor(numRefinements) ~= numRefinements
    error('Refinement must be a power of 2: 1, 2, 4, 8... etc. ');
end

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

% We don't really need this - we already have the boxes from earlier
%boxes = h5read(outputFilename,['/level_0/boxes']); 

% Don't think we want to change the processors, but get this anyway
procs = h5read(outputFilename,['/level_0/Processors']);

% % Refinement
% We have everything now, so refine stuff
newDx =  oldDx/refineFactor;
newProcs = procs;

newProbDomain = oldProbDomain;
c = class(newProbDomain.hi_i);

for i=1:numRefinements
    newProbDomain.hi_i = cast((newProbDomain.hi_i*2)+1, c);
    newProbDomain.hi_j = cast((newProbDomain.hi_j*2)+1, c);
end

% Refine data
for i=1:length(components)
    thisComp = components(i);
    oldData = components(i).Data;
    
    if strcmp(thisComp.Type, 'FArrayBox')
         newData = NaN*ones(Nx*Ny*refineFactor, 1);
  
        for j=1:components(i).NumComps
            dataLength = Nx*Ny;
            offset = max(1, 1+dataLength*(j-1));

            thisDirData = oldData(offset:(offset+dataLength-1));
            reshaped = reshape(thisDirData, Nx, Ny);

            refined = resizem(reshaped, refineFactor, interpolationMethod);

            
            refinedDataLength = double(dataLength)*(refineFactor^double(dim));
            refinedOffset = max(1, 1+refinedDataLength*(j-1));

            unshaped = reshape(refined, refinedDataLength, 1);
            newData(refinedOffset:(refinedOffset+refinedDataLength-1)) = unshaped;

        end
       
        components(i).Data = newData;
        
    elseif strcmp(thisComp.Type, 'unknown') && thisComp.NumComps == 1
        %Assume this is a fluxbox.
        
        offset = 1;
        refinedOffset = 1;
        
        for j = 1:dim
            % Flux box has an extra component in this direction
            thisNx = Nx; thisNy = Ny;
            if j == 1
               thisNx = thisNx + 1;
            else
               thisNy = thisNy + 1;
            end
            
            dataLength = thisNx*thisNy;
            refinedDataLength = double(dataLength)*(refineFactor^double(dim));
            
            thisDirData = oldData(offset:(offset+dataLength-1));
            reshaped = reshape(thisDirData, thisNx, thisNy);

            refined = resizem(reshaped, refineFactor, interpolationMethod);

            unshaped = reshape(refined, refinedDataLength, 1);
            newData(refinedOffset:(refinedOffset+refinedDataLength-1)) = unshaped;
            
            offset = offset + dataLength; % Get this ready for the next direction
            refinedOffset = refinedOffset + refinedDataLength;
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

cd(initialDir);

end