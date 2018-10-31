function memtype =  h5CompoundMemtype(structure)

%Now we need to work out how to write val to HDF5
names = fieldnames(structure);
%types = []; 
sz = NaN*ones(length(names), 1);

 % Need to make sure this is an array of H5ML.id objects
types = H5ML.id.empty(length(names), 0);

for i = 1:length(names)
    dType = class(structure.(names{i}));
    
    % Need to convert matlab class to H5 type. Not a full list.
    intType   = H5T.copy('H5T_NATIVE_INT');
    doubleType = H5T.copy('H5T_NATIVE_DOUBLE');
    H5type = intType; % default option
    
    if strcmp(dType, 'int')
        H5type = intType;
    elseif strcmp(dType, 'double')
        H5type = doubleType;
    end
    
    sz(i) = H5T.get_size(H5type);
    types(end+1) = H5type;
    
end


intType   = H5T.copy('H5T_NATIVE_INT');
%sz(1)     =H5T.get_size(intType);
%sz(2)     =H5T.get_size(intType);
%sz(3)     =H5T.get_size(intType);
%sz(4)     =H5T.get_size(intType);

offset(1)=0;
offset(2:4)=cumsum(sz(1:3));

memtype = H5T.create ('H5T_COMPOUND', sum(sz));

for i = 1:length(names)
    H5T.insert (memtype, names{i},offset(i),types(i));
end
%H5T.insert (memtype, 'lo_i',offset(1),intType);
%H5T.insert (memtype, 'lo_j',offset(2), intType);
%H5T.insert (memtype,  'hi_i',offset(3), intType);
%H5T.insert (memtype,  'hi_j',offset(4), intType);
end