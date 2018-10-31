function h5writeCompound(outputFilename, dsetLocation, newProbDomain)

% Need to extract group and dataset from dsetLocation
C = strsplit(dsetLocation, '/');
dataset = C{end};
group = '';
for i=1:(length(C)-1)
    group = strcat(group, '/', C{i});
end

memtype = h5CompoundMemtype(newProbDomain);

fid = H5F.open(outputFilename, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
gid  = H5G.open(fid, group);
dset = H5D.open(gid, dataset);
H5D.write (dset, memtype, 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', newProbDomain);

end