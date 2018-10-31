function h5writeattCompound(file, loc, att, structure)


fid = H5F.open(file, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
gid  = H5G.open(fid, loc);


memtype = h5CompoundMemtype(structure);

attr_id = H5A.open(gid, att);
H5A.write (attr_id, memtype, structure);

end