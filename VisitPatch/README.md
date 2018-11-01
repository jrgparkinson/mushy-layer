# Visit Patch
Chombo support in Visit is provided via the `/src/databases/Chombo/` folder, within which the avtChomboFileFormat class describes how to open Chombo files.

I have made some small changes to this class in order to open files where a series of FArrayBoxes have been written to the file separately rather than being written as one multicomponent FArrayBox. To be clear, the standard class expects data that looks like:

```
level_0
    FArrayBox 'data':
        component 0: Some field
        component 1: Another field
        etc.
level_1
    FArrayBox 'data':
        component 0: Some field
        component 1: Another field
        etc.
```
        
Whilst my patch also let's you read data like:

```
level_0 
    FArrayBox 'some field'
    FArrayBox 'Another field'
    
level_1
    FArrayBox 'some field'
    FArrayBox 'Another field'
```

## Brief summary of patch
* Add Protected field `bool m_isCheckpointFile;`
* `avtChomboFileFormat(...)`: initialise new field to `False` by default.
* `InitializeReader()`: if the existing check determines that no data exists, assume this is a checkpoint file (set `m_isCheckpointFile=true`)
* `getVar(...)`: change the definition of `nvars` and `varIdx` if we have a checkpoint file 
```c
if (m_isCheckpointFile)
   {
	nVars = 1;
	varIdx = 0;
}
```
then get this data if the code failed to find data in the format it was expecting:
```c
if (m_isCheckpointFile)
        {
        	std:string varname_string(varname);
		std::string name = varname_string + ":datatype=0";
        	data = H5Dopen(level_id, name.c_str());
        }
```
* Do the same for `GetVectorVar(...)`. I haven't tested this yet though, so probably doesn't work.
