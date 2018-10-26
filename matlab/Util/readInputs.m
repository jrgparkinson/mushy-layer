function inputs = readInputs(file)
fileID = fopen(file,'r');

inputs = struct;

if exist(file, 'file') == 2
    %fprintf(['Loading', file, '\n']);
    file_contents = fileread(file);
else
    return;
end


% Parse file_contents
regex = '\.(.[^\n\.]*)=(.[^\n]*)[\n]';
 [~,matches,~]  = regexp(file_contents, regex, 'match', 'tokens', 'tokenExtents');
                
               
for i=1:length(matches)
   key = matches{i}{1};
   val = matches{i}{2};
   
   inputs.(key) = val;
    
end




end