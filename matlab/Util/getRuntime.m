function [ time ] = getRuntime( folder )
time = NaN;

timeFile = [folder, '/time.table.0'];
if exist(timeFile, 'file') == 2
    %fprintf(['Loading', file, '\n']);
    file_contents = fileread(timeFile);
    fid = fopen(timeFile, 'r');
else
    return;
end

pattern = ['main (\d+\.?\d+).*'];
                [mat,tok,ext]  = regexp(file_contents, pattern, 'match', ...
                    'tokens', 'tokenExtents');
                

temp = tok{1};          
time = str2num(temp{1});

stop = 0;
end

