function diags = getDiagnostics(file)
diags = struct();

if exist(file, 'file') == 2
    %fprintf(['Loading', file, '\n']);
    file_contents = fileread(file);
    fid = fopen(file, 'r');
else
    return;
end


count = 1;
while isempty(file_contents) && count < 10
   % fprintf('Waiting for diagnostics file to not be empty');
    pause(1.0);
    %fclose(fid);
    file_contents = fileread(file);
    %fid = fopen(file, 'r');
    count = count + 1;
end

if count == 10
    return;
end

% In some odd cases we may have multiple header lines
%numHeaderLines = count(file_contents, 'time'); % matlab 2017b+
% Try both capitalised and not:
match1 = strfind(file_contents, 'time','ForceCellOutput',true);
match2 = strfind(file_contents, 'Time','ForceCellOutput',true);
numHeaderLines = max(length(match1{1}), length(match2{1}));


% Get everything up to the first end line
header = textscan(fid, '%[^\n]', 1);

% Count total number of lines
n = 0;
tline = fgetl(fid);
while ischar(tline)
    tline = fgetl(fid);
    n = n+1;
end
fclose(fid);
numDataLines = n - numHeaderLines;

% Return this to original state
%fid2 = fopen(file, 'r');

%fclose(fid);

% header is a 1x1 cell so need to get first entry
h = header{1};

if length(h) < 1
    err = 0;
end

%Collapse delimiters = false so that we still get blank entries
bits = strsplit(h{1}, ',', 'CollapseDelimiters', false);

if numDataLines > 0
    data = csvread(file, numHeaderLines, 0);
else
    data = NaN*ones(1, length(bits));
    
end

for i = 1:length(bits)
    % Ignore data sets where the header is blank
    try
        if ~ (strcmp(bits{i}, ''))
            bits{i} = strrep(bits{i}, ' ', '_');
            diags.(bits{i}) = squeeze(data(:, i));
        end
    catch e
       fprintf('Could not add field to structure, identifier:\n%s \n',e.identifier);
       fprintf(1,'Error message:\n%s \n',e.message); 
    end
end



end