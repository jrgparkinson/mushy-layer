function makePlotDataFolder(output_dir, plot_prefix)
frames = [];
times = [];
Si = [];
Vi = [];
Sbox = [];
dim = 2; subcycled = true;

lockFile = fullfile(output_dir, '/lock');

fid = fopen( lockFile, 'wt' );
fprintf( fid, 'Currently processing this directory');
fclose(fid);

diagFile = fullfile(output_dir, '/diags.mat');
diagFileExists = (exist(diagFile, 'file')==2);

% For now skip folders we're currently processing/processed
if diagFileExists
    load(diagFile);
end

% Get all the fames
% Default option:
%parseFrames = 0:100:500000;

if strcmp(plot_prefix, '')
   plot_prefixes = getChomboPrefixes(output_dir, false);
   plot_prefix = plot_prefix{1};
end

% Cleverer option:
parseFrames = [];
searchTerm = fullfile(output_dir, [plot_prefix, '*']);
d = dir(searchTerm);
for d_i = 1:length(d)
    thisFilename = d(d_i).name;
    
    thisFrame = str2num(thisFilename(end-13:end-8));
    fprintf('Parsing %s, frame = %d \n', thisFilename, thisFrame);
    parseFrames(end+1) = thisFrame;
end

%for frame = 100:100:110000
for frame_i = 1:length(parseFrames)
    frame = parseFrames(frame_i);
    
    plotFile = [output_dir, '/toPlot',num2str(frame),'.mat'];
    
    fprintf('Considering frame %d \n', frame);
    
    %if exist(plotFile, 'file') == 2
    %   continue
    %end
    
    % Skip frames already done
    if sum(frames==frame) > 0
        fprintf('   Skipping \n');
        continue
        
    end
    
    
    output = MushyLayerOutput(dim, frame, output_dir, plot_prefix, subcycled);
    %output = MushyLayerOutput(dim, frame, output_dir, plot_prefix, subcycled);
    if length(output.levelArray) > 0
        porosity = output.dataForComp(output.components.Porosity);
        Sl = output.dataForComp(output.components.Liquidconcentration);
        S = output.dataForComp(output.components.Bulkconcentration);
        T = output.dataForComp(output.components.Temperature);
        [X,Z] = output.grid();
        t = output.t;
        save(plotFile, 'porosity', 'Sl', 'T', 'X','Z','t');
        
        sprintf('%s \n', plotFile);
        
        
        % Add diagnostics
        times(end+1) = t;
        
        %iceCells = porosity<1;
        
        Si(end+1) = mean(mean(S(porosity<1)));
        Vi(end+1) = sum(sum(porosity<1))/(size(porosity,1)*size(porosity,2));
        Sbox(end+1) = mean(mean(S));
        
        frames(end+1) = frame;
        
        fprintf('Si=%1.3f, Vi=%1.2f, Sbox=%1.2f\n', Si(end), Vi(end), Sbox(end))
        
        % save at intermediate times too
        if length(times) > 0
            save(diagFile, 'times', 'Si', 'Vi', 'Sbox', 'frames');
        end
        
        
    end
    
    
end % end loop over frames

if length(times) > 0
    fprintf('Saving final diags file\n');
    save(diagFile, 'times', 'Si', 'Vi', 'Sbox', 'frames');
end

% Free up directory
fprintf('Deleting lock file\n');
delete(lockFile);
end