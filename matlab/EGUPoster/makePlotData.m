function makePlotData(subfolders, base_dir,  plot_prefix)

fprintf('Loaded makePlotData\n');

if nargin < 3
    plot_prefix = 'plt';
end

if nargin < 2
    base_dir = getDataDir('middleton/');
end

if nargin < 3
    
   % base_dir =  getDataDir('middleton/'); %'/media/parkinsonjl/FREECOM HDD/';

    ending3 = '0';
    ending2 = 'restart2';
    ending1 = 'restart';

    subfolders = {['CR1.179RaC80Le200KozenyPermeabilitypts256-S3.5-TB-10.0-R0.13-',ending1,'/'], ...
        [ 'CR1.179RaC80Le200KozenyPermeabilitypts256-S3.5-TB-15.0-R0.13-',ending1,'/'], ...
        [ 'CR1.179RaC80Le200KozenyPermeabilitypts256-S3.5-TB-20.0-R0.13-',ending1,'/'], ...
        [ 'CR1.179RaC80Le200KozenyPermeabilitypts256-S3.5-TB-10.0-R0.13-',ending2,'/'], ...
        [ 'CR1.179RaC80Le200KozenyPermeabilitypts256-S3.5-TB-15.0-R0.13-',ending2,'/'], ...
        [ 'CR1.179RaC80Le200KozenyPermeabilitypts256-S3.5-TB-20.0-R0.13-',ending2,'/'], ...
        [ 'CR1.179RaC80Le200KozenyPermeabilitypts256-S3.5-TB-10.0-R0.13-',ending3,'/'], ...
        [ 'CR1.179RaC80Le200KozenyPermeabilitypts256-S3.5-TB-15.0-R0.13-',ending3,'/'], ...
        [ 'CR1.179RaC80Le200KozenyPermeabilitypts256-S3.5-TB-20.0-R0.13-',ending3,'/'], ...
        [ 'T-10Outflow/'], ...
        [ 'T-15Outflow/'], ...
        [ 'T-20Outflow/']         };

end

for i=1:length(subfolders)
    full_folders{i} = [base_dir, subfolders{i}];
end

dim = 2; subcycled = true;

% Three folders to process
for axi = 1:length(full_folders)

    % frame = frames(f_i);
    %frame = frames(axi);
    output_dir = full_folders{axi};
    
    fprintf('Output dir: %s\n', output_dir);
    
    frames = [];
    times = [];
    Si = [];
    Vi = [];
    Sbox = [];
    
    diagFile = [output_dir, '/diags.mat'];
     
    diagFileExists = (exist(diagFile, 'file')==2);
    
    % If diag file exists, check it hasn't been edited recently (this is the sign of a folder already being processed)
    lockFile = [output_dir, '/lock'];
    lockFileExists = (exist(lockFile, 'file')==2);
    
    if lockFileExists
        fprintf('Directory locked, skipping\n');
        continue
    end
    
    
    try 
    
        fid = fopen( lockFile, 'wt' );
        fprintf( fid, 'Currently processing this directory');
        fclose(fid);
            
            
        
        
        % For now skip folders we're currently processing/processed
        if diagFileExists
            load(diagFile);
        end
        
        % Get all the fames
        % Default option:
        %parseFrames = 0:100:500000;
        
        % Cleverer option:
        parseFrames = [];
        d = dir([output_dir,'plt*']);
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
            [X,Z] = output.grid();
            t = output.t;
            save(plotFile, 'porosity', 'Sl','X','Z','t');
            
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
        
    catch
        fprintf('Error, skipping\n');
    
    end
    
    
    
end % end loop over folders

end