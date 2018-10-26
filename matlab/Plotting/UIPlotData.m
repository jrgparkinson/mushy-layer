% GUI to display the spatial structure of fields for
% the chosen (CR, Ra) pairs

function UIPlotData

%clear all;
close all;

% Initial plotfile
outputDir = pwd;
outputDir = ['/home/parkinsonjl/convection-in-sea-ice/MushyLayer/', ...
    'execSubcycle/'];


%prefix = 'VariableMesh2SubcycleRefluxFreestream0.45-convectionDB-32-ref2-';
%frame = 325;

% Initialise everything to be blank/boring
plot_prefixes = {'None'};
frames = [0];
compNames = {''};
xmin = 0; ymin = 0;
xmax = 0; ymax = 0;
ml=NaN;
mlBiLinear = NaN;
frame = 0;
prefix = plot_prefixes{1};
comp = 0;

setSliceDefaults();

% Make GUI
f = figure();
set(f, 'Position', [100 200 1600 1000]);

axLeft = axes('Parent', f, 'position', [0.1 0.2 0.35 0.65]);
axRight = axes('Parent', f, 'position', [0.6 0.2 0.35 0.65]);

topRowY = 850;

outputDirButton = uicontrol('Style','pushbutton',...
    'Callback',@chooseFolder,...
    'String', 'Folder', ...
    'Position', [50 topRowY+25 100 25]);


popupPrefixes = [175 topRowY 300 50];
popupPrefixes = uicontrol('Style', 'popup',...
    'String', plot_prefixes,...
    'Position', popupPrefixes,...
    'Callback', @setPrefix);
%popupPrefixes.Value = find(plot_prefixes==prefix);

popupFrames = [500 topRowY 100 50];
popupFrames = uicontrol('Style', 'popup',...
    'String', num2cell(frames),...
    'Position', popupFrames,...
    'Callback', @setFrame);
%popupFrames.Value = find(frames==frame);

reloadButton = uicontrol('Style','pushbutton',...
    'Callback',@reload,...
    'String', 'Load', ...
    'Position', [700 topRowY+25 100 25]);



popupComp = [950 topRowY 200 50];
popupComp = uicontrol('Style', 'popup',...
    'String', compNames,...
    'Position', popupComp,...
    'Callback', @setComp);
%popupComp.Value = comp;



 txtXSlice = uicontrol('Style','text',...
     'Position',[1200 topRowY + 50 0 0],...
     'String','x slice:'); %$\mathcal{C}$:

sldXSlice = uicontrol('Style', 'slider',...
    'Min',xmin,'Max',xmax,'Value',xslice,...
    'Position', [1200 topRowY 100 50],...
    'Callback', @changeXSlice);

txtYSlice = uicontrol('Style','text',...
     'Position',[1300 topRowY + 50 0 0],...
     'String','y slice:'); %$\mathcal{C}$:
 
sldYSlice = uicontrol('Style', 'slider',...
    'Min',ymin,'Max',ymax,'Value',yslice,...
    'Position', [1300 topRowY 100 50],...
    'Callback', @changeYSlice);



% getPlotPrefixes();
% if length(plot_prefixes) > 0
%     prefix = plot_prefixes{1};
% else
%     prefix = '';
% end
% 
% 
% getFrames();
% if length(frames) > 0
%     frame = frames(1);
% else
%     frame = NaN;
% end


  


if ~isnan(frame)
    getData();
    updatePlot();
end 




    function setSliceDefaults()
       xslice = (xmin+xmax)/2;
        yslice = (ymin+ymax)/2; 
    end


    function getPlotPrefixes()
        
        plot_prefixes = getChomboPrefixes(outputDir);
        
        set(popupPrefixes, 'String', plot_prefixes);
        
        if length(plot_prefixes) > 0
            prefix = plot_prefixes{1};
        end
        
        
        
        temp =0;
    end

    function getFrames()
        % Default option
        frames = [];
        
        if prefix
            
            files = dir([outputDir, prefix, '*']);
            for i=1:length(files)
                filename = files(i).name;
                thisFrame = regexprep(filename,'\.\dd\.hdf5','');
                thisFrame = regexprep(thisFrame,prefix,'');
                frameNum = str2num(thisFrame);
                frameNum2str = num2str(frameNum);
                if frameNum == str2num(frameNum2str)
                    frames(end+1) = str2num(thisFrame);
                end
            end
            set(popupFrames, 'String', num2cell(frames));

            fprintf('New frames: %d \n', frames);
            
       
        end
        
        if length(frames) == 0
            frames = [0];
        end
    end

    function getData()
        fprintf('Loading plotfile states... \n');
        
        ml = MushyLayerOutput(2,frame,outputDir,prefix, true, 'nearest');
        mlBiLinear = MushyLayerOutput(2,frame,outputDir,prefix, true, 'bilinear');
        
        if length(ml.levelArray) > 0
            comp = 1; %ml.components.xAdvectionvelocity; %ml.components.Temperature;
            compNames = fieldnames(ml.components);
            
            set(popupComp, 'String', compNames);
            popupComp.Value = comp;

            % Now we have data, can do stuff that depends on it
            getMinMax();
            setSliceDefaults();

        end

        
    end

    function getMinMax()
        [X,Y] = ml.grid();
        
        xmin = min(min(X));
        xmax = max(max(X));
        
        ymin = min(min(Y));
        ymax = max(max(Y));
        
        
    end

    function updatePlot()
        
        if length(ml.levelArray) == 0
            fprintf('No data to updatePlot() \n');
            return
        end
        
        cla(axLeft);
        cla(axRight);
        
        data = ml.dataForComp(comp);
        
        [X,Y] = ml.grid();
        %data = data;
        dataBiLinear = (mlBiLinear.dataForComp(comp));
        X = X.'; Y=Y.';
        x = squeeze(X(:, 1)); y = squeeze(Y(1, :));
        xLength = x(end)-x(1);
        yLength = y(end)-y(1);
         
        pcolor(X,Y,data, 'Parent', axLeft);
        title(axLeft, compNames{comp});
      %  box on;
        
        pbaspect(axLeft, [xLength, yLength, 1]);
        
        c = colorbar('peer',axLeft,'Tag','colorbar1');
        
        xlab=xlabel(axLeft, '$x$'); ylab=ylabel(axLeft, '$y$');
        
       
        hold on
        
        %plot(axLeft, [xslice xslice], [ymin ymax]);
        %plot(axLeft, [xmin xmax], [yslice yslice]);
        l = line(axLeft, [xmin xmax], [yslice yslice]);
        l.Color = [0 0 0];
        l2 = line(axLeft, [xslice xslice], [ymin ymax]);
        l2.Color = [0 0 0];
        %l.LineWidth = 10;
        
        hold off
        
        leg = {};
        
        if ~isnan(xslice)
            [slice_x, slice_data] = getSlice(x, y, dataBiLinear, 1, xslice);
            plot(slice_x, slice_data, 'Parent', axRight);
            leg{end+1} = sprintf('x=%1.3f', xslice);
        end
        
        if ~isnan(yslice)
            [slice_y, slice_data] = getSlice(x, y, dataBiLinear, 2, yslice);
            hold on;
            plot(slice_y, slice_data, 'Parent', axRight);
            hold off;
            leg{end+1} = sprintf('y=%1.3f', yslice);
        end
        
        legend(axRight, leg);
        
        xlabel('Distance');
        ylabel(compNames{comp});
        
        
        drawnow;
         
        temp = 0;
        
        
        
        
    end


    function [slice_x, slice_data] = getSlice(x,y, data, dir, loc)
        
        if dir == 1
            xNearest = min(abs(loc-x));
            XI = find(abs(loc-x)==xNearest);
            XI = min(XI);
            slice_data = squeeze(data(XI, :));
            slice_x = x;
        elseif dir == 2
            yNearest = min(abs(loc-y));
            YI = find( min(abs(loc-y))==yNearest);
            YI = min(YI);
            slice_data = squeeze(data(:, YI));
            slice_x = y;
        end
    end


    function setComp(source, event)
        comp = source.Value;
        fprintf('Component: %d \n', comp);
        
        updatePlot();
    end

    function chooseFolder(hObject,eventdata,handles)
        outputDir = uigetdir(outputDir,'Select output directory');
        outputDir = [outputDir, '/'];
        getPlotPrefixes();
        getFrames();
        drawnow;
    end

    function reload(hObject,eventdata,handles)
        getData();
        updatePlot();
    end

    function setFrame(source, event)
        newFrame = source.Value;
        fprintf('Frame: %d \n', frames(newFrame));
        frame = frames(newFrame);
    end

 function setPrefix(source, event)
     
        prefix = plot_prefixes{source.Value};
        
        fprintf('Prefix: %s \n', prefix);
        
        getFrames();
    end

    function changeXSlice(source, event)
        xslice = source.Value;
        fprintf('x slice: %1.3f \n', xslice);
        updatePlot();
    end

    function changeYSlice(source, event)
        yslice = source.Value;
        fprintf('y slice: %1.3f \n', yslice);
        updatePlot();
    end


% sldMaxPorosity = uicontrol('Style', 'slider',...
%     'Min',0,'Max',0.97,'Value',maxPorosity,...
%     'Position', txtMaxPorosity.Position + [0 -30 -20 0],...
%     'Callback', @changeMaxPorosity);
%
% txtInfo = uicontrol('Style','text',...
%     'Position',sldMaxPorosity.Position + [0 -70 0 0],...
%     'String',''); %$\mathcal{C}$:
%
%

end