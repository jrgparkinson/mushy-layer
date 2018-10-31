classdef ChomboOutput < handle
    properties
        t
        amr_data
        levelArray
        components
        problemDomain
        output_dir
        output_prefix
        frame
        filename
    end
    
    methods
        function obj = ChomboOutput(dim,Frame,dir,outputprefix, interpolationMethod)
            if nargin < 5
                interpolationMethod = 'nearest';
            end
             INFINITE = 1e300;
            
           fname = ChomboOutput.getFilename(dir, outputprefix, dim, Frame);
           
           obj.output_dir = dir;
           obj.output_prefix = outputprefix;
           obj.frame = Frame;
           obj.filename = fname;
           
           % If we still don't have a file, let's quit
            if ~exist(fname)
                amr = {};
                obj.t = [];
                disp(' ');
                disp(['Frame ',num2str(Frame),' (',fname,') does not exist ***']);
                disp(' ');
                return;
            end 
            
            %disp(['Reading data from ',fname]);
            
            %get time from the 0th level
            obj.t = h5readatt(fname,'/level_0', 'time');
            
            num_comps = double(h5readatt(fname, '/', 'num_components'));
            num_levels = h5readatt(fname, '/', 'num_levels');
            
            %read in the component names
            for comp_i = 0:(num_comps-1)
                comp_name = strcat('component_', num2str(comp_i));
                strName = h5readatt(fname, '/', comp_name);
                %name = cellstr(strName);
                %comp_names(comp_i+1) = name;
                
                structName = strrep(strName, ' ', '');
                obj.components.(structName) = comp_i+1;
            end
            %obj.components = comp_names;
           
            
            % These are the low values for the level 0 grid.
            xlower = 0;
            ylower = 0;
            zlower = 0;
            
            domainExtent = h5readatt(fname, '/level_0', 'prob_domain');

            fileinfo = hdf5info(fname);
            
            % This is the maximum number of levels created, and may be smaller
            max_level = length(fileinfo.GroupHierarchy.Groups) - 2;
            num_levels = max_level+1;
            
            refRatios = ones(num_levels, 1);
            
            clear amr;
            ng = 1;
            for level = 0:max_level,
                
                % read parameters for this level
                gstring = ['/level_',num2str(level)];
                
                ghost_vect = h5readatt(fname,'/level_0/data_attributes','outputGhost');
                num_ghost_x = ghost_vect.intvecti;
                num_ghost_y = ghost_vect.intvectj;
                
                if dim == 3
                    num_ghost_z = ghost_vect.intvectk;
                end
                
                %dx = hdf5read(fname,[gstring,'/dx']);
                dx = h5readatt(fname, gstring, 'dx');
                dy = dx;
                if (dim == 3)
                    dz = dx;
                end;
                
                if level==0
                    dxCoarse = dx;
                end
                
                % Dimensions of boxes at this level.
                % This is where we really read everything in.
                
                data = hdf5read(fname,[gstring,'/data:datatype=0']);
                data(data > INFINITE) = nan;
                didx = 0;
                
                refRatios(level+1) = h5readatt(fname,'/level_0/','ref_ratio');
                
                boxes = hdf5read(fname,[gstring,'/boxes']);
                for j = 1:length(boxes),
                    clear amrdata;
                    amrdata.gridno = ng;
                    amrdata.level = level + 1;  % Chombo levels start at 0;
                    amrdata.dx = dx;
                    amrdata.dy = dx;
                    if (dim == 3)
                        amrdata.dz = dx;
                    end;
                    
                    thisBox = boxes(j);
                    
                    if (dim == 2)
                        
                        lo_i = double(boxes(j).Data{1}) - num_ghost_x;
                        lo_j = double(boxes(j).Data{2}) - num_ghost_y;
                        hi_i = double(boxes(j).Data{3}) + num_ghost_x;
                        hi_j = double(boxes(j).Data{4}) + num_ghost_y;
           
                    elseif (dim == 3)
                        lo_i = double(boxes(j).Data{1});
                        lo_j = double(boxes(j).Data{2});
                        lo_k = double(boxes(j).Data{3});
                        hi_i = double(boxes(j).Data{4});
                        hi_j = double(boxes(j).Data{5});
                        hi_k = double(boxes(j).Data{6});
                    end;
                    amrdata.mx = hi_i - lo_i + 1;
                    amrdata.my = hi_j - lo_j + 1;
                    if (amrdata.mx == 1 || amrdata.my == 1)
                        % Can't plot these tiny grids...
                        % continue;
                    end;
                    if (dim == 3)
                        amrdata.mz = hi_k - lo_k + 1;
                        if (amrdata.mz == 1)
                            continue;
                        end
                    end;
                    
                    amrdata.xlow = xlower + dx*double(lo_i);
                    amrdata.ylow = ylower + dy*double(lo_j);
                    if (dim == 3)
                        amrdata.zlow = zlower + dz*lo_k;
                    end;
                    
                    len = amrdata.mx*amrdata.my;
                    if (dim == 3)
                        len = len*amrdata.mz;
                    end;

                    %The data for this box
                    box_data = data(didx+1:didx+len*num_comps);
                    
                    amrdata.data = reshape(box_data, len, num_comps)';
                    didx = didx + len*num_comps;
                    
                    amr(ng) = amrdata;
                    ng = ng + 1;
                end;  % end loop over boxes on this level
            end;  % end loop on levels
            
            %Create an object to store problem domain information

            obj.problemDomain = ProblemDomain(num_comps, num_levels, max_level, dim, num_ghost_x, num_ghost_y, interpolationMethod, domainExtent, dxCoarse, refRatios);

  
            %Now we have the AMR data, sort it into objects
            for level=1:(max_level+1)
                box_i = 1;
                
                for ngrid = 1:(ng-1)
                    if amr(ngrid).level == level
                       boxArray(box_i) = amr(ngrid); 
                       box_i = box_i + 1;
                    end
                end
                
                chomboLevel = ChomboLevel(level, boxArray, obj.problemDomain);
                
                levels(level) = chomboLevel;
                clear boxArray;
            end
            
            obj.levelArray = levels;
            
        end
        
        function amrPolyshapes = getMeshes(obj, domainBox)
            if nargin < 2
                dom = obj.problemDomain;
                levDx = dom.dxCoarse;
                xPoints = [dom.xlo-levDx dom.xlo-levDx dom.xhi+levDx dom.xhi+levDx];
                yPoints = [dom.ylo-levDx dom.yhi+levDx dom.yhi+levDx dom.ylo-levDx];
               
                domainBox = polyshape(xPoints, yPoints);
            end
            
            amrPolyshapes = [polyshape];
            
            for l = 2:length(obj.levelArray)
                levelPolyshapes = obj.levelArray(l).getMeshes(domainBox);
                
                amrPolyshapes(l) = levelPolyshapes;
            end 
        end
        
        % Return grid covering the whole domain at the finest level
        % (including ghost cells if specified)
        function [X,Y] = grid(obj, ghostCells)
            if nargin < 2
                ghostCells = false;
            end
            
            % Extent of coarsest grid without ghost cells
            [x_min, y_min, x_max, y_max] = obj.levelArray(1).grid_extent_no_ghost();

            %Enlarge slightly for finest grid
            dx_coarse = obj.levelArray(1).dx;
            dx_fine = obj.levelArray(end).dx;
            
            x_min = x_min - dx_coarse/2 + dx_fine/2;
            y_min = y_min - dx_coarse/2 + dx_fine/2;
            x_max = x_max + dx_coarse/2 - dx_fine/2;
            y_max = y_max + dx_coarse/2 - dx_fine/2;
            
            if ghostCells
                %Add ghost cells back in
                ghostx = double(obj.problemDomain.num_ghost_x);
                ghosty = double(obj.problemDomain.num_ghost_y);

                x_min = x_min - ghostx*dx_fine; y_min = y_min - ghosty*dx_fine;
                x_max = x_max + ghostx*dx_fine; y_max = y_max + ghosty*dx_fine;
            end
                
            %Generate fine grid
            [X,Y] = meshgrid(x_min:dx_fine:x_max, y_min:dx_fine:y_max);
        end
        
        % Return grid covering the whole domain at the specified level
        function [X,Y] = domainGrid(obj, level)
            
            if nargin < 2
                level = 1;
            end
         
            % Extent of coarsest grid without ghost cells
            [x_min, y_min, x_max, y_max] = obj.levelArray(1).grid_extent_no_ghost();

            %Enlarge slightly for finest grid
            dx_coarse = obj.levelArray(1).dx;
            dx_level = obj.levelArray(level).dx;
            
            x_min = x_min - dx_coarse/2 + dx_level/2;
            y_min = y_min - dx_coarse/2 + dx_level/2;
            x_max = x_max + dx_coarse/2 - dx_level/2;
            y_max = y_max + dx_coarse/2 - dx_level/2;

            %Generate level grid
            [X,Y] = meshgrid(x_min:dx_level:x_max, y_min:dx_level:y_max);
        end
        
        % Combine the data on all levels
        % global refinement refines all levels
        function compositeGrid = data(obj, ghostCells, global_refinement)
            if nargin < 3
                global_refinement = 1; % default is no global refinement
            end
            %[X, Y] = obj.grid();
 
            dx_fine = obj.finest_dx();
            
            compositeGrid = [];
            
            for lev = 1:length(obj.levelArray)
                chomboLevel = obj.levelArray(lev);
                
                % Get this level's data at the finest level
                refinement = round(chomboLevel.dx/dx_fine) * global_refinement;             
                
                %Add this data to the combined grid
                if lev==1
                    %For the coarsest level, just take all the data
                   levelData = chomboLevel.getData(ghostCells, refinement);
                   compositeGrid = levelData; 
                else
                    %For every other level, just get the boxes of data on
                    %the level
                    
                    for box_i = 1:length(chomboLevel.boxArray)
                        if ghostCells
                           printf('warning: chomboOutput.compositeGrid cannot handle ghostCells yet on level > 0'); 
                        end
                        
                        box = chomboLevel.boxArray(box_i);
                        box_data = box.data_no_ghost(refinement);
                        [lo_i, lo_j, hi_i, hi_j] = box.get_extent_no_ghost(refinement);
                    
                        temp = compositeGrid(:, lo_i : hi_i, lo_j : hi_j);
                        
                        %box_data = nan*box_data;
                        
                        compositeGrid(:, lo_i : hi_i, lo_j : hi_j) = box_data;
                        
                        
                    end % end loop over boxes
                end
            end   % end loop over levels
            
        end
        
        function num_comps = num_comps(obj)
            num_comps = length(obj.components);
        end
        
        function dataForComp = dataForComp(obj, comp, includeGhost)   
            if nargin < 3
                includeGhost = false;
            end
            
            % return a blank array if there isn't any data
            if length(obj.levelArray) < 1
                dataForComp = NaN*ones(3);
                return;
            end
            
            data = obj.data(includeGhost);    
            
            dataComp = data(comp,:,:);
            dataForComp = squeeze(dataComp);
        end
        
        function dx = finest_dx(obj)
             dx = obj.levelArray(end).dx;
        end
        
        function c = plot(obj, var)
%             figure();
%             data = obj.dataForComp(var);
%             [X, Y] = obj.grid();
%             [c,h] = contourf(X, Y, data, 100);
%             set(h, 'linestyle', 'none'); %remove contour lines
%             colorbar();
            if length(obj.levelArray()) > 0
                    [X, Y] = obj.grid();
                    flatOutput = ChomboOutputFlat(obj.data(false), X, Y, obj.finest_dx());
                    c = flatOutput.plot(var);
            end
        end
        
        
        % return the difference between this and another output, 
        % in the form of a mushy layer output defined on this grid
        function diff = subtract(obj, other_output)
            includeGhost = false;
            [X, Y] = obj.grid();
            this_data = obj.data(includeGhost);
            other_data = other_output.data(includeGhost);
            
            % Check data sets have same number of components
            this_num_comps = size(this_data, 1);
            other_num_comps = size(other_data, 1);
            if this_num_comps ~= other_num_comps
                diff = ChomboOutputFlat(this_data.*NaN, X, Y, obj.finest_dx());
                return;
            end
            
            % Try rescaling the coarsest grid with nearest neighbour
            % interpolation if the grids aren't the same size
            s_this = size(squeeze(this_data(1, :, :))); s_other = size(squeeze(other_data(1, :, :)));
            if (s_this ~= s_other)
               
                refinement = obj.finest_dx() / other_output.finest_dx();
                
                % refinement > 1 means obj is coarser
                if (refinement > 1)
                    this_data = obj.data(includeGhost, round(refinement));
                else
                    % refinement < 1 means other_output is coarser
                    other_data = other_output.data(includeGhost, round(1/refinement));
                end   
            end
            
            diff_data = this_data - other_data;
            diff = ChomboOutputFlat(diff_data, X, Y, obj.finest_dx());
        end
        
        
        
        
        % check if we have some data
        function exist = exists(obj)
            if length(obj.levelArray) < 1            
                exist = false;
            else 
                exist = true;
            end
        end

        
    end
    
    methods(Static)
        function fname = getFilename(dir, outputprefix, dim, Frame)
            %Check dir exists!
            if 7~=exist(dir,'dir')
                disp('output directory does not exists!');
            end
            
           
            %dim = h5readatt(filename, '/Chombo_global', 'SpaceDim');
            
            
            % Try to find this file. If a frame is specified, try different
            % formats (0001, 00001, 000001 etc.)
            fname = '';
            
            offset = 10;
            while offset < 1e10 && ~exist(fname)
                if Frame > -1
                    nstr = num2str(Frame+offset);
                    nstr = nstr(2:end);
                    fname = fullfile(dir,[outputprefix, nstr,'.',num2str(dim),'d.hdf5']);
                else
                    fname = fullfile(dir, [outputprefix,'.',num2str(dim),'d.hdf5']);
                end
            
                offset = offset * 10;
            end
            
            
        end
        
        
    end %end static methods
end


