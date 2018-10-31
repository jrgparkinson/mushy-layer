classdef ChomboLevel < handle
    properties
        boxArray
        dim
        problemDomain
        X
        Y
        level
        dx
        dy
        xlow
        ylow
        xhi
        yhi
        mx
        my
    end
    methods
        function obj = ChomboLevel(level, boxArray, problemDomain)
            obj.problemDomain = problemDomain;
            obj.level = level;
            
            for i = 1:length(boxArray)
                chomboBoxes(i) = ChomboBox(boxArray(i), obj.problemDomain);
            end
            obj.boxArray = chomboBoxes;
            
            obj.dx = chomboBoxes(1).dx;
            obj.dy = chomboBoxes(1).dy;
            
            obj.xlow = chomboBoxes(1).xlow;
            obj.ylow = chomboBoxes(1).ylow;
            obj.xhi = chomboBoxes(1).xhi;
            obj.yhi = chomboBoxes(1).yhi;
            
            for i=2:length(boxArray)
                obj.xlow = min(obj.xlow, chomboBoxes(i).xlow);
                obj.ylow = min(obj.ylow, chomboBoxes(i).ylow);
                obj.xhi = max(obj.xhi, chomboBoxes(i).xhi);
                obj.yhi = max(obj.yhi, chomboBoxes(i).yhi);
            end
            
            obj.mx = 1 + (obj.xhi - obj.xlow)/(obj.dx);
            obj.my = 1 + (obj.yhi - obj.ylow)/obj.dy;
            
            [obj.X, obj.Y] = meshgrid(obj.xlow:obj.dx:obj.xhi, ...
                obj.ylow:obj.dy:obj.yhi);
            
        end
        
        function polyshapes = getMeshes(obj, domainBox)
            % lev = obj.levelArray(l);
            polyshapes = [];
            levDx = 0.5*obj.dx;
            
            for b = 1:length(obj.boxArray)
                thisBox = obj.boxArray(b);
                
                %levDx=0.75*levDx;
                
                xPoints = [thisBox.xlow-levDx thisBox.xlow-levDx thisBox.xhi+levDx thisBox.xhi+levDx];
                yPoints = [thisBox.ylow-levDx thisBox.yhi+levDx thisBox.yhi+levDx thisBox.ylow-levDx];
                
                thisPoly = polyshape(xPoints, yPoints);
                %polygons(end+1) = thisPoly;
                
                if b==1
                    mergedPolyshape = thisPoly;
                else
                    mergedPolyshape = union(mergedPolyshape, thisPoly);
                end
                
                %plot(thisPoly);
                %
            end
            
            % For now just one polyshape per level
            % Need to account for possibility that the polyshapes
            % aren't contiguous
            
            % Restrict shape to be within domain
            
            mergedPolyshape = intersect(domainBox, mergedPolyshape);
            
            
            polyshapes = mergedPolyshape;
        end
        
        function plot(obj, comp, includeGhost)
            if (nargin < 3)
                includeGhost = false;
            end
            
            if includeGhost
                plotX = obj.X;
                plotY = obj.Y;
            else
                [plotX, plotY] = obj.grid_no_ghost();
            end
            
            plot_data = obj.dataForComp(comp, includeGhost);
            
            contourf(plotX, plotY, plot_data);
            colorbar();
        end
        
        function [i_offset, j_offset] = levelOffset(obj)
            i_offset = obj.problemDomain.lo_i_offset;
            j_offset = obj.problemDomain.lo_j_offset;
            
            % If level isn't the coarsest, scale offset appropriately
            for lev = 2:obj.level
                i_offset = i_offset*obj.problemDomain.refRatio(lev);
                j_offset = j_offset*obj.problemDomain.refRatio(lev);
            end
            
        end
        
        % included for convenience as I keep forgetting to type
        % grid_no_ghost
        function  [X, Y] = grid(obj)
            [X, Y] = obj.grid_no_ghost();
        end
        
        function [X, Y] = grid_no_ghost(obj)
            
            X = obj.X(obj.problemDomain.num_ghost_x+1:(end-obj.problemDomain.num_ghost_x), ...
                obj.problemDomain.num_ghost_y+1:(end-obj.problemDomain.num_ghost_y));
            
           
            
            Y = obj.Y(obj.problemDomain.num_ghost_x+1:(end-obj.problemDomain.num_ghost_x), ...
                obj.problemDomain.num_ghost_y+1:(end-obj.problemDomain.num_ghost_y));
        end
        
        function [x_min, y_min, x_max, y_max] = grid_extent_no_ghost(obj)
            [X, Y] = obj.grid_no_ghost();
            x_min = X(1,1); y_min = Y(1,1);
            x_max = X(end, end); y_max = Y(end, end);
        end
        
        function [x_min, y_min, x_max, y_max] = grid_extent(obj)
            x_min = obj.X(1,1); y_min = obj.Y(1,1);
            x_max = obj.X(end, end); y_max = obj.Y(end, end);
        end
        
        function [lo_i, lo_j, hi_i, hi_j] = gridIndex_extent_no_ghost(obj, refinement)
            if nargin < 2
                refinement = 1;
            end
            
            %First get extent for no refinement
            numBoxes = length(obj.boxArray);
            [lo_i, lo_j, hi_i, hi_j] = obj.boxArray(1).get_extent_no_ghost();
            
            for box_i = 2:numBoxes
                box = obj.boxArray(box_i);
                
                [lo_i_new, lo_j_new, hi_i_new, hi_j_new] = box.get_extent_no_ghost();
                
                lo_i = min(lo_i, lo_i_new); lo_j = min(lo_j, lo_j_new);
                hi_i = max(hi_i, hi_i_new); hi_j = max(hi_j, hi_j_new);
            end
            
            %Now scale indexes with refinement
            lo_i = lo_i*refinement; lo_j=lo_j*refinement;
            hi_i = hi_i*refinement; hi_j=hi_j*refinement;
        end
        
        function [lo_i, lo_j, hi_i, hi_j] = gridIndex_extent(obj, refinement)
            %First get extent for no refinement
            numBoxes = length(obj.boxArray);
            [lo_i, lo_j, hi_i, hi_j] = obj.boxArray(1).get_extent();
            
            for box_i = 2:numBoxes
                box = obj.boxArray(box_i);
                
                [lo_i_new, lo_j_new, hi_i_new, hi_j_new] = box.get_extent();
                
                lo_i = min(lo_i, lo_i_new); lo_j = min(lo_j, lo_j_new);
                hi_i = max(hi_i, hi_i_new); hi_j = max(hi_j, hi_j_new);
            end
            
            %Now scale indexes with refinement
            lo_i = lo_i*refinement; lo_j=lo_j*refinement;
            hi_i = hi_i*refinement; hi_j=hi_j*refinement;
        end
        
        %matlab safe version
        function [lo_i, lo_j, hi_i, hi_j] = gridIndex_extent_matlab(obj, refinement)
            [lo_i, lo_j, hi_i, hi_j] = obj.gridIndex_extent(refinement);
            ghostx = obj.problemDomain.num_ghost_x;
            ghosty = obj.problemDomain.num_ghost_y;
            
            lo_i = lo_i + 1 + ghostx;   hi_i = hi_i + 1 + ghostx;
            lo_j = lo_j + 1 + ghosty;   hi_j = hi_j + 1 + ghosty;
        end
        
        % All data on this level, including ghost cells, at the refinement
        % specified
        function refinedData = data(obj, refinement)
            if nargin < 2
                refinement = 1;
            end
            
            
            num_x_interior_coarse = obj.problemDomain.domainExtent.hi_i+1 - obj.problemDomain.domainExtent.lo_i;
            num_y_interior_coarse = obj.problemDomain.domainExtent.hi_j+1 - obj.problemDomain.domainExtent.lo_j;
            dxCoarse = obj.problemDomain.dxCoarse;
            
            num_x_interior_lev = num_x_interior_coarse * (dxCoarse/obj.dx);
            num_y_interior_lev = num_y_interior_coarse * (dxCoarse/obj.dx);
            
            ghostx = obj.problemDomain.num_ghost_x;
            ghosty = obj.problemDomain.num_ghost_y;
            
            num_x =  num_x_interior_lev + 2*ghostx;
            num_y =  num_y_interior_lev + 2*ghosty;
            
            % this 'data' object covers the entire domain on this level
            data =  nan * ones(obj.problemDomain.num_comps, ...
                int16(num_x), ...
                int16(num_y));
            
            
            %Fill with data from each box
            numBoxes = length(obj.boxArray);
            for box_i = 1:numBoxes
                box = obj.boxArray(box_i);
                box_data = box.data;
                [lo_i, lo_j, hi_i, hi_j] = box.get_extent_matlab();
                %                 temp = data(:, lo_i:hi_i, lo_j:hi_j);
                data(:, lo_i:hi_i, lo_j:hi_j) = box_data();
            end
            
            
            
           % theta = squeeze(data(6, :,:));
            
            %Now refine it (but only refine interior points)
            %Do nearest neighbour interpolation onto the finest grid
            
            %             refinedData = nan* ones(obj.problemDomain.num_comps, ...
            %                  (obj.mx-2*ghostx)*refinement + 2*ghostx, ...
            %                  (obj.my-2*ghosty)*refinement + 2*ghosty);
            refinedData = nan* ones(obj.problemDomain.num_comps, ...
                num_x_interior_lev*refinement + 2*ghostx, ...
                num_y_interior_lev*refinement + 2*ghosty);
            
           
            
            for var_i = 1:obj.problemDomain.num_comps
              [loi, loj, hii, hij] =  obj.gridIndex_extent_no_ghost(1);
              % Need domain extent on this level
              
               %loi = obj.problemDomain.domainExtent.lo_i + 1;
               %loj = obj.problemDomain.domainExtent.lo_j + 1;
               %hii = obj.problemDomain.domainExtent.hi_i + 1;
               %hij = obj.problemDomain.domainExtent.hi_j + 1;
               
               [i_offset, j_offset] = obj.levelOffset();
             
               if hii > size(data, 2) || hij > size(data, 3) || loi < 1 || loj < 1
                  blah = 0; % catch this
               end
               
                %componentData = squeeze(data(var_i, ...
                 %   loi:hii,...
                %    loj:hij));
                
                % Just get all the data on this level for this variable
                % componentData spans the entire domain, including
                % uncovered regions
                componentData = squeeze(data(var_i, ...
                    :,...
                    :));
                
                %MATLAB likes to be annoying and replace NaN values by
                %1 when doing resizem if the other values are between 0
                %and 1 themselves, this is a way around that. This is a
                %bit ridiculous.
                largeNumber = 1e5;
                componentData = componentData * largeNumber;
                
                %First refine everything. newData will cover the entire
                %domain, including uncovered regions
                newData = resizem(componentData, refinement, obj.problemDomain.interpolationMethod);
                
                % Get back our actual values
                newData = newData / largeNumber;
                
                trim_ghostx = (refinement-1)*ghostx;
                trim_ghosty = (refinement-1)*ghosty;
                
                
               % [loi_ref, loj_ref, hii_ref, hij_ref] = obj.gridIndex_extent_matlab(refinement);
               % hii_ref = hii_ref; hij_ref = hij_ref;
%                testing =  newData(1+trim_ghostx:end-trim_ghostx,  1+trim_ghosty:end- trim_ghosty);
%                test2 = squeeze(refinedData(var_i,loi_ref:hii_ref,loj_ref:hij_ref));
               
                
%                if size(testing, 2) ~= size(test2 , 2) || ...
%                    size(testing, 1) ~= size(test2, 1)
%                    balh = 0; % catch this
%                end
                
                %Now remove the extra ghost cells we've created
                test = squeeze(refinedData(var_i,:,:));
                refinedData(var_i,:,:) = newData(1+trim_ghostx:end-trim_ghostx,  1+trim_ghosty:end- trim_ghosty);
            end
            
            %theta = squeeze(refinedData(6, :,:));
            %tempVar = 0;
        end
        
        
        function data_no_ghost = data_no_ghost(obj, refinement)
            if (nargin < 2)
                refinement = 1;
            end
            
            data = obj.data(refinement);
            
            ghostx = obj.problemDomain.num_ghost_x;
            ghosty = obj.problemDomain.num_ghost_y;
            
            num_x_interior_coarse = obj.problemDomain.domainExtent.hi_i+1 - obj.problemDomain.domainExtent.lo_i;
            num_y_interior_coarse = obj.problemDomain.domainExtent.hi_j+1 - obj.problemDomain.domainExtent.lo_j;
            dxCoarse = obj.problemDomain.dxCoarse;
            
            num_x_interior_lev = num_x_interior_coarse * (dxCoarse/obj.dx);
            num_y_interior_lev = num_y_interior_coarse * (dxCoarse/obj.dx);
            
            data_no_ghost = nan*ones(obj.problemDomain.num_comps, ...
                refinement*int16(num_x_interior_lev), ...
                refinement*int16(num_y_interior_lev));
            
            %Now remove ghost cells
            
            for var_i = 1:obj.problemDomain.num_comps
                data_no_ghost(var_i, :,:) = data(var_i, ghostx+1:end-ghostx, ghosty+1:end-ghosty);
            end
        end
        
        function dat = getData(obj, includeGhost, refinement)
            if includeGhost
                dat = obj.data(refinement);
            else
                dat = obj.data_no_ghost(refinement);
            end
        end
        
        function dataForComp = dataForComp(obj, comp, includeGhost, refinement)
            if nargin < 4
                refinement = 1;
            end
            if (nargin < 3)
                includeGhost = false;
            end
            
            data = obj.getData(includeGhost, refinement);
            dataForComp = squeeze(data(comp,:,:))';
        end
        
        function patches = get_patches(obj)
            half_dx = obj.dx/2;
            
            
            
            for j = 1:length(obj.boxArray)
                [Xbox,Ybox] = obj.boxArray(j).grid_no_ghost();
                x_min = Xbox(1,1)-half_dx; y_min = Ybox(1,1)-half_dx;
                x_max = Xbox(1,end)+half_dx; y_max = Ybox(end, 1)+half_dx;
                width = x_max-x_min; height = y_max-y_min;
                patches(j, :) = [x_min, y_min, x_max, y_max];
                
            end
            
            % Loop over and merge patches
            
            prev_patches = [];
            
            while isequal(prev_patches, patches) == 0
                prev_patches = patches;
                
                for j=1:size(patches, 1)
                    % Loop through all other patches to look for neighbours
                    topLeftj = [patches(j, 1) patches(j, 4)];
                    topRightj = [patches(j, 3) patches(j, 4)];
                    bottomLeftj = [patches(j, 1) patches(j, 2)];
                    bottomRightj = [patches(j, 3) patches(j, 2)];
                    for i=1:size(patches, 1)
                        topLefti = [patches(i, 1) patches(i, 4)];
                        topRighti = [patches(i, 3) patches(i, 4)];
                        bottomLefti = [patches(i, 1) patches(i, 2)];
                        bottomRighti = [patches(i, 3) patches(i, 2)];
                        
                        
                        if ~isequal(patches(i, :), patches(j, :))
                            % Check if patches share a common side
                            if (isequal(topLeftj, topRighti) && isequal(bottomLeftj, bottomRighti)) ||  ...
                                    (isequal(topLeftj , bottomLefti) && isequal(topRightj,bottomRighti)) || ...
                                    (isequal(topRightj , topLefti) && isequal(bottomRightj, bottomLefti)) || ...
                                    (isequal(bottomLeftj , topLefti) && isequal(bottomRightj , topRighti) )
                                
                                patches(j, 1) = min(patches(j, 1), patches(i,1));
                                patches(j, 2) = min(patches(j, 2), patches(i,2));
                                
                                patches(j, 3) = max(patches(j, 3), patches(i,3));
                                patches(j, 4) = max(patches(j, 4), patches(i,4));
                                
                                patches(i, :) = -1;
                                
                                % Need to re define this
                                topLeftj = [patches(j, 1) patches(j, 4)];
                                topRightj = [patches(j, 3) patches(j, 4)];
                                bottomLeftj = [patches(j, 1) patches(j, 2)];
                                bottomRightj = [patches(j, 3) patches(j, 2)];
                            end
                        end
                    end
                    
                    
                    
                end
                
                
                
            end
        end
        
    end
end
