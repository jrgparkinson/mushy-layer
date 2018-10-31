classdef ChomboBox < handle
    properties
        gridno
        level
        dx
        dy
        dz
        mx
        my
        mz
        mx_interior
        my_interior
        xlow
        ylow
        zlow
        xhi
        yhi
        zhi
        data
        problemDomain
        X
        Y
    end
    methods
        function obj = ChomboBox(box, problemDomain)
            obj.problemDomain = problemDomain;
            
            obj.gridno = box.gridno;
            obj.level = box.level;
            obj.dx = box.dx;
            obj.dy = box.dy;
            
            obj.mx = double(box.mx);
            obj.my = double(box.my);
            
            obj.mx_interior = double(obj.mx-2*obj.problemDomain.num_ghost_x);
            obj.my_interior = double(obj.my-2*obj.problemDomain.num_ghost_y);
            
            obj.xlow = double(box.xlow) + obj.dx/2;
            obj.ylow = double(box.ylow) + obj.dy/2;
            
            obj.data = box.data;
            
            if obj.problemDomain.dim == 3
                obj.dz = box.dz;
                obj.mz = box.mz;
                obj.zlow = box.zlow + obj.dz/2;
            end
            
            
            %Reshape data into (ncomps) x (nx) x (ny)
            new_data = nan * ones(obj.problemDomain.num_comps, obj.mx, obj.my);
            for ncomp = 1:obj.problemDomain.num_comps
                old_data = obj.data(ncomp, :);
                new_data(ncomp, :, :) = reshape(obj.data(ncomp, :), obj.mx, obj.my);
            end
            
            obj.data = new_data;
            
            
            obj.xhi = obj.xlow + obj.dx*(obj.mx-1);
            obj.yhi = obj.ylow + obj.dy*(obj.my-1);
            [obj.X, obj.Y] = meshgrid(obj.xlow : obj.dx : obj.xhi, ...
                obj.ylow : obj.dy : obj.yhi);


        end
        
        function dataForComp = dataForComp(obj, comp)
            data_no_ghost = obj.data_no_ghost();
            dataForComp = squeeze(data_no_ghost(comp,:,:))';
        end
        
        function plot(obj, comp)
            [X_no_ghost, Y_no_ghost] = obj.grid_no_ghost();
            this_comp_data = obj.dataForComp(comp);
            
            figure(1);
            contourf(X_no_ghost, Y_no_ghost, this_comp_data);
            colorbar();
        end
        
        function [X,Y] = grid_no_ghost(obj)
           X = obj.X(obj.problemDomain.num_ghost_x+1:(end-obj.problemDomain.num_ghost_x), ...
           obj.problemDomain.num_ghost_y+1:(end-obj.problemDomain.num_ghost_y));
           Y = obj.Y(obj.problemDomain.num_ghost_x+1:(end-obj.problemDomain.num_ghost_x), ...
           obj.problemDomain.num_ghost_y+1:(end-obj.problemDomain.num_ghost_y));
        end
        
        function data = data_no_ghost(obj, refinement)
            if nargin < 2
                refinement = 1;
            end
            
            ghostx = obj.problemDomain.num_ghost_x;
            ghosty = obj.problemDomain.num_ghost_y;
            
           data = obj.data(:, ... 
               ghostx+1 : (end-ghostx), ... 
               ghosty+1 : (end-ghosty));
           
           
           %Now refine if necessary
           if refinement > 1
               
               % Get all data including ghost cells for accurate
               % interopolation at boundaries
               data = obj.data;
           
               refinedData = nan*ones(obj.problemDomain.num_comps, int16(obj.mx-2*ghostx)*refinement, int16(obj.my-2*ghosty)*refinement);
           
                for var_i = 1:obj.problemDomain.num_comps
                    componentData = squeeze(data(var_i,:,:));
                    
                    %MATLAB likes to be annoying and replace NaN values by
                    %1 when doing resizem if the other values are between 0
                    %and 1 themselves, this is a way around that. This is a
                    %bit ridiculous.
                    largeNumber = 1e5;
                    componentData = componentData * largeNumber;
                    
                    %First refine everything
                    newData = resizem(componentData, refinement, obj.problemDomain.interpolationMethod);
                    
                    % Get back our actual values
                    newData = newData / largeNumber;
                     
                    temp2 = newData(ghostx*2+1:end-2*ghostx,...
                        ghosty*2+1:end-2*ghosty);
                    temp = refinedData(var_i,:,:);
                    refinedData(var_i,:,:) = newData(ghostx*2+1:end-2*ghostx,...
                        ghosty*2+1:end-2*ghosty);
                end
                
                data = refinedData;
           end
           
           
        end
        
    
       
        
        %convert xlow, ylow etc. into lo_i, lo_j
        function [lo_i, lo_j, hi_i, hi_j] = get_extent(obj, refinement)
            if nargin < 2
                refinement = 1;
            end
            
            % Get the domain offset on this level
            i_offset_lev = obj.problemDomain.lo_i_offset;
            j_offset_lev = obj.problemDomain.lo_j_offset;
            
            for lev_i = 2:obj.level
                i_offset_lev = i_offset_lev*obj.problemDomain.refRatio(lev_i);
                j_offset_lev = j_offset_lev*obj.problemDomain.refRatio(lev_i);
            end
            
            lo_i = (obj.xlow-obj.dx/2)/obj.dx + i_offset_lev;
            lo_j = (obj.ylow-obj.dy/2)/obj.dy + j_offset_lev;
            hi_i = lo_i + obj.mx-1;
            hi_j = lo_j + obj.my-1;
            
            % Remove ghost cells
            ghostx = obj.problemDomain.num_ghost_x;
            ghosty = obj.problemDomain.num_ghost_y;
            
            lo_i = lo_i + ghostx;   hi_i = hi_i - ghostx;
            lo_j = lo_j + ghosty;   hi_j = hi_j - ghosty;
            
            new_mx_interior = obj.mx_interior*refinement;
            new_my_interior = obj.my_interior*refinement;
            %Scale with refinement
            if refinement ~= 1
                lo_i = lo_i * refinement;
                lo_j = lo_j * refinement;

%                 hi_i = hi_i * refinement - 1;
%                 hi_j = hi_j * refinement - 1;
                hi_i = lo_i + new_mx_interior-1;
                hi_j = lo_j + new_my_interior-1;
            end
            
            % Add ghost cells back on
            lo_i = lo_i - ghostx;   hi_i = hi_i + ghostx;
            lo_j = lo_j - ghosty;   hi_j = hi_j + ghosty;
        end
        
        %Like get_extent but MATLAB safe (so indexes are all > 0)
        function [lo_i, lo_j, hi_i, hi_j] = get_extent_matlab(obj, refinement)
            if nargin < 2
                refinement = 1;
            end
            
            [lo_i, lo_j, hi_i, hi_j] = obj.get_extent(refinement);
            
            % Add 1 and add ghost cells
            ghostx = obj.problemDomain.num_ghost_x;
            ghosty = obj.problemDomain.num_ghost_y;
            
            lo_i = lo_i + 1 + ghostx;   hi_i = hi_i + 1 + ghostx;
            lo_j = lo_j + 1 + ghosty;   hi_j = hi_j + 1 + ghosty;
        end
        
        function [lo_i, lo_j, hi_i, hi_j] = get_extent_no_ghost(obj, refinement)
            if nargin < 2
                refinement = 1;
            end
            
            [lo_i, lo_j, hi_i, hi_j] = obj.get_extent(refinement);
            matlab_offset = 1; 
            lo_i = lo_i + obj.problemDomain.num_ghost_x + matlab_offset;
            lo_j = lo_j + obj.problemDomain.num_ghost_y + matlab_offset;
            hi_i = hi_i - obj.problemDomain.num_ghost_x + matlab_offset;
            hi_j = hi_j - obj.problemDomain.num_ghost_y + matlab_offset;
        end
        
        function [lo_i, lo_j, hi_i, hi_j] = get_extent_no_ghost_matlab(obj, refinement)
            if nargin < 2
                refinement = 1;
            end
            
            [lo_i, lo_j, hi_i, hi_j] = obj.get_extent_no_ghost(refinement);
           
            % Convert to Matlab co-ords:
            % add ghost+1 to both coords
            
            %Pretty sure get_extent_no_ghost is already matlab safe
            
%             lo_i = lo_i + obj.problemDomain.num_ghost_x + 1;
%             lo_j = lo_j + obj.problemDomain.num_ghost_y + 1;
%             hi_i = hi_i + obj.problemDomain.num_ghost_x + 1;
%             hi_j = hi_j + obj.problemDomain.num_ghost_y + 1;
            
        end
        
    end
    methods (Static)
        function [lo_i, lo_j, hi_i, hi_j] = make_extent_matlab_safe(lo_i, lo_j, hi_i, hi_j)
           %Need to make sure we don't have indices < 1
           
            i_offset = max(1 - lo_i, 0);
            j_offset = max(1 - lo_j, 0);
            
            lo_i = lo_i + i_offset;     hi_i = hi_i + i_offset;
            lo_j = lo_j + j_offset;     hi_j = hi_j + j_offset;
        end
    end
end