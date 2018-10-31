classdef Box < handle
    properties
        lo_i
        hi_i
        lo_j
        hi_j
        dx
        values
        values_comp
        resized_data
        resized_data_comp
        num_components
    end
    
    methods
        function box = Box(lo_i, hi_i, lo_j, hi_j, dx, data_comp, num_components)
           box.lo_i = lo_i; box.lo_j = lo_j;  box.hi_i = hi_i;  box.hi_j = hi_j;
           box.dx = dx;
           box.num_components = num_components;
           
           s = size(data_comp); num_comp = s(1);
           for comp = 1:num_comp
               new = reshape(data_comp(comp, :, :), [hi_i+1-lo_i, hi_j+1-lo_j]).';
               vals_comp(comp, :, :) = new;
           end
           box.values_comp = vals_comp;

           box.values = reshape(data_comp(1, :, :), [hi_i+1-lo_i, hi_j+1-lo_j]);
               
        end
        
        % Helper function to get around matlab indexing
        function [lo_i, hi_i, lo_j, hi_j] = get_grid_extent_plus_one(obj, dx)
            if nargin < 2
                dx = obj.dx;
            end
            [lo_i, hi_i, lo_j, hi_j] = get_grid_extent(obj, dx);
            
             lo_i = lo_i + 1;
             lo_j = lo_j + 1;
             hi_i = hi_i + 1;
             hi_j = hi_j + 1;
        end
        
        function [lo_i, hi_i, lo_j, hi_j] = get_grid_extent(obj, dx)
            if nargin < 2
                dx = obj.dx;
            end
            
            ref = obj.dx/dx;
            
            lo_i = ref*obj.lo_i;
            lo_j = ref*obj.lo_j;
                        
            hi_i = obj.hi_i;
            hi_j = obj.hi_j;
            
            for iter = 1:log2(ref)
                hi_i = 2*hi_i + 1;
                hi_j = 2*hi_j + 1;
            end
            
        end
        
        function [X, Y] = get_grid(obj, dx)
            if nargin < 2
                dx = obj.dx;
            end
            [new_lo_i, new_hi_i, new_lo_j, new_hi_j] = obj.get_grid_extent(dx);
            
            x_hi = (new_hi_i + 0.5)*dx; 
            y_hi = (new_hi_j + 0.5)*dx;
            
            x_lo = (new_lo_i + 0.5)*dx; 
            y_lo = (new_lo_j + 0.5)*dx;
            

            [X, Y] = ndgrid(x_lo:dx:x_hi, y_lo:dx:y_hi);
        end
        
        function values = get_values(obj, required_dx, component, interpolation)
           if (nargin < 4)
               interpolation = 'none';
           end
           
           obj.resize_data(required_dx, interpolation);
           values = squeeze(obj.resized_data_comp(component, :, :)).';
        end
        
        %Interpolation could either be a method e.g. cubic, or just leave
        %blank/'none' 
        function resize_data(obj, required_dx, interpolation)

            [Xold, Yold] = obj.get_grid(obj.dx);
            [Xnew, Ynew] = obj.get_grid(required_dx);
            
            % Do this for every component
            %resized_data_comp = nan*ones(obj.num_components, 1);
            for comp = 1:obj.num_components
                if (nargin < 3 || strcmp(interpolation,'none') == true)
                    %Just fill new grid points with same value as surrounding
                    %points
                    old_data = squeeze(obj.values_comp(comp, :, :));
                    refined_data = squeeze(obj.values_comp(comp, :, :));

                    %Keep doubling the grid size
                    for refinements = 1:log2(obj.dx/required_dx)

                        for i = 1:size(old_data,1)
                            for j = 1:size(old_data, 2)                       
                                refined_data(2*i-1:2*i , 2*j-1:2*j) = old_data(i,j).*ones(2,2);
                            end
                        end
                        old_data = refined_data;
                    end

                else
                    F = griddedInterpolant(Xold, Yold, obj.values_comp(comp, :, :), interpolation);
                    refined_data = F(Xnew, Ynew);                 
                end
                
                resized_data_comp(comp, :, :) = refined_data;
            end
            
            obj.resized_data = refined_data;
            
            obj.resized_data_comp = resized_data_comp;
            
        end

    end
    
end