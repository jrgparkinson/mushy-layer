classdef Plotfile < handle
    properties
        filename
        num_levels
        num_components
        level_data
        merged_grid
    end
    
    methods
        function pf = Plotfile(filename)
           pf.filename = filename;
           
           %Read in data from the HDF5 file
           pf.num_levels = h5readatt(filename, '/', 'num_levels');
           pf.num_components = h5readatt(filename, '/', 'num_components');
           
           %Read in the levels. Note that Chombo indexing starts at 0,
           %but matlab starts at 1
           
           %level 0 is special as it's box covers the entire domain - need
           %to know this
           ld_arr(1) = LevelData(pf.filename, 0);
           base_grid = ld_arr(1).get_whole_level_grid();
           
           for level = 2:pf.num_levels
               ld_arr(level) = LevelData(pf.filename, level-1, ld_arr(1));
           end
           pf.level_data = ld_arr;
           
        end
        
        function dx = get_finest_dx(obj)
           finest_level = obj.level_data(end);
           dx = finest_level.dx;
        end
        
        function resize_data(obj, dx)
           for level=obj.level_data
              level.resize_data(dx) 
           end
        end
        
        function plot(obj, component, color_bar_range)
            if (nargin < 3)
                color_bar_range = obj.get_value_range(component);
            end
            
            %color_bar_range = [0,2];
            
            finest_dx = obj.get_finest_dx();

            %Convert every level to exist at this dx
            obj.resize_data(finest_dx);
            
            %Create merged grid
            obj.mergeAll()
            
            %Plot grid
            
            [X, Y] = obj.getWholeDomainGrid(finest_dx);      
            figure;
            
            data = obj.get_merged_grid(component);
            
            pcolor(X,Y,data); shading flat; colorbar; 
            
            % Matlab throws an error if this isn't true. Let it determine
            % it's own range in this case.
            if (color_bar_range(2) > color_bar_range(1))
                caxis(color_bar_range);
            end
            
            %pcolor(X,Y,obj.merged_grid); shading interp;
            
        end
        
        
        
        function mergeAll(obj, dx)
            if nargin < 2
                dx = obj.get_finest_dx();
            end
            
            [X,~] = obj.getWholeDomainGrid(dx);
            
            for comp = 1:obj.num_components
                
                values = nan .* X;

                for ld = obj.level_data            
                    for box = ld.boxes
                       values = LevelData.add_box_to_mesh(X,box,comp, values);
                    end
                end
                
                merged(comp, :, :) = values;
            
            end

            obj.merged_grid = merged;
            
        end
        
        function gr = get_merged_grid(obj, component, dx)
            if nargin < 3
                dx = obj.get_finest_dx();
            end
            
            if size(obj.merged_grid, 1) == 0
                obj.mergeAll(dx);
            end
            gr = squeeze(obj.merged_grid(component, :, :));
        end
        
         function [X,Y] = getWholeDomainGrid(obj, dx)
           %box = obj.level_data(1).boxes(1);
           %[X, Y] = box.get_grid(dx);
           [X, Y] = obj.level_data(1).get_whole_level_grid();
         end
        
         function range = get_value_range(obj, component)
             pf_min = nan; pf_max = nan;
             for ld = obj.level_data
                ld_min = min(ld.raw_data_comp(component,:));      ld_max = max(ld.raw_data_comp(component,:));
                pf_min = min(pf_min, ld_min);   pf_max = max(pf_max, ld_max);
             end
             
             range = [pf_min, pf_max];
         end
       
         %Calculate the element by element absolute difference between two
         %files
         function diff = abs_diff(pf1, pf2)
             grid1 = pf1.get_merged_grid();
             grid2 = pf2.get_merged_grid();
             diff = abs(grid2-grid1);
         end
        
    end
end