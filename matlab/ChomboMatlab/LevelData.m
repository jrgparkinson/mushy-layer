classdef LevelData < handle
    properties
        offset
        boxes
        data
        raw_data  
        raw_data_comp
        dx
        dt
        ref_ratio
        X
        Y
        comp_names
    end
    
    properties (Constant)
        LO_I = 'lo_i'; 
        HI_I = 'hi_i'; 
        LO_J = 'lo_j'; 
        HI_J = 'hi_j';
    end
    
    methods     
        
        function ld = LevelData(filename, level, level0)
            
            %Read all the data in
            ld.offset = h5read(filename,strcat('/level_' , num2str(level) , '/data:offsets=0'));
            box_limits = h5read(filename,strcat('/level_' , num2str(level) , '/boxes'));
            ld.raw_data = h5read(filename,strcat('/level_' , num2str(level) , '/data:datatype=0'));
            
            %Attributes
            ld.dx = h5readatt(filename, strcat('/level_' , num2str(level)), 'dx');
            ld.dt = h5readatt(filename, strcat('/level_' , num2str(level)), 'dt');
            ld.ref_ratio = h5readatt(filename, strcat('/level_' , num2str(level)), 'ref_ratio');
            
            %Get the different components
            num_components = h5readatt(filename, '/', 'num_components');
            
            component_names = cell(num_components, 1);    
            for comp=1:num_components
                name = h5readatt(filename, '/', strcat('component_', int2str(comp-1)));
                component_names{comp} = [name];
            end
            ld.comp_names = component_names;
            
            num_boxes = length(box_limits.(LevelData.LO_I));
            
            %Parse the data
            % The data is simply a 1d list which covers every component and
            % every box in turn
            
            %We can have multiple boxes at one level
            box_data = ld.raw_data;
            
            %We also have multiple 
            num_vals = length(box_data);
            vals_per_box = num_vals/num_boxes;
            vals_per_comp = num_vals/num_components;
            vals_per_comp_per_box = vals_per_box/num_components;
            
            
            box_data_comp = nan*ones(num_components, vals_per_comp);
            
            for box_iter=1:num_boxes
                box_offset = (box_iter-1) * vals_per_box;
                for comp=1:num_components
                    source_from = box_offset + (comp-1)*vals_per_comp_per_box + 1;
                    source_to = box_offset + comp*vals_per_comp_per_box;
                    
                    dest_from = (box_iter-1)*vals_per_comp_per_box + 1;
                    dest_to = dest_from + vals_per_comp_per_box - 1;

                    box_data_comp(comp, dest_from:dest_to) =  box_data(source_from:source_to).';
                end
            end
            
            ld.raw_data_comp = box_data_comp;
            
            
            
            start_index = 1; %The index to start extracting values from
            for iter = 1:num_boxes
                
                
                lo_i = double(box_limits.(LevelData.LO_I)(iter));
                hi_i = double(box_limits.(LevelData.HI_I)(iter));
                
                lo_j = double(box_limits.(LevelData.LO_J)(iter));
                hi_j = double(box_limits.(LevelData.HI_J)(iter));
                
                %work out how many values to take for this box
                num_values = (hi_i + 1 - lo_i) * (hi_j + 1 - lo_j);
                
                %Create empty vector to store the data for each component
                this_box_data_comp = nan*ones(num_components, num_values);
                
                %Take the values for this box from the large array
                for comp=1:num_components
                    temp = box_data_comp(comp, start_index:(start_index+num_values-1)).';
                    this_box_data_comp(comp, :) = temp;
                end
                
                %this_box_data = box_data(1:num_values);
                %box_data = box_data(num_values+1:end);
                
                box_arr(iter) = Box(lo_i, hi_i, lo_j, hi_j, ld.dx, this_box_data_comp, num_components);
                
                start_index = start_index + num_values;
            end
            
            ld.boxes = box_arr;
            
            if (nargin < 3)
                %This is level = 0
                [X, Y] = ld.get_whole_level_grid();
            else
                %This is level > 0
                [X, Y] = level0.get_whole_level_grid(ld.dx);
            end
            
            %The base grid covers the whole domain at the dx of this level
            ld.X = X; ld.Y = Y;
            
        end
 
        
        %plot the data on this level. need to be careful because the box
        %layout may be disjoint. We have the base grid, so just need to plot
        %the value we have in the relevant places on it.
        function plot(obj, component)
            values = nan .* obj.X;
            
            for box = obj.boxes
               %[lo_i, hi_i, lo_j, hi_j] = box.get_grid_extent_plus_one();
               %values(lo_i:hi_i, lo_j:hi_j) = box.values;
               values = LevelData.add_box_to_mesh(obj.X, box, component, values);
            end
            
            contourf(obj.X, obj.Y, values);
        end
        
        %Plot a single box from this level on the same mesh as this level uses
        function plot_box(obj, box, component)         
            values = LevelData.add_box_to_mesh(obj.X, box, component);
            contourf(obj.X, obj.Y, values);
        end
        
        function [X, Y] = get_whole_level_grid(obj, dx)
            if nargin < 2
                dx = obj.dx;
            end
            
            %Iterate over all boxes and merge them to create one
            % Grid which covers them all
            
            % Start with first box
            [lo_i, hi_i, lo_j, hi_j] = obj.boxes(1).get_grid_extent(dx);
           
            for box = obj.boxes
                 [new_lo_i, new_hi_i, new_lo_j, new_hi_j] = box.get_grid_extent(dx);
                 lo_j = min(lo_j, new_lo_j);
                 lo_i = min(lo_i, new_lo_i);
                 hi_j = max(hi_j, new_hi_j);
                 hi_i = max(hi_i, new_hi_i);
            end
            
                      
            x_hi = (hi_i + 0.5)*dx; 
            y_hi = (hi_j + 0.5)*dx;
            
            x_lo = (lo_i + 0.5)*dx; 
            y_lo = (lo_j + 0.5)*dx;
            

            [X, Y] = ndgrid(x_lo:dx:x_hi, y_lo:dx:y_hi);
        end
        
   
        function resize_data(obj, dx)
           for box = obj.boxes
              box.resize_data(dx) 
           end
        end
                  
    end
    
    methods(Static)
        %Take any box and add it to a particular mesh
        function values = add_box_to_mesh(X, box, component, values)
            if (nargin < 4)
                values = nan .* X;
            end
            
            required_dx = X(2,1) - X(1,1);
            [lo_i, hi_i, lo_j, hi_j] = box.get_grid_extent_plus_one(required_dx);
            
            %Hard code interpolation
            interpolation = 'none'; %'cubic';
            %component = 1;
            new_values = box.get_values(required_dx, component, interpolation);
            
            values(lo_i:hi_i, lo_j:hi_j) = new_values;
            
        end
    end
    
end