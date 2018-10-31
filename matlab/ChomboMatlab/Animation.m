classdef Animation < handle
    properties
        plotfiles
        num_files
        frames
        slice_frames
    end
    
    methods

        function anim = Animation(base_name, num_files)
            anim.num_files = num_files;
            
            for j = 0:num_files
                filename = strcat(base_name,sprintf('%02d',j), '.2d.hdf5');
                plotfile_arr(j+1) = Plotfile(filename);
            end
            
            anim.plotfiles = plotfile_arr;
           
        end
        
        function range = get_value_range(obj)
           an_min = nan; an_max = nan;
             for pf = obj.plotfiles
                pf_range = pf.get_value_range();
                an_min = min(an_min, pf_range(1));
                an_max = max(an_max, pf_range(2));
             end 
             
             range = [an_min, an_max];
        end
        
        %Type: 0 - whole grid
        % 1 - a slice at the x index specified
        function make_movie(obj, type, index)
            if (nargin < 2)
                type = 0;
            end
           F(obj.num_files) = struct('cdata',[],'colormap',[]);
           for j = 1:obj.num_files
               plotfile = obj.plotfiles(j);
               
               if (type == 0) 
                color_bar_range = obj.get_value_range();
                plot(plotfile, color_bar_range);
               elseif (type == 1)                    
                    grid = plotfile.get_merged_grid;
                    if (index < size(grid, 1))
                        slice = grid(index,:);
                        plot(slice)
                        axis([0,256,0,2])
                    end
               else
                   error('Animation::make_movie Invalid plottig type specified');
               end
               F(j) = getframe(gcf);
           end
           
           obj.frames = F;
        end
        
        function play_movie(obj,n_times)
            fig = figure;
            movie(fig,obj.frames,n_times)
        end
        
       
        
    end
    
end