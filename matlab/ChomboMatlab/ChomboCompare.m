% For now, exact_file must be single level

classdef ChomboCompare < handle
    properties
        exact_file
        computed_file
        diff_output
        exact_func
    end
    
    properties (Constant)
        L1_ERR = 1
        L2_ERR = 2
        MAX_ERR = 3
        err_title = {'L1 error', 'L2 error', 'Max error'}
    end
    
    methods
        function obj = ChomboCompare(exact_file)
            if nargin > 0
                obj.exact_file = exact_file;
            end
        end
        
        
        %Compare the component of our computed file with the analytic
        %function provided and return an output file where this component
        %has been differenced
        function diff_output = diff_exact_func(obj, chOutput_computed, exact_func, computedComp_i,t, param)
            
            % Initialise our output
            diff_output = chOutput_computed;
            for lev = 1:length(diff_output.levelArray)
                
                diffLevel = diff_output.levelArray(lev);
                
                computed_dx_lev = diffLevel.dx;
                
                % Get all computed data on this level
                %computed_level = diffLevel.data_no_ghost();
                computed_level = chOutput_computed.levelArray(lev).data_no_ghost();
                
                %[Xlev, Ylev] = chOutput_computed.levelArray(lev).grid_no_ghost();
                %Get a grid spanning the entire domain at this level of
                %refinement
                [Xlev, Ylev] = chOutput_computed.domainGrid(lev);
                
                % Fill the entire domain at the refinement of this level
                % with the exact solution
                if nargin > 4
                    exact_data_comp = exact_func(Xlev, Ylev, t, param);
                else
                    exact_data_comp = exact_func(Xlev, Ylev);
                end
                
                % Calculate the difference between coarsened exact and
                % computed on this level
                computed_level_comp = squeeze(computed_level(computedComp_i, :, :));
                diff_level_comp = exact_data_comp - computed_level_comp;
                
                % Replace the data in each box of diff_output with the
                % actual difference between computed and exact solutions
                for box_i = 1:length(diffLevel.boxArray)
                    new_diff_box = diffLevel.boxArray(box_i);
                    ghostx = new_diff_box.problemDomain.num_ghost_x;  ghosty = new_diff_box.problemDomain.num_ghost_y;
                    
                    % Number of cells to take from this box. Default is
                    % to take all
                    num_x = new_diff_box.mx; num_y=new_diff_box.my;
                    lo_i_box = 1; lo_j_box = 1;
                    
                    % Get the extent of this box relative to the whole
                    % level
                    
                    
                    lo_i_lev = (new_diff_box.xlow + new_diff_box.dx/2)/new_diff_box.dx;
                    lo_j_lev = (new_diff_box.ylow + new_diff_box.dy/2)/new_diff_box.dy;
                    
                    % Get the offset to account for domain bottom left
                    % corner not being (0,0)
                    [i_offset, j_offset] = diffLevel.levelOffset();
                    lo_i_lev = lo_i_lev + i_offset;
                    lo_j_lev = lo_j_lev + j_offset;
                    
                    %Add some ghost cells if needed
                    if lo_i_lev < 0
                        lo_i_lev = lo_i_lev + ghostx;
                    end
                    
                    if lo_j_lev < 0
                        lo_j_lev = lo_j_lev + ghosty;
                    end
                    
                    
                    probDomain = new_diff_box.problemDomain;
                    
                    % remove some edge ghost cells if appropriate
                    
                    if new_diff_box.xlow < probDomain.xlo
                        num_x = num_x - ghostx;
                        lo_i_box = lo_i_box + ghostx;
                    end
                    
                    if new_diff_box.ylow < probDomain.ylo
                        num_y = num_y - ghosty;
                        lo_j_box = lo_j_box + ghosty;
                    end
                    
                    if new_diff_box.xhi > probDomain.xhi
                        num_x = num_x - ghostx;
                        %lo_i_box = lo_i_box - ghostx;
                    end
                    
                    if new_diff_box.yhi > probDomain.yhi
                        num_y = num_y - ghosty;
                        %lo_j_box= lo_j_box-ghosty;
                    end
                    
                    % Get interior extent of this box
                    hi_i_box = lo_i_box + num_x - 1;
                    hi_j_box = lo_j_box + num_y - 1;
                    
                    
                    %lo_i_lev = 1; lo_j_lev = 1;
                    hi_i_lev = lo_i_lev + num_x - 1;
                    hi_j_lev = lo_j_lev + num_y - 1;
                    
                    new_diff_box.data(computedComp_i, lo_i_box:hi_i_box, lo_j_box:hi_j_box) = diff_level_comp(lo_i_lev:hi_i_lev, lo_j_lev:hi_j_lev);
                    test = squeeze(new_diff_box.data(computedComp_i, :, :));
                    
                    % Save our updated box
                    diffLevel.boxArray(box_i) = new_diff_box;
                end
                
                dat = diffLevel.data_no_ghost();
                test  = squeeze(dat(computedComp_i, :, :));
                temp = 0;
                
                
                
                
                % Make sure we save our new data
                diff_output.levelArray(lev) = diffLevel;
                
                test = diffLevel.dataForComp(6, false, 1);
                temp=0;
            end
            
            
            obj.diff_output = diff_output;
        end
        
        % Return a chombo output object
        % By default, compare component 1 vs components 1, ... component n
        % vs component n
        % Optionally specify pairs of components to compare e.g. component 5 vs component 7
        % and place the difference in the place of the first component in diff_output
        function diff_output = diff(obj, chOutput_computed, componentsComputed, componentsExact)
            
            
            
            includeGhost = false;
            [X, Y] = obj.exact_file.grid();
            exact_data = obj.exact_file.data(includeGhost);
            
            if length(chOutput_computed.levelArray) < 1
                diff_output = ChomboOutputFlat(exact_data.*NaN, X, Y, obj.exact_file.finest_dx());
                return;
            end
            
            computed_data = chOutput_computed.data(includeGhost);
            
            % Default: compare like for like
            if nargin < 3
                componentsComputed = 1:length(chOutput_computed.components);
                componentsExact = componentsComputed;
            end
            
            if length(componentsComputed) ~= length(componentsExact)
                disp('ChomboCompare error: list of components to compare are not the same length');
                diff_output = ChomboOutputFlat(exact_data.*NaN, X, Y, obj.chOutput_computed.finest_dx());
                return;
            end
            
            
            % Check data sets have same number of components
            if size(exact_data, 1) ~= size(computed_data, 1)
                % This is a problem
                printf('ChomboCompare error: computed and exact files have different numbers of components');
                diff_output = ChomboOutputFlat(exact_data.*NaN, X, Y, obj.chOutput_computed.finest_dx());
                return;
            end
            
            % Our final difference object will have the same grids as
            % computed_data, but all the values will be subtracted from the
            % 'exact' solution interpolated onto these grids
            diff_output = chOutput_computed;
            
            % Remember that exact data only has one level
            %  exact_data = obj.exact_file.levelArray(1).data();
            
            
            % Iterate over levels, each of which has a different dx so
            % needs different refinement
            for lev = 1:length(diff_output.levelArray)
                
                % Get exact data to compare against
                % Two possibilities we allow:
                % 1) exact data only has one (fine level) - interpolate and
                % subtract
                % 2) exact data matches level structure of computed data -
                % do like for like comparison
                
                if length(obj.exact_file.levelArray) == 1
                    exact_data = obj.exact_file.levelArray(1).data_no_ghost();
                    exact_dx = obj.exact_file.finest_dx();
                else
                    exact_data = obj.exact_file.levelArray(lev).data_no_ghost();
                    exact_dx = obj.exact_file.levelArray(lev).dx;
                end
                
                diffLevel = diff_output.levelArray(lev);
                
                computed_dx_lev = diffLevel.dx;
                
                %refinement = computed_dx_lev/obj.exact_file.finest_dx();
                refinement = computed_dx_lev/exact_dx;
                
                % Get all computed data on this level
                %computed_level = diffLevel.data_no_ghost();
                computed_level = chOutput_computed.levelArray(lev).data_no_ghost();
                
                test = squeeze(computed_level(6, :, :));
                temp=0;
                
                % Iterate over components
                for comp_i = 1:length(componentsComputed)
                    computedComp_i = componentsComputed(comp_i);
                    exactComp_i = componentsExact(comp_i);
                    % Interpolate exact data for this component down to
                    % this level if required
                    exact_data_comp = squeeze(exact_data(exactComp_i, :, :));
                    
                    if refinement == 1
                        % Don't need to do any interpolating in this case
                        interpolated_exact_comp = exact_data_comp;
                    elseif mod(refinement, 2) ~= 0
                        disp('ChomboCompare error: expect refinement to be a factor of 2. Stuff might not work.');
                        diff_output = ChomboOutputFlat(computed_data.*NaN, X, Y, chOutput_computed.finest_dx());
                        return;
                    else
                        spacing = refinement;
                        
                        num_points = size(exact_data_comp);
                        
                        x = (1+(spacing-1)/2):spacing:(num_points(2)-(spacing-1)/2);
                        y = (1+(spacing-1)/2):spacing:(num_points(1)-(spacing-1)/2);
                        [sample_X, sample_Y] = meshgrid(x, y);
                        
                        try
                            interpolated_exact_comp = interp2(exact_data_comp, sample_X, sample_Y, 'spline');
                        catch e
                            
                            fprintf(1,'Failed to do interpolation in ChomboCompare.diff() \n%s \n',e.identifier);
                            fprintf(1,'There was an error! The message was:\n%s \n',e.message);
                            interpolated_exact_comp = NaN*sample_X;
                        end
                    end
                    
                    % Calculate the difference between coarsened exact and
                    % computed on this level
                    test = squeeze(computed_level(6, :, :));
                    temp=0;
                    
                    computed_level_comp = squeeze(computed_level(computedComp_i, :, :));
                    diff_level_comp = interpolated_exact_comp - computed_level_comp;
                    
                    % Replace the data in each box of diff_output
                    for box_i = 1:length(diffLevel.boxArray)
                        new_diff_box = diffLevel.boxArray(box_i);
                        ghostx = new_diff_box.problemDomain.num_ghost_x;  ghosty = new_diff_box.problemDomain.num_ghost_y;
                        
                        % Number of cells to take from this box. Default is
                        % to take all
                        num_x = new_diff_box.mx; num_y=new_diff_box.my;
                        lo_i_box = 1; lo_j_box = 1;
                        
                        % Get the extent of this box relative to the whole
                        % level
                        lo_i_lev = int32((new_diff_box.xlow + new_diff_box.dx/2)/new_diff_box.dx);
                        lo_j_lev = int32((new_diff_box.ylow + new_diff_box.dy/2)/new_diff_box.dy);
                        
                        %Add some ghost cells if needed
                        if lo_i_lev < 0
                            lo_i_lev = lo_i_lev + ghostx;
                        end
                        
                        if lo_j_lev < 0
                            lo_j_lev = lo_j_lev + ghosty;
                        end
                        
                        
                        probDomain = new_diff_box.problemDomain;
                        
                        % remove some edge ghost cells if appropriate
                        
                        if new_diff_box.xlow < probDomain.xlo
                            num_x = num_x - ghostx;
                            lo_i_box = lo_i_box + ghostx;
                        end
                        
                        if new_diff_box.ylow < probDomain.ylo
                            num_y = num_y - ghosty;
                            lo_j_box = lo_j_box + ghosty;
                        end
                        
                        if new_diff_box.xhi > probDomain.xhi
                            num_x = num_x - ghostx;
                            %lo_i_box = lo_i_box - ghostx;
                        end
                        
                        if new_diff_box.yhi > probDomain.yhi
                            num_y = num_y - ghosty;
                            %lo_j_box= lo_j_box-ghosty;
                        end
                        
                        % Get interior extent of this box
                        hi_i_box = lo_i_box + num_x - 1;
                        hi_j_box = lo_j_box + num_y - 1;
                        
                        
                        %lo_i_lev = 1; lo_j_lev = 1;
                        hi_i_lev = lo_i_lev + num_x - 1;
                        hi_j_lev = lo_j_lev + num_y - 1;
                        
                        test = diff_level_comp(lo_i_lev:hi_i_lev, lo_j_lev:hi_j_lev);
                        
                        new_diff_box.data(computedComp_i, lo_i_box:hi_i_box, lo_j_box:hi_j_box) = diff_level_comp(lo_i_lev:hi_i_lev, lo_j_lev:hi_j_lev);
                        test = squeeze(new_diff_box.data(computedComp_i, :, :));
                        
                        % Save our updated box
                        diffLevel.boxArray(box_i) = new_diff_box;
                    end
                    
                    dat = diffLevel.data_no_ghost();
                    test  = squeeze(dat(computedComp_i, :, :));
                    temp = 0;
                    
                    
                end
                
                % Make sure we save our new data
                diff_output.levelArray(lev) = diffLevel;
                
                test = diffLevel.dataForComp(6, false, 1);
                temp=0;
            end
            
            
            obj.diff_output = diff_output;
            
        end
        
        function comp_diff = comp_diff(obj, computed_file, comp)
            diff = obj.diff(computed_file);
            comp_diff = diff.dataForComp(comp);
        end
        
        function [L1, L2, Max, Sum] = err(obj, computed_file, comp)
            %comp_diff = obj.comp_diff(computed_file, comp);
            
            if length(computed_file.levelArray) < 1
                L1 = NaN; L2 = NaN; Max = NaN; Sum=NaN;
                return;
            end
            
            comp_diff = obj.diff_output.dataForComp(comp);
            
            L1 = 0; L2 = 0; Max = 0; Sum = 0;
            
            
            % Starting from the finest level, iterate over all boxes and
            % calculate the error metric, then set the cells we've already
            % counted to NaN to avoid double counting
            
            volume_weighting =computed_file.finest_dx()^2;
            
            for level_i = length(computed_file.levelArray):-1:1
                computedLevel = computed_file.levelArray(level_i);
                refinement = computedLevel.dx/computed_file.finest_dx();
                
                %lev_dx = computedLevel.dx;
                
                L1_lev = 0; L2_lev = 0; Max_lev = 0; Sum_lev = 0;
                
                for box_i = 1:length(computedLevel.boxArray)
                    computedBox = computedLevel.boxArray(box_i);
                    
                    [lo_i, lo_j, hi_i, hi_j] = computedBox.get_extent_no_ghost_matlab(refinement);
                    % Need to subtract ghost cells
                    %ghostx = computedBox.problemDomain.num_ghost_x;  ghosty = computedBox.problemDomain.num_ghost_y;
                    %lo_i = lo_i + ghostx - 1;
                    %lo_j = lo_j + ghosty - 1;
                    %hi_i = hi_i - ghostx - 1;
                    %hi_j = hi_j - ghosty - 1;
                    
                    num_x = hi_i+1-lo_i;      num_y = hi_j+1-lo_j;
                    
                    subBox = comp_diff(lo_i:hi_i, lo_j:hi_j);
                    % volume_weighted_phi = subBox*(lev_dx*lev_dx);
                    %volume_weighted_phi_squared = subBox.*subBox*(lev_dx*lev_dx);
                    
                    % Find all not-Nan cells
                    % num_cells = sum(sum(~isnan(subBox), 2));
                    
                    max_subBox_err = max(max(abs(subBox)));
                    
                    sum_over_box = double(nansum(nansum(abs(subBox))));
                    L2_box = nansum(nansum(subBox.*subBox));
                    % sum_over_box = double(nansum(nansum(abs(volume_weighted_phi))));
                    % L2_box = nansum(nansum(volume_weighted_phi_squared));
                    
                    
                    L1_lev = L1_lev + double(sum_over_box*volume_weighting);
                    L2_lev = L2_lev + double(L2_box*volume_weighting);
                    Max_lev = max(Max_lev, max_subBox_err);
                    Sum_lev = Sum_lev + double(nansum(nansum(subBox)) * volume_weighting);
                    
                    
                    % Remove these cells from calculations on lower levels
                    comp_diff(lo_i:hi_i, lo_j:hi_j) = NaN.*ones(num_x, num_y);
                    
                    temp = 0;
                end
                
                
                L1 = L1 + L1_lev;
                L2 = L2 +   L2_lev;
                Max = max(Max,   Max_lev);
                Sum = Sum +   Sum_lev;
                
                
                
            end
            
            
            %Assume one level
                      %  comp_diff = abs(obj.diff_output.dataForComp(comp));
                      %  sum_rows = sum(comp_diff);
                      %  sum_lev = sum(sum_rows);
                      %  err = sum_lev / (length(comp_diff)^2);
            %
            
            L2 = sqrt(L2);
            
            temp = 0;
            
        end
    end
end



