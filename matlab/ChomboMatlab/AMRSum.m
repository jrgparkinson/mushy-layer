function [L1, L2, Max, Sum] = AMRSum(ml, comp)

if length(ml.levelArray) < 1
    L1 = NaN; L2 = NaN; Max = NaN; Sum=NaN;
    return;
end

dat = ml.dataForComp(comp);

L1 = 0; L2 = 0; Max = 0; Sum = 0;

% Starting from the finest level, iterate over all boxes and
% calculate the error metric, then set the cells we've already
% counted to NaN to avoid double counting

volume_weighting =ml.finest_dx()^2;

for level_i = length(ml.levelArray):-1:1
    computedLevel = ml.levelArray(level_i);
    refinement = computedLevel.dx/ml.finest_dx();
    
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
        
        subBox = dat(lo_i:hi_i, lo_j:hi_j);
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
        dat(lo_i:hi_i, lo_j:hi_j) = NaN.*ones(num_x, num_y);
        
        temp = 0;
    end
    
    
    L1 = L1 + L1_lev;
    L2 = L2 +   L2_lev;
    Max = max(Max,   Max_lev);
    Sum = Sum +   Sum_lev;
    
    
    
end

L2 = sqrt(L2);

end