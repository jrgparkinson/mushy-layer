classdef ProblemDomain < handle
    properties
        dim
        num_ghost_x
        num_ghost_y
        num_comps
        num_levels
        max_level
        interpolationMethod %For details of different interpolation methods see documentation for resizem()
        domainExtent
        xhi
        xlo
        yhi
        ylo
        dxCoarse

        lo_i_offset % To account for situations where the bottom corner isn't at (0,0)
        lo_j_offset
        refRatio %vector of ratios. Item n contains ratio between n and n-1. First element is therefore redundant.

    end
    
    methods

        function obj = ProblemDomain(num_comps, num_levels, max_level, dim, num_ghost_x, num_ghost_y, interpolationMethod, domainExtent, dxCoarse, refRatio)

            obj.num_comps = num_comps;
            obj.num_levels = num_levels;
            obj.max_level = max_level;
            obj.dim = dim;
            obj.num_ghost_x = num_ghost_x;
            obj.num_ghost_y = num_ghost_y;
            obj.interpolationMethod = interpolationMethod;
            obj.domainExtent = domainExtent;
            obj.dxCoarse = dxCoarse;
            
            obj.xlo = (double(obj.domainExtent.lo_i) + 1/2) * obj.dxCoarse;
            obj.ylo = (double(obj.domainExtent.lo_j) + 1/2) * obj.dxCoarse;
            obj.xhi = obj.xlo + (double(obj.domainExtent.hi_i+1 - obj.domainExtent.lo_i))*obj.dxCoarse;
            obj.yhi = obj.ylo + (double(obj.domainExtent.hi_j+1 - obj.domainExtent.lo_j))*obj.dxCoarse;
            
            obj.refRatio = refRatio;
            
            obj.lo_i_offset = 0;
            obj.lo_j_offset = 0;
            
%             if obj.domainExtent.lo_i < 0
%                 obj.lo_i_offset = - obj.domainExtent.lo_i;
%             end
%             
%             if obj.domainExtent.lo_j < 0
%                 obj.lo_j_offset = - obj.domainExtent.lo_j;
%             end
            
            obj.lo_j_offset =  - obj.domainExtent.lo_j;
            obj.lo_i_offset =  - obj.domainExtent.lo_i;

        end
        
        
        function shape = domainPolyshape(obj, dx)
            dxShift = 0; %dx-obj.dxCoarse; %;
            xcoords = [obj.xlo+dxShift obj.xlo+dxShift obj.xhi-dxShift obj.xhi-dxShift];
            ycoords = [obj.ylo+dxShift obj.yhi-dxShift obj.yhi-dxShift obj.ylo+dxShift];
            
            %xcoords = [0.5 0.5 1.0 1.0];
            %ycoords = [1 2  2 1];
            
            shape = polyshape(xcoords, ycoords);
           
            
        end
    end
end