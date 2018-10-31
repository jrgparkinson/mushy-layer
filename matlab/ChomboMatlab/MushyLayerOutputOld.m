% This is a class to wrap up a generic ChomboOutput class with specific MushyLayer 
% functionality
classdef MushyLayerOutputOld < handle
    properties
        chomboOutput
        mlComps
    end
    
    methods
        function obj = MushyLayerOutputOld(dim, frame, output_dir, plot_prefix)
            obj.mlComps = MushyLayerComponents();
            obj.chomboOutput = ChomboOutput(dim, frame, output_dir, plot_prefix);
        end
        
        %This currently only calculates based on the base level
        function [Nu_av, Nu_base] = nusseltNumber(obj)
            levels = obj.chomboOutput.levelArray;
            
            %check we have some data
            if length(levels) < 1
                Nu_av = nan; Nu_base = nan;
                return
            end
            
            coarseLevel = levels(1);
            
           dx = coarseLevel.dx;

            theta =coarseLevel.dataForComp(obj.mlComps.theta);
            Vz = coarseLevel.dataForComp(obj.mlComps.zFluidVel);
            Vx = coarseLevel.dataForComp(obj.mlComps.xFluidVel);
            [X, Y] = coarseLevel.grid_no_ghost();
            X_arr = X(1,:);
            Y_arr = Y(:, 1);
            
            [~, dtheta_dz] = gradient2order(theta, dx,dx);
            
            L = X_arr(end) - X_arr(1);

            integrand = (-dtheta_dz + Vz.*theta)/L;
            
%             figure(1);
%             %contourf(X, Z, theta);
%             %colorbar();
%             plot(Y_arr, theta(:, 29));
%             title('theta(y) cross section');
%             
%             figure(10);
%             quiver(X, Y, Vx, Vz);
%             title('Fluid velocity');
% 
%             figure();
%             contourf(X, Y, integrand);
%             colorbar();
%             title('Local Nu');

            %Calculate Nu(y) by integrating in the x direction for each y point
            
            Nu = trapz(integrand.')*dx;

%             figure(3);
%             plot(Y(:, 1), Nu);
%             title('Nu(y)');

            %Now average over each y level
            Nu_av = mean(Nu);
%              figure();
%             plot(Y_arr, Nu);
%             title('Nu(y)');
%             xlabel('y'); ylabel('Nu');
            
        
            
           % integrand = - dtheta_dz.'/L;
            integrand = (- dtheta_dz.' + (Vz.*theta).' )/(L);
            Nu = trapz(integrand)*dx;
            
%             figure();
%             plot(X_arr, integrand(1,:));
%             title('Nu(y=0)');
%             xlabel('x'); ylabel('Nu');
            
           Nu_base = Nu(1);
           

        end
        
        function plot(obj, var)
            obj.chomboOutput.plot(var)
        end
        
    end
    
    
    
end