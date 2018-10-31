classdef ChomboOutputFlat
    properties
        data
        X
        Y
        finest_dx
    end
    
    methods
        function obj = ChomboOutputFlat(data, X, Y, finest_dx)
            obj.data = data;
            obj.X = X;
            obj.Y = Y;
            obj.finest_dx = finest_dx;
        end
        
        function c = plot(obj, var)

            dat = obj.dataForComp(var);
            %             [c,h] = contourf(obj.X, obj.Y, dat, 100);
            %             set(h, 'linestyle', 'none'); %remove contour lines
            %             colorbar();
            x = [obj.X(1, 1), obj.X(1, end)];
            y = [obj.Y(1,1), obj.Y(end, 1)];
            
            dx = obj.X(1,2) - obj.X(1,1);
            %dx = dx*4;
            %x(1) = x(1)  - dx/2;
            %x(end) = x(end) + dx/2;
            
            %y(1) = y(1)  - dx/2;
            %y(end) = y(end) + dx/2;
            
            
            im = imagesc(x, y, dat);
            c = colorbar();
            axis equal;
            
            % the dx/2's are vital for making things line up
            axis([x(1)-dx/2 x(2)+dx/2 y(1)-dx/2 y(2)+dx/2]);
            set(gca, 'YDir', 'normal');
          
        end
        
        function dataForComp = dataForComp(obj, comp)
            dataForComp = squeeze(obj.data(comp,:,:))';
        end
        
    end
    
end