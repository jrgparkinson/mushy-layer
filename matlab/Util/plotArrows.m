% Adds an arrow to each line of a contour plot
% Usage:
%   [C, h] = contour(X, Y, psi);
%   plotArrows (C, h, arrowColor, [minimum height]);
% 
% minimum height: if specified, minimum height is lowest vertical point where 
% arrows will be placed, as a fraction of the total height of the contour.
% This option is basically there to stop placing arrows near the domain
% boundary.
% 
%

function plotArrows(contour, handle, arrowColor, atHeight)

if nargin < 4
    atHeight = 0.0;
end

   lengthGrow = 0;
  
    headWidth = 10;
    headLength = 10;
   
   
    
    xcNew(1, :) = [NaN];
    ycNew(1, :) = [NaN];
    
    contour_i = 1;
    
    j = 1;
    for i=2:length(contour)
        if find(handle.TextList ==  contour(1,i))
            contour_i = contour_i+1;
            j = 1;
        else
            xcNew(contour_i, j) = contour(1,i);
            ycNew(contour_i, j) = contour(2,i);
            j=j+1;
        end
    end
    
    numContours = contour_i;
    
    for contour_i =1:numContours
        
        xc = squeeze(xcNew(contour_i, :));
        yc = squeeze(ycNew(contour_i, :));
        
        yc(yc==0) = NaN;
        xc(xc==0) = NaN;
    
        plottedOnLine = false;
        
        for i=3:(length(xc)-1)
            
            percentHeight = (yc(i)-min(yc))/(max(yc)-min(yc));
            if (percentHeight > atHeight ...
                    && ~plottedOnLine)
                if xc(i) < xc(i+1)
                    xinit = xc(i) - lengthGrow;
                    xfinal = xc(i+1) + lengthGrow;
                else
                    xinit = xc(i+1) - lengthGrow;
                    xfinal = xc(i) +lengthGrow;
                end
                
                % Do this to get arrow data
                hq = quiver(xinit,yc(i),(xfinal-xinit),(yc(i+1)-yc(i)),0);
                
                hq.Visible = 'off';
                
                % Get arrow data
                Uquiv = hq.UData;
                Vquiv = hq.VData;
                Xquiv = hq.XData;
                Yquiv = hq.YData;
                
                xSign = 1;
                ySign = 1;
                
                if xc(i+1) < xc(i)
                    xSign = -1;
                end
                
                if yc(i+1) < yc(i)
                    ySign = -1;
                end
                   
                
                ah = annotation('arrow',...
                    'headStyle','cback1',...
                    'HeadLength',headLength,'HeadWidth',headWidth, ...
                    'Color', arrowColor, 'LineStyle', 'none');
                set(ah,'parent',gca);
                
                xStart = Xquiv;
                
                if xSign == -1
                    xStart = Xquiv + Uquiv;
                end
                
                set(ah,'position',[xStart Yquiv xSign*Uquiv ySign*Vquiv]);

                plottedOnLine = true;
            end
            
            
            
        end
    end
    
end