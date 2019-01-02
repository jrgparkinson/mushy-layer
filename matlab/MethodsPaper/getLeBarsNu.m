% Nusselt number from benchmark problem
% given by Le Bars & Worster (2006)

function Nu = getLeBarsNu(chi, Da, Ra)

Nu = NaN;
if chi == 0.4
   if Da == 1e-6
       if Ra == 1e7
           Nu = 1.08; 
       elseif Ra == 1e8
           Nu = 3.07;
       elseif Ra == 1e9
           Nu = 12.9;
       end
           
   elseif Da == 1e-2
       if Ra == 1e3
           Nu = 1.01;
       elseif Ra == 1e4
           Nu = 1.41;
       elseif Ra == 1e5
           Nu = 3.17;
       elseif Ra == 5e5
           Nu = 5.24;
       end
           
   end
   
elseif chi == 0.9
    if Da == 1e-6
       if Ra == 1e7
           Nu = 1.08;
       elseif Ra == 1e8
           Nu = 3.08;
       elseif Ra == 1e9
           Nu = 13.15;
       end
   elseif Da == 1e-2
        if Ra == 1e3
           Nu = 1.02;
       elseif Ra == 1e4
           Nu = 1.67;
       elseif Ra == 1e5
           Nu = 4.09;
       elseif Ra == 5e5
           Nu = 6.89;
       end
   end
end


end