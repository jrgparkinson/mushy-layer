function [theta, porosity] = analyticSoln(zGrid, params)
%Grids
theta = zGrid;
porosity = zGrid;
zMushCalc = [];

% Find the top BC and height of mushy layer
lower = params.thetaBottom;
upper = 2*lower;

% initial guess
params.thetaInterface = 1; % Strictly: ( params.lewis * params.ThetaInf - params.thetaInf) / (params.lewis -1 );

tolerance = 1e-8;

thetaBottomGuess = 0;

while((params.thetaBottom-thetaBottomGuess)^2 > tolerance)
    
  
    params.thetaInf =  (upper+lower)/2;
    hGuess = zMush(params.thetaInterface, params);
    thetaBottomGuess = params.thetaInf + 	(1-params.thetaInf)*exp(-params.V*(1-hGuess));
    
    % If we guessed too high, reduce thetaInf
    if (thetaBottomGuess > params.thetaBottom)
        
        upper = params.thetaInf;
        
    else
        
        % We guessed too low, increase thetaInf
        lower = params.thetaInf;
        
    end
end

params.thetaInterface = ( params.lewis * params.ThetaInf - params.thetaInf) / (params.lewis -1 );
h = zMush(params.thetaInterface, params);


thetaDefined = [];
thetaGrid = params.thetaTop:0.0005:params.thetaBottom;
zMushCalc = [];

%Find z(theta) for theta < thetaInterface
for T_i = 1:length(thetaGrid)
    T = thetaGrid(T_i);
    z = zMush(T, params);
    
    if z < h && isreal(z)
        zMushCalc(end+1) = zMush(T, params);
        thetaDefined(end+1) = T;
    end
end


%Combine solutions in two different regions
for z_i =1:length(zGrid)
    z = zGrid(z_i);
    
     if params.thetaBottom > params.thetaTop
         z = 1-z;
     end
    
    if z<=h
        %Interpolate
        %thetaTotal[z_i] =
        theta(z_i) = interp1(zMushCalc,thetaDefined,z, 'spline');
        
       
        porosity(z_i) = 1 - (1-theta(z_i))/(params.concRatio - theta(z_i));
    else
        theta(z_i) = params.thetaInf + (params.thetaInterface - params.thetaInf)*exp(-params.V * (z-h));
        
        porosity(z_i) = 1;
    end
    
   

    
end

  
    
end

