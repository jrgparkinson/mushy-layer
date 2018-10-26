function [z] = zMush(theta, params)


if (params.stefan ==  0)
    z = 0;
    return;
end

    A = 0.5*(params.concRatio + params.thetaInf + params.stefan);
    B = sqrt(A^2 - params.concRatio*params.thetaInf - params.stefan*params.thetaInterface);
    alpha = A+B;
    beta = A-B;
    
    
    z = (1/params.V) *( ((alpha-params.concRatio)/(alpha-beta)) * log((alpha)/(alpha-theta)) + ...
        ((params.concRatio - beta)/(alpha-beta)) * log((beta)/(beta-theta)) );
    
   
    
end