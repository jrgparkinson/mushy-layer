function [z] = zMushWorster(CompRatio, thetaInf, Stefan, theta)

A = 0.5*(CompRatio + thetaInf + Stefan);
B = sqrt(A^2 - CompRatio*thetaInf);
alpha = A+B;
beta = A-B;
    
%Worster (1991) non-dimensional equation
z = ((alpha-CompRatio)/(alpha-beta)) * log((alpha+1)/(alpha-theta)) + ...
        ((CompRatio - beta)/(alpha-beta)) * log((beta+1)/(beta-theta));
       
end