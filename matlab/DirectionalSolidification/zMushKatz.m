function [zKatz] = zMushKatz(CompRatio, thetaInf, Stefan, diffusivity, nonDimV, height, theta)

A = 0.5*(CompRatio + thetaInf + Stefan);
B = sqrt(A^2 - CompRatio*thetaInf - Stefan);
alpha = A+B;
beta = A-B;
    

zKatz =(1/nonDimV)* ( ((alpha-CompRatio)/(alpha-beta)) * log((alpha)/(alpha-theta)) + ...
                      ((CompRatio - beta)/(alpha-beta)) * log((beta)/(beta-theta)) );
         
end