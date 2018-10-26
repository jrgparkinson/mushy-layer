H = 0.2;
baseL = 0.2;
% Must be even
Nx = 164; Ny = 128;

L = linspace(0.8, 1.2, 9)*baseL;

for i=1:(40)
   newNx = Nx - 2*i;
    %for j = 1:(Ny/2)
        newNy = Ny; % + j*2;
        
        newL = baseL*newNx/newNy;
        
        foundL = false;
        
        
           fprintf('(Nx, Ny) = (%d, %d), L = %f \n', newNx, newNy, newL); 
           
     
    
end