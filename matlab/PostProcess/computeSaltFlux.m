% Compute the vertical salt flux =(1/Le) chi*dSl/dz - U_z*S_l - V*S
function saltFlux = computeSaltFlux(S, Uz, Sl, chi, dz, V, Le)

[~, dSldz] = gradient(Sl, dz);

F = (1/Le)*chi.*dSldz - Uz.*Sl - V*S ; 

HorizAv = mean(F.');  % Averaged over x direction

saltFlux = mean(HorizAv);
%saltFlux = median(HorizAv);

temp = 0;

end