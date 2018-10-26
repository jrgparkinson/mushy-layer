
function A = saltEqOp(Sstar, H, params)

A=0*Sstar;

Sl = liquidSalinity(H, Sstar, params);
chi = porosity(H, Sstar, params);

for i=2:length(Sstar)-1

chiPlusHalf = 0.5*( chi(i)+chi(i+1) );
chiMinusHalf = 0.5*( chi(i)+chi(i-1) );

A(i) = Sstar(i) - params.dt* ( (params.V/(2*params.dx))*(Sstar(i-1)-Sstar(i+1)) ... 
    + (1/(params.Le*params.dx^2))*( chiPlusHalf*(Sl(i+1)-Sl(i)) -  chiMinusHalf*(Sl(i)-Sl(i-1)) )    );

end

end