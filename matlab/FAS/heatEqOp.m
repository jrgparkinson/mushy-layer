
function A = heatEqOp(Hstar, S, params)

A=0*Hstar;

T = temperature(Hstar, S, params);

for i=2:length(Hstar)-1
   
A(i) = Hstar(i) - params.dt* ( (params.V/(2*params.dx))*(Hstar(i-1)-Hstar(i+1)) ... 
    + (1/params.dx^2)*(T(i+1)+T(i-1)-2*T(i))  );

end

end