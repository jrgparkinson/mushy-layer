% Nonlinear Jacobi for now as don't need to update f and df after each cell
% update
function phi = nonlinearJacobiSmoothing(phi, rhs, phi_const, nonlinearEq, params)

[ f, df ] = nonlinearEq( phi, rhs, phi_const, params);

%H = H - f./df;

for i=2:length(phi)-1
    phi(i) = phi(i) - f(i)/df(i);
end

end