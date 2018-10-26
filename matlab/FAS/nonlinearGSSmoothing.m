function phi = nonlinearGSSmoothing(phi, rhs, phi_const, nonlinearEq, params)

for i=2:length(phi)-1
    % Just calculate f and df at the point we're interested in
    [ f, df ] = nonlinearEq( phi(i-1:i+1), rhs(i-1:i+1), phi_const(i-1:i+1), params);
    phi(i) = phi(i) - f(2)/df(2);
end

end

