function [T, Sl, psi, a, Sa, x] = analyticSolChannel(Le, Ra, dSldz, T0, S0)

linearTEq = false;

if linearTEq

beta = Le*dSldz/6;


%Solve quadratic for a^2 (A1 l^2 + A2 l + A3 == 0 where l = a^2)
A1 = 6*beta^2*Ra^2*T0 + 3*beta^2*Ra - 6*beta^2*Ra^2;
A2 = beta*Ra*T0 + 9*beta - 3*S0*beta*Ra + 2*beta*Ra*T0;
A3 = S0-T0;
a2 = (-A2 - sqrt(A2^2 - 4*A1*A3))/(2*A1);
a = sqrt(a2);
%a = sqrt((Sa - S0)/(3*beta - 2*beta*Ra*(Sa-S0)));
%a1 = beta*Ra*


%a = 0.02;

Sa = (S0 + 3*beta*a^2 + 2*beta*Ra*T0*a^2)/(1+2*beta*Ra*a^2);
%Sa = T0 + (6*a^2*beta)/(3*a^2*beta*Ra-1)

% Integration constants (unknown)

alpha = (Sa - T0)/a;

x = linspace(0,1.0*a,400);

psi = -(1/6)*Ra*alpha*(x.^3-3*a^2*x);
T = T0 + 0.5*alpha*x.^2;
%Sl = S0 + Le*dSldz*( Ra*A*(a*x.^2 - (1/3)*x.^3) + Ra*A*0.1*x.^2 );
Sl = S0 + (1/24)*dSldz*Le*(Ra*alpha*x.^4 + 6*(2-a^2 * Ra * alpha)*x.^2 );
%Sl = S0 + alpha*(0.5*a^2*x.^2 - (1/12)*x.^4);


psi = -(1/2)*Ra*alpha*(x.^2-2*a*x);
T = T0 + alpha*x;
Sl = S0 + beta*(Ra*alpha*x.^3 + 3*(1 - a*Ra*alpha)*x.^2 );

else
    beta = Le*dSldz/24;

    %gamma = 1/a;
    gamma = 200;
    
    a = sqrt(24*beta + gamma)/sqrt(8*beta*gamma*Ra);
    a = 0.012;
    Sa = (-576 * beta^2 - 144*beta*gamma - 5*gamma^2 + 64*beta*gamma*Ra*S0)/ ...
        (64*beta*gamma*Ra);
    Sa = S0-0.1;
    
    x = linspace(0,1.0*a,400);
    
    T = T0 + 0.5*gamma*x.^2;
  
    psi = (1/6)*(Ra*gamma)*( 3*a^2*x - x.^3);
    Sl = S0 + beta*(Ra*gamma*x.^4 + 12*x.^2 - 6*Ra*gamma*a^2*x.^2);
    
    
end

end
