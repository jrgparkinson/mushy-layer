function [ phiNew, newDt ] = multigridSolver( phi, rhs, phi_const, x, nonlinearEq, opFunc, smoother, params )
%MULTIGRIDSOLVER Summary of this function goes here
%   Detailed explanation goes here
% Reference: Feng Wei Yang (thesis) September 2014 "Multigrid solution methods for nonlinear
% time-dependent system"

% Setup initial details
N_finest = size(phi) - 2;

rhs_finest = rhs;
phiConst_finest = phi_const;

%Initial guest
phiNew = rhs_finest;
resid = 10;
Vcycles = 0;

while resid > params.multigridMaxResidual
    
    prev_resid = resid;
    
    [phiNew, error, resid] = vcycle(phiNew, rhs_finest, phiConst_finest, ...
        N_finest, x, nonlinearEq, opFunc,  smoother, params);
    
    Vcycles = Vcycles+1;
    
    if Vcycles > params.maxNumCycle
        error = 1;
        fprintf('Max num vcycles reached \n');
    end
    
    if prev_resid/resid < 2
       error = 1;
       fprintf('Solver hang \n'); 
    end
    
    % Try halving dt if we get an error
    if error > 0
        params.dt = params.dt/2;
        % Go back to original guess
        phiNew = rhs_finest;
        resid = 10;
        Vcycles = 0; 
    end
    
    
    %fprintf('Vcycle complete, residual = %1.10f \n', resid);
    
end

newDt = params.dt;

end

function [phi_star, error, v_resid] = vcycle(phi_star, rhs, phi_const, N, x, nonlinearEq, opFunc, smoother, params)

error = 0;
phi_star_init = phi_star;
v_resid = 10;

% Do smoothing
for smooth = 1:params.numSmooth
    phi_star = smoother(phi_star, rhs, phi_const, nonlinearEq, params);
end

%Calculate residual
residual = rhs - opFunc(phi_star, phi_const,params);

%Restrict residual and fine grid soln
phi_star_coarse = restrict(phi_star, x);
residual_coarse = restrict(residual, x);
phiConst_coarse = restrict(phi_const, x);
x_coarse = restrict(x, x);
N = N/2;

%Compute modified RHS
%modified_rhs = rhs_coarse + residual_coarse - residualFunc(phi_star_coarse,rhs_coarse, phiConst_coarse,params);
%modified_rhs = rhs_coarse;
rhs_coarse = residual_coarse + opFunc(phi_star_coarse, phiConst_coarse,params);
%rhs_coarse = restrict(rhs, x);

% If fine grid
if N < 2
    %Solve exact problem
    residual = 10;
    numIter = 0;
    phi_guess = phi_star_coarse;
    while residual >  params.newtonResidLimit
        [ f, df ] = nonlinearEq( phi_guess, rhs_coarse, phiConst_coarse, params);
        phi_newguess = phi_guess - f./df;
        
        residual = max(phi_newguess - phi_guess);
        
        phi_guess = phi_newguess;
        
        numIter = numIter + 1;
        if numIter > params.maxIter
            disp('Nonlinear solve reached max iter');
            break;
        end
    end
    
    if residual < params.newtonResidLimit
        phi_improved = phi_guess;
        %error_coarse = Hguess - H_star_coarsest;
    else
        error = error + 1;
        return;
    end
    
else
    [phi_improved,error, ~] = vcycle(phi_star_coarse, rhs_coarse, phiConst_coarse, N, x_coarse, nonlinearEq, opFunc,  smoother, params);
end



%Compute the error
e = phi_improved - phi_star_coarse;

% interpolate error
e_fine = interpolate(e, x_coarse);

%Perform correction
phi_star= phi_star + e_fine;


% Do smoothing
for smooth = 1:params.numSmooth
    phi_star = smoother(phi_star, rhs, phi_const, nonlinearEq, params);
end


v_resid = max(abs(phi_star - phi_star_init));


end




function phi_fine = interpolate(phi_coarse, x_coarse)

N_coarse = length(x_coarse);
N_coarse_interior = N_coarse-2;

N_fine_interior = 2*N_coarse_interior;
N_fine = N_fine_interior + 2;

phi_fine = zeros(N_fine, 1);

for i =1:N_coarse_interior
    phi_fine(2*i-1) = phi_coarse(i);
    phi_fine(2*i) = 0.5*(phi_coarse(i)+phi_coarse(i+1));
end


end


% Calculate volume weighted average
function phi_coarse = restrict(phi_fine, x_fine)

N_fine = length(x_fine);
N_fine_interior = N_fine-2;

N_coarse_interior = N_fine_interior/2;
N_coarse = N_coarse_interior + 2;


phi_coarse = zeros(N_coarse, 1);

%Preserve BCs
phi_coarse(1) = phi_fine(1);
phi_coarse(end) = phi_fine(end);

for i=2:N_coarse-1
    phi_coarse(i) = (1/3)*( phi_fine(2*i) +  phi_fine(2*i-1) + phi_fine(2*i-2));
end

end
