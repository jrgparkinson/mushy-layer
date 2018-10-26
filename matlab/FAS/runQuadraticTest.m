%clear all;
close all;


plotInterval = 10;
enforceAnalyticSoln = false;

% The current implementation of Jacobi smoothing is much quicker
smoothing = @nonlinearJacobiSmoothing;
%smoothing = @nonlinearGSSmoothing;

Nvals = [32];

for N_i=1:length(Nvals)
    
    
    params = struct;
    
    % Physical constants
    params.stefan = 5.7;
    params.pc = 0; %Should be zero for this problem
    params.cp = 1;
    params.concRatio = 1.4;
    params.Te = 0;
    params.Se = 1;
    params.V = 1;
    params.Le = 30;
    
     % Grid and time stepping
    params.dt = 1e-3;
    params.dtGrowth = 1.0; %1 = no growth, greater than 1 means frow at each timestep
    params.N = Nvals(N_i); %Num interior points
    
    params.t_max = 1e5*params.dt;
    params.x_max = 1; params.x_min = 0;
    params.dx = (params.x_max-params.x_min)/(params.N+1);
    t = 0:params.dt:params.t_max;
    x = params.x_min:params.dx:params.x_max;
    
    params.steadyStateLimit = 1;
    
    %Before we start, get the analytic solution
    params.thetaInterface = 1; % This is determined by our nondimensionalisation
    params.thetaInf = 1.4;
    [h, thetaAnalytic, chiAnalytic] = directionSolidificationAnalyticSol(x, params);
    
    % Flip analytic soln
    %thetaAnalytic = fliplr(thetaAnalytic);
    %chiAnalytic = fliplr(chiAnalytic);
    
    %params.thetaTop =params.thetaInf + (params.thetaInterface-params.thetaInf)*exp(-params.V*(1-h));
    params.thetaBottom = 1.39; 
    
    % Initial and boundary conditions
    params.S0 = 1;
    params.BC_S = [params.S0 1];
    
    params.H0 = params.thetaBottom + params.stefan;
    porosity_eutectic = 1-params.Se/params.concRatio;
    params.He = 0 + params.stefan*porosity_eutectic;
    params.BC_H = [params.H0 params.He]; %Format is x_lo, x_hi, z_lo, z_hi etc...
    
    
   
    
    % Nonlinear solver
    params.maxIter = 10;
    params.maxNumCycle = 70;
    params.numSmooth = 8;
    params.newtonResidLimit = 1e-10;
    params.multigridMaxResidual = 1e-10;
    
  
    % Make initial enthalpy
    H = params.H0*ones(length(x), length(t));
    S = params.S0*ones(length(x), length(t));
    
    %Enforce analytic solution
    if enforceAnalyticSoln
        H(:, 1) = thetaAnalytic + params.stefan*chiAnalytic;
    end
    
    %Enforce BCs
    H(1, :) = params.BC_H(1);
    H(end, :) = params.BC_H(2);
    
    S(1, :) = params.BC_S(1);
    S(end, :) = params.BC_S(2);
    
    
    
    
    
    figure(1);
    %hold on;
    
    for t_i = 1:(length(t)-1)
        
        
        
        %Solve multigrid for H^n+1
        %H(:, t_i+1) = multigridSolve(H(:, t_i), @operator, params);
        
        % Linear bottom solve
        %A = nonlinearDiffusionOperator(H(:, t_i), S(:, t_i), params);
        %H(:, t_i+1) = A\H(:, t_i);
        
        % Nonlinear bottom solve
        %     residual = 1;
        %     Hguess = H(:, t_i);
        %     while residual > 1e-5
        %         [ f, df ] = heatEqNonLinear( Hguess, H(:, t_i), S, params);
        %         Hnewguess = Hguess - f./df;
        %
        %         residual = max(Hnewguess - Hguess);
        %
        %         Hguess = Hnewguess;
        %     end
        %     H(:, t_i+1) = Hguess;
        
        
        % Do first salt update
        [S(:, t_i+1), params.dt] = multigridSolver( S(:, t_i), S(:, t_i), H(:, t_i), x, ...
            @saltEqNonLinear, @saltEqOp, @nonlinearJacobiSmoothing,  params );
        
        Sdiff = 1;
        niter = 0;
        
        while Sdiff > 1e-5
            
            % Nonlinear multigrid solve using updated S
            [H(:, t_i+1), params.dt] = multigridSolver( H(:, t_i), H(:, t_i), S(:, t_i+1), x, ...
                @heatEqNonLinear, @heatEqOp, smoothing, params );
            
            %Now do S solve again with new H and check it's fairly similar to before
            S_old = S(:,t_i+1);
            
            [S(:,t_i+1), params.dt] = multigridSolver( S(:, t_i), S(:, t_i), H(:, t_i+1), x, ...
                @saltEqNonLinear, @saltEqOp, smoothing, params );
            
            Sdiff = max(abs(S(:,t_i+1)-S_old));
            niter = niter + 1;
        end
        
        
        
        if mod(t_i-1,plotInterval) == 0 && plotInterval > 0
            chi = porosity( H(:, t_i+1),  S(:, t_i+1), params);
            theta = temperature( H(:, t_i+1),  S(:, t_i+1), params);
            Sl = liquidSalinity( H(:, t_i+1),  S(:, t_i+1), params);
            
            clf;
            hold on;
            plot(  x, theta);
            plot(  x, chi, '--');
            plot(  x, S(:, t_i+1), ':');
            
            plot(x, thetaAnalytic, '-');
            plot(x, chiAnalytic, '--');
            
            
            hold off;
            
            xlabel('x');
            ylabel('value');
            title(['t = ',num2str(t(t_i))]);
            
            legend({'Temperature','porosity','bulk salinity'});
            axis([0.5 1 0 params.thetaInf]);
            
            pause(0.01);
        end
        
        
        dHdt = max(abs(H(:, t_i+1) - H(:, t_i)))/params.dt;
        dSdt = max(abs(S(:, t_i+1) - S(:, t_i)))/params.dt;
        
        if mod(t_i-1,10) ==0
            fprintf('Step %d, N outer iterations = %d,  dH/dt = %1.3f, dSdt = %1.3f, new dt = %1.8f \n', t_i, niter, dHdt, dSdt, params.dt);
        end
        
        % Try increasing dt slightly
        params.dt = params.dt*params.dtGrowth;
        
        if (dHdt < params.steadyStateLimit && dSdt <  params.steadyStateLimit)
            disp('Steady state reached');
            break;
        end
        
        
        
    end
    
    theta = temperature( H(:, t_i+1),  S(:, t_i+1), params);
    
    err = theta - thetaAnalytic.';
    maxErr = max(abs(err));
    avErr = sum(abs(err))/params.N;
    % After reaching steady state, calculate the error and save it somewhere
    save(['solidification-N',num2str(params.N),'.mat']);
    
    fprintf('N=%d, maxErr = %1.6f, avErr = %1.6f \n',params.N,maxErr,avErr);
    
    
    
end

