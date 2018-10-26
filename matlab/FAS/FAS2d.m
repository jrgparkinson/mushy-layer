function  FAS2d
close all;

params.gamma = 1000.0;
params.numSmooth = 4;
params.num_ghost = 1;
params.alpha = -1;

[X,Y] = getGrids(16, params.num_ghost);
params.dx =  X(1, 2)-X(1, 1);

rhs_with_ghosts = -params.alpha*2*((X-X.*X) + (Y-Y.*Y)) + params.gamma*((X-X.^2).*(Y-Y.^2)).*exp((X-X.^2).*(Y-Y.^2));
rhs = rhs_with_ghosts(2:end-1, 2:end-1);


u = 0*X;
exactSol = (X-X.*X).*(Y-Y.*Y);

%u = exactSol;

h = figure();
set(h, 'Position', [200 200 1000 700]);

initres = residual(u, rhs, params, params.dx);
fprintf('Initial residual norm: %1.2e \n', max(max(abs(initres))) );


for i=1:10
    uold = u;
    u_corr = cycle(u, rhs, params);
    %u = relax(u, rhs, params);
    u = uold + u_corr;
    
    
    res = residual(u, rhs, params, params.dx);
    sumRes = sum(sum(abs(res)));
    L1 = sumRes/(size(res,1)*size(res,2));
    fprintf('Iteration %d, residual norm: %1.2e \n', i, max(max(abs(res))) );
    
    if L1 < 1e-10
        fprintf('Converged \n');
        break;
    end
    m=2; n=2;
    
    subplot(m, n, 1);
    %pcolor(X, Y, u);
    imagesc(u(2:end-1, 2:end-1));
    title(['Iter: ', num2str(i)]);
    colorbar();
    
    subplot(m, n, 2);
    %pcolor(X, Y, exactSol);
    imagesc(exactSol(2:end-1, 2:end-1));
    title(['Exact']);
    colorbar();    
    
    subplot(m, n, 3);
    %pcolor(X, Y, u-exactSol);
    %imagesc(u(2:end-1, 2:end-1)-exactSol(2:end-1, 2:end-1));
    imagesc(res(2:end-1, 2:end-1));
    title(['Residual']);
    colorbar();
    
    subplot(m, n, 4);
    %pcolor(X, Y, f);
    imagesc(rhs(2:end-1, 2:end-1));
    title(['RHS']);
    colorbar();
    
    drawnow;
    
    temp=0;
    
end



end




function u_corr = cycle(u, rhs, params)

%fprintf('Cycle on grid %d x %d \n', size(u));
N_mesh = size(u,1)-2;
[X,Y] = getGrids(N_mesh, params.num_ghost);
dx = X(1, 2)-X(1, 1);

ng = params.num_ghost;

% first check we're not on the bottom level
% remember we have two ghost cells in each direction
if min(size(u)) < 5
    uold = u;
    
    u = bottomSolve(u, rhs, params, dx);
    
    
    u_corr = u - uold;
else
    
    % first relax
    uold = u;
    
    u = relax(u, rhs, params, dx);
    
    % compute fine grid residual
    resid = residual(u,rhs,params,dx);
    
    % setup for coarser solve(s)
    u_coarse = restrict(u, 1);
    %residual_coarse = restrict(resid(ng:end-ng, ng:end-ng), 0);
    residual_coarse_ghost = restrict(resid, 1);
    residual_coarse = residual_coarse_ghost(1+ng:end-ng, 1+ng:end-ng);
    
    %wrong_rhs_coarse = restrict(rhs);
    
    % rhs on coarse level is NOT restrict(rhs)
    % but is restrict(residual + op(coarse solution))
    op_u_coarse = applyOp(u_coarse, params);
    rhs_coarse = residual_coarse +  ...
        op_u_coarse(1+ng:end-ng, 1+ng:end-ng);
    
    % do coarser solve(s)
    u_corr_coarse = cycle(u_coarse, rhs_coarse, params);
    
    % get back solution at this resolution
    %    diff = u-prolong(u_coarse);
    %    sumDiff = sum(sum(diff));
    
    u_corr_coarse = fillBCs(u_corr_coarse);
    u = u + prolong(u_corr_coarse, ng);
    
    % final relaxation
    % params.numSmooth = 4;
    u = relax(u, rhs, params,  dx);
    
    u_corr = u-uold;
    
    temp = 0;
    %params.numSmooth = 4;
end

end


function u_coarse = restrict(u_fine, num_ghost)

% fineInterior = u_fine(2:end-1, 2:end-1);
%
% coarseInterior = resizem(fineInterior, 0.5);
%
% coarseSize = size(coarseInterior)+2;
%
% u_coarse = NaN*ones(coarseSize);
% u_coarse(2:end-1, 2:end-1) = coarseInterior;


Nfine_interior = size(u_fine, 1) - 2*num_ghost ;
Ncoarse_interior = round(Nfine_interior/2);
% [Xcoarse,Ycoarse] = meshgrid(linspace(0,1,Ncoarse), linspace(0,1,Ncoarse));
%    [Xfine,Yfine] = meshgrid(linspace(0,1,Nfine), linspace(0,1,Nfine));
%
%

[Xfine,Yfine] = getGrids(Nfine_interior, num_ghost);
[Xcoarse,Ycoarse]= getGrids(Ncoarse_interior,  num_ghost);
u_coarse = interp2(Xfine,Yfine,u_fine,Xcoarse,Ycoarse);

%u_coarse = fillBCs(u_coarse);

end



function u = fillBCs(u)

% BCs are 0 on all boundaries


if (size(u, 1) > 4)
    % 2nd order
    a = -2; b = 1/3;
 u(1, :) = a*u(2, :) + b*u(3, :);
 u(:, 1) = a*u(:, 2) + b*u(:, 3);
 u(end, :) = a*u(end-1, :) + b*u(end-2, :);
 u(:, end) = a*u(:, end-1) + b*u(:, end-2);
else
% 1st order
u(1, :) = -u(2, :);
u(:, 1) = -u(:, 2);
u(end, :) = -u(end-1, :);
u(:, end) = -u(:, end-1);
end

u(1, :) = -u(2, :);
u(:, 1) = -u(:, 2);
u(end, :) = -u(end-1, :);
u(:, end) = -u(:, end-1);

% 0th order

% u(1, :) = 0;
% u(:, 1) = 0;
% u(end, :) = 0;
% u(:, end) = 0;

end

function u_fine = prolong(u_coarse, ng)

%coarseInterior = u_coarse(1+ng:end-ng, 1+ng:end-ng);

Ncoarse_interior = size(u_coarse, 1)-2*ng;
Nfine_interior = Ncoarse_interior*2;

[Xfine,Yfine] = getGrids(Nfine_interior, ng);
[Xcoarse,Ycoarse]= getGrids(Ncoarse_interior, ng);

u_fine = interp2(Xcoarse,Ycoarse,u_coarse,Xfine,Yfine, 'cubic');

end

function u = relax(u, rhs, params, dx)


for iter=1:params.numSmooth
    
    u =fillBCs(u);
    
    for i=2:size(u, 1) -1
        
        for j=2:size(u, 2)-1
            lap = params.alpha*(u(i+1, j) + u(i-1, j) + u(i, j+1) + u(i, j-1) - 4*u(i,j))/(dx*dx);
            nl = params.gamma*u(i,j)*exp(u(i,j));
            res = lap+nl-rhs(i-params.num_ghost,j-params.num_ghost);
            
            denom = -params.alpha*4/(dx*dx) + params.gamma*(1+u(i,j))*exp(u(i,j));
            
            u(i,j) = u(i,j) - res/denom;
            %u(i,j) = u(i,j) + 0.25*(lap*dx*dx - dx*dx*f(i,j));
            
        end
    end
end


end

function r= residual(u, rhs, params, dx)
r = 0*u;

nonlinearterm = 0*u;
lapterm  = 0*u;

u =fillBCs(u);

for i=2:size(u, 1) -1
    for j=2:size(u, 2)-1
        lap = params.alpha*(u(i+1, j) + u(i-1, j) + u(i, j+1) + u(i, j-1) - 4*u(i,j))/(dx*dx);
        nl = params.gamma*u(i,j)*exp(u(i,j));
        r(i,j) = rhs(i-params.num_ghost,j-params.num_ghost) - (lap+nl);
        
        
        nonlinearterm(i,j) = nl;
        lapterm(i,j) = lap;
        % denom = -4/(dx*dx) + params.gamma*(1+u(i,j))*exp(u(i,j));
        
        % u(i,j) = u(i,j) - res/denom;
        %u(i,j) = u(i,j) + 0.25*(lap*dx*dx - dx*dx*f(i,j));
        
    end
end

temp = 0;

end

function  op = applyOp(u, params)

u =fillBCs(u);

Ncoarse = size(u, 1) - 2;
[X,Y] = getGrids(Ncoarse, params.num_ghost);
dx = X(1, 2) - X(1,1);

op = 0*u;

for i=2:size(u, 1) -1
    for j=2:size(u, 2)-1
        lap = params.alpha*(u(i+1, j) + u(i-1, j) + u(i, j+1) + u(i, j-1) - 4*u(i,j))/(dx*dx);
        nl = params.gamma*u(i,j)*exp(u(i,j));
        op(i,j) = lap+nl;
        
        % denom = -4/(dx*dx) + params.gamma*(1+u(i,j))*exp(u(i,j));
        
        % u(i,j) = u(i,j) - res/denom;
        %u(i,j) = u(i,j) + 0.25*(lap*dx*dx - dx*dx*f(i,j));
        
    end
end

end

function u = bottomSolve(u, rhs, params, dx)
u = relax(u, rhs, params, dx);
end

function [X,Y] = getGrids(Ninterior,  num_ghost)
dx = 1/Ninterior;

N = Ninterior+2*num_ghost;

minVal = dx/2-num_ghost*dx;
maxVal = 1-dx/2 + num_ghost*dx;

x = linspace(minVal,maxVal, N);
y = linspace(minVal,maxVal, N);

%x = linspace(0, 1, N);
%y = linspace(0, 1, N);

[X, Y] = meshgrid(x,y);

dxNew = x(2)-x(1);
% if dxNew ~= dx
%     fprintf('Error generating grids!');
% end

end

function myPlot(field)

figure();
imagesc(field(2:end-1, 2:end-1));
c = colorbar();
c.Ticks =  c.Limits;

end