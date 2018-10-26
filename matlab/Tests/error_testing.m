res = [32,64,128,256];

% err = [0.2, 0.05, 0.01261, 0.003154]; %div(U^*) - 2nd order
% err = [0.01653, 0.004876, 0.001325, 0.0003470]; % U_x (2 levels) - 2nd order
%err = [0.3287, 0.131,0.05629, 0.02453]; % U.grad(theta) (with 2 levels) - 1st order
%err = [0.1418, 0.03709, 0.00939, 0.002360]; % U.grad(theta) (with 1 level) - 2nd order
%err = [0.1494, 0.03776, 0.01183, 0.005552]; % U.grad(theta) 2 levels with analytic U - 2nd order tending to 1st order at higher resolution
err = [0.01322, 0.003332, 0.0008347, 0.0002088]; %U.grad(theta) 1 level with analytic U fluxbox. 2nd order
%err = [0.03239, 0.01007, 0.004275, 0.001788]; %U.grad(theta) 2 levels with analytic U fluxbox. roughly 1st order
err = [0.01340, 0.003332, 0.0009712, 0.0003439]; %U.grad(theta) 2 levels with analytic U fluxbox and second order theta ghost cells - roughly second order

err = [0.02820, 0.007853, 0.002075, 0.0005350]; %U flux box from calculated U FArrayBox (2 levels) . 2nd order. 
err = [0.2286, 0.09517, 0.04038, 0.01784]; %U.grad(theta) 2 levels. Roughly 1st order

err = [0.05716, 0.02379, 0.01, 0.004460]; %U.grad(theta) 2 levels. Roughly 1st order

err = [0.0064, 0.001614, 0.0004]; %U.grad(theta) 2 levels enforced analytic grad(P). 2nd order.

err = [0.01675, 0.005124, 0.001338]; %U.grad(theta) 2 levels enforced P . 2nd order.
err = [0.02458, 0.01192, 0.006794]; %U.grad(theta) 2 levels enforced div(U^*) . 1st order.

figure(1);
plot(res, err)

figure(2);
plot(log(res), log(err))