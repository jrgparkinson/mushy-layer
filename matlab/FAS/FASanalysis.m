clear all; close all;

set(0, 'defaultlinelinewidth',3);
set(0, 'defaultaxeslinewidth',3);
set(0, 'defaultpatchlinewidth',3);
set(0, 'defaultAxesFontSize',22);
set(0, 'defaultTextFontSize',18);

Narr=[4,8,16,32,64,128,256,512];

maxErrs = [];
avErrs = [];
expectedErr = NaN*Narr;

for N_ii=1:length(Narr)
   N = Narr(N_ii);
   
   filename = ['solidification-N',num2str(N),'.mat'];
   load(filename);
   
   
    maxErrs(end+1) = maxErr;
    avErrs(end+1) = avErr;
    
    if N_ii == 2
        expectedErr(N_ii) = maxErr;
    elseif N_ii > 2
        expectedErr(N_ii) = expectedErr(N_ii-1)/4;
    end
    
    fprintf('Read %s \n', filename);
end



figure(1);
hold on;
box on; grid on;
plot(log(1./Narr), log(avErrs));
plot(log(1./Narr), log(expectedErr), '--');
xlabel('log($$\Delta x$$)', 'Interpreter','latex'); 

%ylabel('max$$(H-H_{analytic})$$', 'Interpreter','latex');
ylabel('log$$\sum (H-H_{analytic}) \Delta x$$', 'Interpreter','latex');

legend({'FAS','2nd order convergence'}, 'Location', 'southeast');


hold off;


print('FASconvergence', '-dpng', '-r300');


