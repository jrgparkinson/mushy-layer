%Load all output from Ra250 run and see if Nu(t) is converging
Ra = 50;
plot_prefix = ['bm2-Ra', num2str(Ra), '-32pts'];
output_dir = '/Users/parkinsonj/Documents/OneDrive/Documents/Oxford/Data/Ra50-32-0/';
 

Nu_calt =  containers.Map({100,200,250}, [2.651, 3.813 4.199]);
 % Value from Caltagirone (1975)
times = 0:100:17100;

Nu_t_av = [];
Nu_t_base = [];
Nu_t_calt = [];

for t_i=1:length(times)

    Frame = times(t_i);
    output = MushyLayerOutput(dim, Frame, output_dir, plot_prefix);
  
    [Nu_t_av(t_i), Nu_t_base(t_i)] = output.nusseltNumber();
    Nu_t_calt(t_i) = Nu_calt(Ra);
end

figure(1);
plot(times, Nu_t_calt, times, Nu_t_av, times, Nu_t_base);
legend('Nu Caltagirone 1975', 'Nu(t) average','Nu(t) base');
title('Convergence of Nu(t)');