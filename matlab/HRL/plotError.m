function plotError(error, final_plot, plot_title, flags_list, grid_res, use_Nu_caltagirone)
    if nargin < 6
        use_Nu_caltagirone = false;
    end
    
    NUSSELT = -1; NUSSELT_ERROR = -2; 
    num_flags = length(flags_list);

    hold on;
    for i=1:num_flags
        if final_plot == NUSSELT
            plot( log(grid_res), error(i, :), 'x-');
        elseif final_plot == NUSSELT_ERROR
            plot( log(grid_res), log(error(i, :)), 'x-');
        elseif final_plot > 0
            plot( log10(grid_res), log10(error(i, :)), 'x-');
        end
    end
    xlabel('log(N_x)');
    %axis([2 6 -7 -1]);
    if final_plot == NUSSELT
        title('nusselt number');
        ylabel('Nu');
    elseif final_plot == NUSSELT_ERROR
        title('nusselt error convergence');
        if use_Nu_caltagirone
            ylabel('log(Nu - Nu_{caltagirone})');
        else
            ylabel('log(Nu - Nu_{finest 1 level})');
        end
    elseif final_plot > 0
        title('theta error convergence');
        err_name = ChomboCompare.err_title(final_plot);
        ylabel([err_name{1}, '(theta - theta_{finest 1 level})']);
    end


    leg = [];
    keySet = keys(plot_title);
    for i = 1:num_flags
        leg{1, i} = plot_title(flags_list{1,i});
    end

    legend(leg);
    %axis equal;
    grid on;
    hold off;

end