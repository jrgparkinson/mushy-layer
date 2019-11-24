from mushyLayerRunUtils import write_inputs, get_executable_name, get_final_plot_file
import getopt
import sys
import os
import re
import numpy as np
from util import shared_storage

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# plt.style.use('classic')

def format_field_name(field):

    formatted_names = {'Porosity': 'Porosity, $\chi$',
                       'xDarcy velocity': '$x-$velocity',
                       'yDarcy velocity': '$y-$velocity',
                       'T err': '$\theta$'}

    if field in formatted_names.keys():
        return formatted_names[field]
    else:
        return field

def latexify(fig_width=None, fig_height=None, columns=1):
    """Set up matplotlib's RC params for LaTeX plotting.
    Call this before plotting a figure.

    Parameters
    ----------
    fig_width : float, optional, inches
    fig_height : float,  optional, inches
    columns : {1, 2}
    """

    # code adapted from http://www.scipy.org/Cookbook/Matplotlib/LaTeX_Examples

    # Width and max height in inches for IEEE journals taken from
    # computer.org/cms/Computer.org/Journal%20templates/transactions_art_guide.pdf

    assert (columns in [1, 2])

    if fig_width is None:
        fig_width = 3.39 if columns == 1 else 6.9  # width in inches

    if fig_height is None:
        golden_mean = (np.sqrt(5) - 1.0) / 2.0  # Aesthetic ratio
        fig_height = fig_width * golden_mean  # height in inches

    max_height_inches = 8.0
    if fig_height > max_height_inches:
        print("WARNING: fig_height too large:" + str(fig_height) +
              "so will reduce to" + str(max_height_inches) + "inches.")
        fig_height = max_height_inches

    # Need the mathsrfs package for \mathscr if text.usetex = True
    font_size = 9

    params = {'backend': 'ps',
              'text.latex.preamble': ['\\usepackage{gensymb}', '\\usepackage{mathrsfs}'],
              'axes.labelsize': font_size,  # fontsize for x and y labels (was 10)
              'axes.titlesize': font_size,
              'legend.fontsize': font_size,  # was 10
              'xtick.labelsize': font_size,
              'ytick.labelsize': font_size,
              'lines.markersize': 3,
              'lines.linewidth': 1,
              'text.usetex': True,
              'figure.figsize': [fig_width, fig_height],
              'font.family': 'serif'
              }

    matplotlib.rcParams.update(params)


def chombo_compare_analysis(data_folder):
    """
    Create inputs files and run them for doing chombo compare on all simulations in this directory
    """

    all_folders = [x for x in os.listdir(data_folder) if os.path.isdir(os.path.join(data_folder,x))]
    # print(all_folders)

    # Firstly compute richardson errors
    uniform_folders = [x for x in all_folders if 'Uniform' in x]

    # This doesn't sort properly
    uniform_folders = sorted(uniform_folders)

    uniform_resolution = []
    for folder in uniform_folders:
        res = get_folder_resolution(folder)
        uniform_resolution.append(res)

    # Sort by resolution
    uniform_folders = [x for _,x in sorted(zip(uniform_resolution, uniform_folders))]

    print(uniform_folders)

    for i in range(0, len(uniform_folders)-1):
        # for run_folder in uniform_folders:
        this_folder = os.path.join(data_folder, uniform_folders[i])
        next_folder = os.path.join(data_folder, uniform_folders[i+1])

        run_compare(next_folder, this_folder, 'richardson')

    # Now compute 512 errors
    fine_dir = os.path.join(data_folder, uniform_folders[-1])

    non_fine_folders = [x for x in all_folders if not x==uniform_folders[-1]]

    for folder in non_fine_folders:
        this_folder = os.path.join(data_folder, folder)

        run_compare(fine_dir, this_folder, 'finest')


def get_folder_resolution(folder):
    parts = re.findall('-(\d+)-', folder)
    res = int(parts[-1])
    return res


def get_folder_details(folder):

    if 'Uniform' in folder:
        coarse_nx = get_folder_resolution(folder)
        ref_rat = 0
        max_lev = 0
    else:
        # AMR-Subcycle-Reflux-Freestream0.99-MaxLevel1-ref4-PorousMushyHole-16--0
        # print(folder)
        result = re.findall('.*-MaxLevel(\d+)-.*ref(\d+)-.*-(\d+)-', folder)
        if result:
            parts = result[0]
            max_lev = int(parts[0])
            ref_rat = int(parts[1])
            coarse_nx = int(parts[2])
        else:
            max_lev = float('NaN')
            ref_rat = float('NaN')
            coarse_nx = float('NaN')


    return coarse_nx, ref_rat, max_lev


def run_compare(next_folder, this_folder, err_type):

    print('Running compare between this folder %s \n and next folder %s' % (this_folder, next_folder))

    # check this folder contains some valid files
    valid_files = [x for x in os.listdir(this_folder) if '.hdf5' in x]
    final_computed_plt_file = get_final_plot_file(this_folder)
    if not final_computed_plt_file:
        print('No plot files found in this folder')
        return

    final_exact_plt_file = get_final_plot_file(next_folder)
    if not final_exact_plt_file:
        print('No plot files found in higher resolution folder')
        return

    chombo_dir = os.environ['CHOMBO_HOME']
    compare_dir = os.path.join(chombo_dir, 'lib', 'util', 'ChomboCompare')
    compare_exec = get_executable_name(compare_dir, 'compare2d')
    compare_exec = os.path.join(compare_dir, compare_exec)
    # print('Found executable %s ' % compare_exec)

    this_err_folder = os.path.join(this_folder, 'error-' + err_type)

    if not os.path.exists(this_err_folder):
        os.makedirs(this_err_folder)

    compare_params_file = os.path.join(this_err_folder, 'compare-' + err_type + '.inputs')
    error_file = os.path.join(this_err_folder, 'err-' + err_type + '.2d.hdf5')

    computed_file = os.path.join(this_folder, final_computed_plt_file)
    exact_file = os.path.join(next_folder, final_exact_plt_file)
    compare_params = {'compare.sameSize': 0,
                      'compare.exactRoot': exact_file,
                      'compare.computedRoot': computed_file,
                      'compare.errorRoot': error_file,
                      'compare.doPlots': 1,
                      'compare.HOaverage': 0,
                      'compare.no_average_var': 'T err'}

    write_inputs(compare_params_file, compare_params)
    cmd = 'cd %s \n %s %s \n \n' % (this_err_folder, compare_exec, compare_params_file)

    print(cmd)

    os.system(cmd)

def load_error(folder):


    error_file = os.path.join(folder, 'pout.0')

    return load_error_file(error_file)

def load_error_file(error_file):
    errors = {}

    if not os.path.exists(error_file):
        # print('Cannot find %s' % error_file)
        return errors

    with open(error_file, 'r') as f:
        # file_contents = f.readlines()

        float_format = ['(-?\d\.\d+e[+-]\d+)']*4

        error_line_format = '^([\w\s]+):\s+' + ', '.join(float_format)

        # print(error_line_format)

        for line in f.readlines():

            #print(line)
            match = re.findall(error_line_format, line)

            if match:
                m = match[0]
                errors[m[0]] = {'L1': float(m[1]),
                                'L2': float(m[2]),
                                'Max': float(m[3]),
                                'Sum': float(m[4])}

            #print(errors)

    return errors


def run_chombo_compare(argv):

    figure_number = 8
    # data_folder = '/home/parkinsonjl/mnt/sharedStorage/TestDiffusiveTimescale/PorousMushyHole-t5e-05-hole0.04'
    # data_folder = '/home/parkinsonjl/mnt/sharedStorage/TestDiffusiveTimescale/PorousMushyHole-t5e-05-hole0.03'
    # data_folder = '/home/parkinsonjl/mnt/sharedStorage/TestDiffusiveTimescale/PorousMushyHole-t0.00015-hole0.04-veryGoodUseThis'
    #data_folder = '/home/parkinsonjl/mnt/sharedStorage/TestFinal/FixedPorousHole-1proc-minPorosity0.0-GOOD/'
    field = 'Porosity'
    err_type = 'L2'

    run_analysis = False
    include_richardson = True # for problems with no analytic solution




    figure_number = 6
    # data_folder = '/home/parkinsonjl/mnt/sharedStorage/TestDiffusiveTimescale/FixedPorousHole-1proc'
    data_folder = shared_storage.get_dir('TestFinal/FixedPorousHole-1proc-minPorosity0.0-GOOD/')
    run_analysis = False
    field = 'xDarcy velocity'
    err_type = 'L2'

    # figure_number = 0
    # data_folder = '/home/parkinsonjl/mnt/sharedStorage/TestDiffusiveTimescale/ConvectionDB-cfl0.17/chi0.4-Da1.0e-06-Ra1.0e+09'
    # run_analysis = False
    # field = 'Temperature'
    # err_type = 'L2'

    # data_folder = '/home/parkinsonjl/mnt/sharedStorage/TestDiffusiveTimescale/NoFlow/'
    # run_analysis = True
    # field = 'T err'
    # err_type = 'Max'
    # include_richardson = False

    try:
        opts, args = getopt.getopt(argv, "f:v:e:r:n:a")
    except getopt.GetoptError as err:
        print(str(err))
        print('run_chombo_compare.py -f <folder> -a<run analysis> -v <variable to consider> -e < err type> -r <include richardson errors?> -n <figure number>')
        sys.exit(2)

    for opt, arg in opts:
        if opt in "-f":
            data_folder = str(arg)
        elif opt in "-a":
            run_analysis = True
        elif opt in "-v":
            field = str(arg)
        elif opt in "-e":
            err_type = str(arg)
        elif opt in "-r":
            include_richardson = bool(arg)
        elif opt in "-n":
            figure_number = int(arg)


    # Compute the errors
    if run_analysis:
        chombo_compare_analysis(data_folder)

    # Collate errors and make plots

    # Each subfolder have these two folders:
    # richardson error is computed between consecutive simulations (128 vs 256, 256 vs 512)
    # fine error is always computed vs the finest resolution uniform mesh (64 vs 512, 128 vs 512 etc.)
    # richardson error only exists for uniform meshes
    richardson_error_folder = 'error-richardson'
    fine_error_folder = 'error-finest'

    # Get folders which have error calculations in them
    all_folders = [x for x in os.listdir(data_folder) if os.path.isdir(os.path.join(data_folder, x))
                   and (os.path.exists(os.path.join(data_folder, x, fine_error_folder)) or
                        os.path.exists(os.path.join(data_folder, x, richardson_error_folder)))]

    # Load all folders and put errors into a table indexed on resolution
    # resolution | richardson error | 512 uniform error | 512 amr errors...
    #  8
    #  16
    #  32
    #  ...


    err_data_sets = {}

    # Timings contains entries like [max lev, ref rat, coarse nx, time, ncells]
    timings = []

    richardson_name = 'Single-level Richardson'


    for folder in all_folders:
        this_richardson_err_folder = os.path.join(data_folder, folder, richardson_error_folder)
        this_fine_err_folder = os.path.join(data_folder, folder, fine_error_folder)

        coarse_nx, ref_rat, max_lev = get_folder_details(folder)

        # If we couldn't understand the folder name, skip it
        if np.isnan(coarse_nx):
            print('Skipping %s' % folder)
            continue


        richardson_err = load_error(this_richardson_err_folder)

        if max_lev == 0:

            # Also do richardson error here
            if include_richardson:

                if field not in richardson_err.keys():
                    print('Field %s not found for richardson error in folder %s' % (field, folder))
                    print('Available fields: ' + str(richardson_err.keys()))
                    continue

                if richardson_name in err_data_sets.keys():
                    err_data_sets[richardson_name].append([coarse_nx, richardson_err[field][err_type]])
                else:
                    err_data_sets[richardson_name] =  [[coarse_nx, richardson_err[field][err_type]]]


                data_set_name = 'Single-level 512 difference'
            else:
                data_set_name = 'Uniform'

        elif max_lev == 1:
            data_set_name  = '$n_{ref}$ = %d' % ref_rat
        else:
            data_set_name = '$n_{ref}$ = (%s)' % ','.join(['%d' % ref_rat]*max_lev)

        fine_err = load_error(this_fine_err_folder)

        if len(fine_err) == 0:
            continue

        if field not in fine_err.keys():
            print('Field %s not found for fine error in folder %s' % (field, folder))
            print('Available fields: ' + str(fine_err.keys()))
            continue


        this_err_entry = [coarse_nx, fine_err[field][err_type]]

        if data_set_name in err_data_sets.keys():
            err_data_sets[data_set_name].append(this_err_entry)
        else:
            err_data_sets[data_set_name] = [this_err_entry]


        # Also get and record timing details
        time_file = os.path.join(data_folder, folder, 'time.table.0')
        time = float('NaN')
        ncells = float('NaN')
        if os.path.exists(time_file):
            with open(time_file, 'r') as f:
                for line in f.readlines():
                    matches = re.findall('.*\[0\]main\s+(\d+\.\d+)\s+.*', line)

                    if matches:
                        time = float(matches[0])
                        break

        # [0]main 49.91810 1

        pout_file = os.path.join(data_folder, folder, 'pout.0')
        #total number of points updated = 819200

        with open(pout_file, 'r') as f:
            for line in f.readlines():
                matches = re.findall('.*total number of points updated = (\d+).*', line)

                if matches:
                    ncells = int(matches[0])
                    break

        these_timings = [max_lev, ref_rat, coarse_nx, time, ncells]

        timings.append(these_timings)

    print(err_data_sets)
    print(timings)

    ref_rats = set([x[1] for x in timings])
    # print(ref_rats)

    finest_timing = []
    for r in ref_rats:
        timings_for_ref_rat = [x for x in timings if x[1] == r]

        # Sort by resolution (2nd index)
        timings_for_ref_rat = sorted(timings_for_ref_rat, key=lambda x: x[2])

        # Take the last one
        finest_timing.append(timings_for_ref_rat[-1])

    print(finest_timing)


    latexify(fig_width=6.0, fig_height=2.5)

    # Make left axes wider
    fig, axes = plt.subplots(1, 2) #  gridspec_kw={'width_ratios':[2,1]}

    if include_richardson:
        key_order = ['Single-level Richardson', 'Single-level 512 difference']
    else:
        key_order = ['Uniform']

    key_order.extend(['$n_{ref}$ = 2', '$n_{ref}$ = 4', '$n_{ref}$ = (2,2)'])

    for ds_name in key_order:
        if ds_name not in err_data_sets.keys():
            continue

        ds = err_data_sets[ds_name]

        # Sort by nx
        ds = sorted(ds, key=lambda x:x[0])

        print(ds)

        nx = [x[0] for x in ds]
        err = [x[1] for x in ds]
        axes[0].plot(nx, err, marker='x', label=ds_name)


    # Also add 2nd order
    # Need to pick a data set to base this off

    collated_datasets = [value for values in err_data_sets.values() for value in values]
    # print('Collated datasets: ' + str(collated_datasets))
    all_nx = [x[0] for x in collated_datasets]
    min_nx = np.amin(all_nx)
    max_nx = np.amax(all_nx)*1.5

    all_err = [x[1] for x in collated_datasets]
    max_err = np.amax(all_err)
    init_err = max_err*4

    # print('all nx: ' + str(all_nx))

    #min_nx = np
    #a_ds_name = err_data_sets.keys()[0]
    #richardson_ds = err_data_sets[a_ds_name]
    nx_second_order = [min_nx, max_nx]
    err_second_order = [init_err, init_err*(float(min_nx)/float(max_nx))**2.0]

    print('nx 2nd order:' + str(nx_second_order))
    print('err 2nd order:' + str(err_second_order))

    # for i in range(1,len(nx_second_order)):
    #    err_second_order[i] = err_second_order[i-1] * (float(nx_second_order[i-1])/float(nx_second_order[i]))**2

    axes[0].plot(nx_second_order, err_second_order, linestyle=':', label ='2nd order')

    axes[0].set_xlabel('$1/\Delta x$')
    axes[0].set_ylabel('$L_2$ error (%s)' % format_field_name(field))

    axes[0].set_xscale('log')
    axes[0].set_yscale('log')

    # Make room for the legend
    #xl = axes[0].get_xlim()
    #yl = axes[0].get_ylim()
    #axes[0].set_xlim([xl[0], xl[1]*10])
    #axes[0].set_ylim([yl[0], yl[1] * 10])

    # Should sort out legend ordering
    # could add a title like title="(a) $\leftarrow$", if wanted
    leg_font_size = 8
    axes[0].legend(loc='center left', bbox_to_anchor=(1,0.75), prop={'size': leg_font_size})

    xl = axes[0].get_xlim()
    yl = axes[0].get_ylim()
    axes[0].text(xl[0]*0.9, yl[1]*1.3, '(a)')


    if finest_timing:
        ref_rats_plot = [x[1] for x in finest_timing]
        timings_plot = [x[3] for x in finest_timing]
        ncells_plot = np.array([float(x[4]) for x in finest_timing])

        largest_time = np.amax(timings_plot)
        largest_ncells = float(np.amax(ncells_plot))
        print('AMR performance normalisation. Max time = %.1f, Max num cells = %.2e' % (largest_time, largest_ncells))

        timings_plot = timings_plot/ largest_time
        ncells_plot = ncells_plot / largest_ncells

        # Make these black so they stand out from the other plot
        axes[1].plot(ref_rats_plot, timings_plot, marker='s', color='k', linestyle='-', label='Normalized CPU time')
        axes[1].plot(ref_rats_plot, ncells_plot, marker = 's', color='k', linestyle='--', label='Normalized cells advanced')

        axes[1].set_xlabel('Refinement ratio')
        # axes[1].set_ylabel('') # no y label

        axes[1].legend(loc='center right', bbox_to_anchor=(-0.2,0.25),  prop={'size': leg_font_size})

        axes[1].set_xlim([0, 4])
        axes[1].set_ylim([0, 1])

        axes[1].set_xticks([0, 2, 4])

        xl = axes[1].get_xlim()
        yl = axes[1].get_ylim()
        axes[1].text(-0.05, 1.03, '(b)')

    # set axis positions
    axes[0].set_position([0.1, 0.14, 0.3, 0.77])
    axes[1].set_position([0.77, 0.14, 0.2, 0.77])

    axes[0].tick_params(direction='in', which='both')
    axes[1].tick_params(direction='in', which='both')


    # Finally, save plot

    figure_output_directory = data_folder
    filename = 'Fig%dError-%s-%s.eps' % (figure_number, field, err_type)
    filename = filename.replace(' ', '_') # remove spaces
    figure_full_path = os.path.join(figure_output_directory, filename)
    print('Saving as %s' % figure_full_path)
    plt.savefig(figure_full_path, format='eps')

    plt.show()








if __name__ == "__main__":
    run_chombo_compare(sys.argv[1:])

