# Run all convergence tests for the methods paper
import os
import sys
from colorama import Fore, Style
import getopt
from runAMRConvergenceTest import runTest
from BatchJob import BatchJob
from mushyLayerRunUtils import get_base_output_dir, get_matlab_base_command, read_inputs, get_mushy_layer_dir, check_exec_exists


######################################
# 2) Convection in a fixed porous medium
# Just run for a single grid resolution and compare to previously published values
######################################

def uniform_porous_resolution_specific_params(nz_coarse, ref_rat, max_level, max_refinement):
    mushyLayerBaseDir = get_mushy_layer_dir()

    nx_coarse = nz_coarse

    gridFile = mushyLayerBaseDir + '/grids/leftRight/' + str(nx_coarse) + 'x' + str(nz_coarse)

    params_file = mushyLayerBaseDir + '/params/convergenceTest/convectionDarcyBrinkmanConvTest.parameters'
    params = read_inputs(params_file)

    params['main.refine_thresh'] = str(3.0 / float(nz_coarse))
    params['main.tag_buffer_size'] = str(max(2, int(float(nz_coarse) / 8)))

    # Make sure we don't split up the grids as there's currently a bug in Chombo with higher order advection methods
    params['main.max_grid_size'] = int(nz_coarse*2)

    return nx_coarse, params, gridFile

def get_default_cfl():
    cfl = 0.1
    return cfl

def test_uniform_porous_convection(argv):
    # Default vals:
    cfl = get_default_cfl()
    nz_uniform = 64 # don't do uniform mesh simulations by default
    nz_vm = 64
    chi = 0.4
    bc_accuracy = 2
    u_del_u_method = 0
    advection_method = 1

    fixed_da = -1.0
    fixed_ra = -1.0

    try:
        opts, args = getopt.getopt(argv, "n:v:c:p:b:u:a:d:r:")
    except getopt.GetoptError as err:
        print(str(err))
        print('test_uniform_porous_convection.py -n<num uniform mesh grid points> -v<num variable mesh points>'
              ' -c<cfl> -p <chi> -b <BC accuracy> -u <u del u method> -a <advection method> -d<fixed darcy> -r<fixed Ra>')
        sys.exit(2)

    for opt, arg in opts:
        if opt in "-n":
            nz_uniform = int(arg)
        elif opt in "-c":
            cfl = float(arg)
        elif opt in "-v":
            nz_vm = int(arg)
        elif opt in "-p":
            chi = float(arg)
        elif opt in "-b":
            bc_accuracy = int(arg)
        elif opt in "-u":
            u_del_u_method = int(arg)
        elif opt in "-a":
            advection_method = int(arg)
        elif opt in "-d":
            fixed_da = float(arg)
        elif opt in "-r":
            fixed_ra = float(arg)

    base_output_dir = get_base_output_dir()
    matlab_command = get_matlab_base_command()

    # Check we have the executables needed to run this test
    success = check_exec_exists(os.path.join(get_mushy_layer_dir(), 'execSubcycle'), 'mushyLayer')
    success = success and check_exec_exists(os.path.join(get_mushy_layer_dir(), 'setupNewRun'), 'setupnewrun')

    if not success:
        sys.exit(-1)

    print(Fore.GREEN + 'Setup tests for convection in a fixed uniform porous medium' + Style.RESET_ALL)
    physical_problem = 'convectionDB'
    # Nz_uniform = 128
    # if Nz_vm <= 0:
    #    Nz_vm = Nz_uniform #int(float(Nz_uniform) / 2)

    amr_setup = []
    if nz_uniform > 0:
        amr_setup.append({'max_level': 0, 'ref_rat': 2, 'run_types': ['uniform'], 'Nzs': [int(nz_uniform/2), nz_uniform, 2*nz_uniform]})

    if nz_vm > 0:

        these_nzs = [nz_vm, 2*nz_vm, 4*nz_vm]
        amr_setup.append({'max_level': 1, 'ref_rat': 2, 'run_types': ['variable'], 'Nzs': these_nzs})

    num_procs = [1, 1, 1]
    # chi = 0.4

    # Try and speed things up for now, should eventually make this criteria smaller
    # extra_params = {'main.steady_state': 1e-4}
    extra_params = {'main.BCAccuracy': bc_accuracy,
                    'main.uDeluMethod': u_del_u_method,
                    'main.advectionMethod': advection_method,
                    'main.cfl': cfl,
                    'main.initial_cfl': cfl / 10}

    base_dataFolder = os.path.join(base_output_dir, 'ConvectionDB-cfl' + str(cfl))

    da_ra_vals = [{'Da': 1e-6, 'RaT': [1e7, 1e8, 1e9],
                   'lebars': [1.08, 3.07, 12.9]},
                  {'Da': 1e-2, 'RaT': [1e3, 1e4, 1e5, 5e5],
                   'lebars': [1.01, 1.41, 3.17, 5.24]}]

    if fixed_da > 0.0:
        if fixed_ra > 0.0:
            da_ra_vals = [{'Da': fixed_da, 'RaT': [fixed_ra], 'lebars': [float('NaN')]}]
        elif fixed_da == 1e-6:
            da_ra_vals = [{'Da': 1e-6, 'RaT': [1e7, 1e8, 1e9],
                           'lebars': [1.08, 3.07, 12.9]} ]
        elif fixed_da == 1e-2:
            da_ra_vals = [{'Da': 1e-2, 'RaT': [1e3, 1e4, 1e5, 5e5],
                           'lebars': [1.01, 1.41, 3.17, 5.24]}]

    # [1.01, 1.41, 3.17, 5.24];
    # [1.08, 3.07, 12.9];
    all_job_ids = []

    analysis_command = matlab_command + ' " '

    ra_format = "%1.1e"
    da_format = "%1.1e"
    chi_format = "%1.1f"

    chi_str = chi_format % chi

    for da_ra in da_ra_vals:
        da = da_ra['Da']
        da_str = da_format % da
        nu_lebars = [str(a) for a in da_ra['lebars']]

        # extra_params['main.max_time'] = 0.1 / float(Da)
        #        extra_params['main.max_dt'] = ?

        for ra in da_ra['RaT']:
            extra_params['parameters.rayleighTemp'] = ra
            extra_params['main.plot_interval'] = -1
            extra_params['main.plot_period'] = 0.05
            extra_params['main.checkpoint_interval'] = 2000
            
            # this is for testing purposes
            # extra_params['main.max_step'] = 10
           
            if chi < 1.0:
                extra_params['parameters.darcy'] = da * pow(1 - chi, 2) / pow(chi, 3.0)
            else:
                extra_params['parameters.darcy'] = da

            extra_params['bc.porosityHiVal'] = str(chi) + ' ' + str(chi)  # 0.4 0.4
            extra_params['bc.porosityLoVal'] = str(chi) + ' ' + str(chi)

            # output_dir = 'chi' + str(chi) + '-Da' + str(Da) + '-Ra' + str(extra_params['parameters.rayleighTemp'])
            ra_str = ra_format % ra
            output_dir = "chi" + chi_str + "-Da" + da_str + "-Ra" + ra_str

            this_data_folder = os.path.join(base_dataFolder, output_dir)
            job_ids = runTest(this_data_folder, physical_problem, uniform_porous_resolution_specific_params,
                              amr_setup, num_procs, '', extra_params, restart_from_low_res=True)
            
            all_job_ids = all_job_ids + job_ids

        ra_str_vals = [ra_format % a for a in da_ra['RaT']]
        ra_str = '{\'' + '\',\''.join(ra_str_vals) + '\'}'

        analysis_command = analysis_command + ' compileNu(\'' + base_dataFolder + '\', \'' + chi_str + '\', \'' + da_str + '\', ' + ra_str + ', ' + str(
            nz_uniform) + ', ' + str(nz_vm) + ', [' + ','.join(nu_lebars) + ']);'

        # Now do analysis
    analysis_command = analysis_command + ' exit; " '

    run_analysis_name = 'runAnalysis.sh'

    job_name = physical_problem + '-analysis'
    s = BatchJob(base_dataFolder, job_name, '')

    s.set_dependency(all_job_ids)
    s.set_custom_command(analysis_command)
    # s.write_slurm_file(run_analysis_name)
    s.run_task(run_analysis_name)
    print(Fore.GREEN + 'Submitted analysis job \n' + Fore.RESET)


if __name__ == "__main__":
    test_uniform_porous_convection(sys.argv[1:])
