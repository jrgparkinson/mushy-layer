# Run all convergence tests for the methods paper
import os, sys
from colorama import Fore, Style
import getopt

from runAMRConvergenceTest import runTest
from SlurmTask import SlurmTask
from mushyLayerRunUtils import getBaseOutputDir, getMatlabBaseCommand

######################################
# 2) Convection in a fixed porous medium
# Just run for a single grid resolution and compare to previously published values
######################################

def testUniformPorousConvection(argv):

    # Default vals:
    cfl = 0.4
    Nz_uniform = 128
    Nz_vm = -1

    try:
       opts, args = getopt.getopt(argv,"n:v:c:")
    except getopt.GetoptError as err:
        print(str(err))
        print('testUniformPorousConvection.py -n<num uniform mesh grid points> -v<num variable mesh points> -c<cfl>')
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-n"):
            Nz = int(arg)
        elif opt in ("-c"):
            cfl = float(arg)
        elif opt in ("-v"):
            Nz_vm = int(arg)

    base_output_dir = getBaseOutputDir()
    matlab_command = getMatlabBaseCommand()

    print(Fore.GREEN + 'Setup tests for convection in a fixed uniform porous medium' + Style.RESET_ALL)
    physicalProblem = 'convectionDB'
    #Nz_uniform = 128
    if Nz_vm <= 0:
        Nz_vm = Nz_uniform #int(float(Nz_uniform) / 2)

    AMRSetup = [{'max_level': 0, 'ref_rat': 2, 'run_types': ['uniform'], 'Nzs': [Nz_uniform]},
                {'max_level': 1, 'ref_rat': 2, 'run_types': ['variable'], 'Nzs': [Nz_vm]}]


    num_procs = [4]
    chi = 0.4

    # Try and speed things up for now, should eventually make this criteria smaller
    extra_params = {'main.steady_state': 1e-4}
    #cfl = 0.45
    extra_params['main.cfl'] = cfl
    extra_params['main.initial_cfl'] = cfl / 10
    base_dataFolder = os.path.join(base_output_dir, 'ConvectionDB-cfl' + str(cfl))

    Da_Ra_vals = [{'Da': 1e-6, 'RaT': [1e7, 1e8, 1e9], 'lebars': [1.08, 3.07, 12.9]},
                  {'Da': 1e-2, 'RaT': [1e3, 1e4, 1e5, 5e5], 'lebars': [1.01, 1.41, 3.17, 5.24]}]

    # [1.01, 1.41, 3.17, 5.24];
    # [1.08, 3.07, 12.9];
    all_job_ids = []

    analysis_command = matlab_command + ' " '

    for Da_Ra in Da_Ra_vals:
        Da = Da_Ra['Da']
        NuLebars = [str(a) for a in Da_Ra['lebars']]

        extra_params['main.max_time'] = 0.1 / float(Da)
#        extra_params['main.max_dt'] = ?

        for Ra in Da_Ra['RaT']:
            extra_params['parameters.rayleighTemp'] = Ra

            extra_params['parameters.darcy'] = Da * pow(1 - chi, 2) / pow(chi, 3.0)

            extra_params['bc.porosityHiVal'] = str(chi) + ' ' + str(chi)  # 0.4 0.4
            extra_params['bc.porosityLoVal'] = str(chi) + ' ' + str(chi)

            #output_dir = 'chi' + str(chi) + '-Da' + str(Da) + '-Ra' + str(extra_params['parameters.rayleighTemp'])
            output_dir = "chi%1.1f-Da%1.1e-Ra%1.1e" % (chi, Da, extra_params['parameters.rayleighTemp'])

            # extra_params = {}
            thisDataFolder = os.path.join(base_dataFolder, output_dir)
            job_ids = runTest(thisDataFolder, physicalProblem, AMRSetup, num_procs, '', extra_params)
            all_job_ids = all_job_ids + job_ids

        Ra_str_vals = [str(a) for a in Da_Ra['RaT']]
        Ra_str = '{\'' + '\',\''.join(Ra_str_vals) + '\'}'

        analysis_command = analysis_command + ' compileNu(\'' + base_dataFolder + '\', \'' + str(chi) + '\', \'' + str(
            Da) + '\', ' + Ra_str + ', ' + str(Nz_uniform) + ', ' + str(Nz_vm) + ', [' + ','.join(NuLebars) + ']);'

        # Now do analysis
    analysis_command = analysis_command + ' exit; " '

    runAnalysisName = 'runAnalysis.sh'


    jobName = physicalProblem + '-analysis'
    s = SlurmTask(base_dataFolder, jobName, '')

    s.setDependency(all_job_ids)
    s.setCustomCommand(analysis_command)
    s.writeSlurmFile(runAnalysisName)
    s.runTask(runAnalysisName)
    print(Fore.GREEN + 'Submitted analysis job \n' + Fore.RESET)



if __name__ == "__main__":
    testUniformPorousConvection(sys.argv[1:])