from MushyLayerRun import MushyLayerRun
import os

os.environ["CH_TIMER"] = "1"
dropbox_data_folder = "/home/parkinsonjl/Dropbox/DPhil/Data/"

# flag = 'noStefan-godunov-'
flag = 'fixed-'
output_dir = 'bm1-' + flag
plot_prefix = 'bm1-' + flag
check_prefix = 'bm1-' + flag + 'chk-'
num_proc = 1
max_level = 0
plot_interval = 1
check_interval = 50
parameters = {}
res = [128]

parameters['main.frameAdvectionMethod'] = '1' # 0-godunov, 1-upwind

parameters['parameters.stefan'] = 5.7 #5.7
parameters['parameters.compositionRatio'] = '1.4'
parameters['parameters.K'] = '1'
parameters['parameters.specificHeatRatio'] = '1'
parameters['parameters.lewis'] = 10e200  # (infinite)
parameters['parameters.V'] = 1e-6  # 1e-6
parameters['parameters.deltaTemp'] = 25
parameters['parameters.bottomTemp']= -20
parameters['parameters.topTemp']  = 15  #15
parameters['parameters.eutecticTemp']        = -20
parameters['parameters.height']             = 0.5
parameters['parameters.initialComposition']    = 0.7
parameters['parameters.eutecticComposition']    = 0.2 
parameters['parameters.liquidusSlope']    = 50
parameters['parameters.heatConductivityLiquid']   = 0.5         # W/(m K)
parameters['parameters.specificHeatLiquid']         = 3.5e3     # J/(kg K)
parameters['parameters.liquidDensity']            = 1000        # kg/(m^3)
parameters['parameters.heatConductivityRatio'] = 1
parameters['parameters.heleShawCoolingCoeff'] = 0 #we don't consider a hele-shaw setup yet

parameters['main.problem_type'] = '1'
parameters['main.enforceAnalyticSoln'] = 'false'

parameters['parameters.bottomTemp']   =    -20.0    
parameters['parameters.eutecticTemp']     =   -20.0
parameters['parameters.topTemp']     =   15.0        

parameters['parameters.initialComposition']       = 0.25  

parameters['main.tagging_val'] = 0.2
parameters['main.tagging_scalar_var'] = 5

gridType = MushyLayerRun.NARROW

ml_run = MushyLayerRun(output_dir,  plot_prefix,  check_prefix,  max_level,  num_proc,  plot_interval, check_interval, gridType, parameters)


#ml_run.continueRun('64-0', 'bm1-chk-64pts0200.2d.hdf5')
#ml_run.calcAnalyticSoln([1024])

#parameters['main.printAnalyticSoln'] = 'true'

ml_run.run(res)

#ml_run.calc_analytic_soln(res)

 # Copy data to dropbox so we can analyse it
if ml_run.machine_specific.has_dropbox():
    ml_run.copy_final_output(res, dropbox_data_folder)

# Copy data from gyre servers to cloud
#if ml_run.machine_specific.is_gyre():
ml_run.copy_all_output_to_cloud(res)


#ml_run.processOutput('bm1/32-15/', 'bm1-32pts')

# ml_run.runForFixedGrid([16,32])




