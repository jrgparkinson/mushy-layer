# This test is to check that simulations which have been stopped then restarted 
# produce the same results as those which have been run continuously
# ABANDONDED because chombo compare can't do AMR to AMR comparison

base_dir = '/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle'
continuous_output = base_dir + '/output'
restarted_output = base_dir + '/output-restart'
f

chombo_compare_exec = '/home/parkinsonjl/Chombo-trunk/lib/util/ChomboCompare/compare2d.Linux.64.mpiCC.gfortran.DEBUG.MPI.ex'

restart_timestep = 50
final_timestep = 90

# Run chombo compare on all these files

for plot_i in range(restart_timestep, final_timestep):

    command = (chombo_compare_exec + ' restartInputs.compare' + 
     'compare.exactRoot = ' + continuous_output)
    print(command)
    print('\n')
