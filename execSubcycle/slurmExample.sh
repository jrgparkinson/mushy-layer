#!/bin/bash 
# This is an example script for running some mushy layer code using slurm
# First setup the slurm options, then list the commands to execute
# This batch script will only submit to the gyres, but you can also submit to the atmnodes (shared) and priority ocean queue with --partition=legacy,shared,priority-ocean
# Set your minimum acceptable walltime, format: day-hours:minutes:seconds 
#SBATCH --time=7-0:0:0
#SBATCH --partition=legacy
#SBATCH -x cirrus[1-4],yau[1-6],yau[8-10],yau12
# Set name of job shown in squeue
#SBATCH --job-name pts256CR76.000RaC80
# Request CPU resources
#SBATCH --ntasks=2                  # Number of MPI ranks
#SBATCH --cpus-per-task=1            # Number of cores per MPI rank 
#SBATCH --nodes=1                    # Number of nodes
#SBATCH --ntasks-per-node=2         # How many tasks on each node
#SBATCH --ntasks-per-socket=1        # How many tasks on each CPU or socket (not sure what this really means)
#SBATCH --distribution=cyclic:cyclic # Distribute tasks cyclically on nodes and sockets
# Memory usage (MB)
#SBATCH --mem-per-cpu=4000
#SBATCH --output=/network/group/aopp/oceans/AW002_PARKINSON_MUSH/optimalStates-restructured/Le200/CR76.000/RaC80/pts256-0/sbatch.out   # Standard output and error log

# Could do some pre-processing here if required e.g.
python /run/preprocess/code.py

# Go to the directory where the inputs file is located, and where we want to save the output
cd /network/group/aopp/oceans/AW002_PARKINSON_MUSH/optimalStates-restructured/Le200/CR76.000/RaC80/pts256-0; 
 
# Run the correct version of the code (depends if we're on gyre or not)
hs=$HOSTNAME 
if [[ $hs == *"gyre"* ]]; then 
	# I always forget to build the code on gyre, so adding a line here to make sure that's happened
 /home/parkinsonj/mushy-layer/execSubcycle/buildMushyLayer.sh; 
 
 # Run the code in paralle with 2 processors (gyre version)
 mpirun -np 2 /home/parkinsonj/mushy-layer/execSubcycle/mushyLayer2d.Linux.64.mpiCC.gfortran.OPT.MPI.GYRE.ex /network/group/aopp/oceans/AW002_PARKINSON_MUSH/optimalStates-restructured/Le200/CR76.000/RaC80/pts256-0/inputs 
else

# Run the code in paralle with 2 processors (non-gyre version)
 mpirun -np 2 /home/parkinsonj/mushy-layer/execSubcycle/mushyLayer2d.Linux.64.mpiCC.gfortran.OPT.MPI.ex /network/group/aopp/oceans/AW002_PARKINSON_MUSH/optimalStates-restructured/Le200/CR76.000/RaC80/pts256-0/inputs
fi

# Do some post processing
cd /home/parkinsonj/phd/; source env/bin/activate; python /home/parkinsonj/phd/python/analysis_postprocess/postProcessSingleFolder.py -R 80.000000 -C 76.000000 -P 256 -e 0 -r -a
