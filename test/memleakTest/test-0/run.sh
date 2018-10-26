#!/bin/bash 
# Set your minimum acceptable walltime, format: day-hours:minutes:seconds 
#SBATCH --time=1-0:0:0
#SBATCH --partition=shared,priority-ocean
# Set name of job shown in squeue
#SBATCH --job-name memleakTest
# Request CPU resources
#SBATCH --ntasks=1                  # Number of MPI ranks
#SBATCH --cpus-per-task=1            # Number of cores per MPI rank 
#SBATCH --nodes=1                    # Number of nodes
#SBATCH --ntasks-per-node=1         # How many tasks on each node
#SBATCH --ntasks-per-socket=0.0        # How many tasks on each CPU or socket (not sure what this really means)
#SBATCH --distribution=cyclic:cyclic # Distribute tasks cyclically on nodes and sockets
# Memory usage (MB)
#SBATCH --mem-per-cpu=5000
#SBATCH --output=/home/parkinsonjl/convection-in-sea-ice/MushyLayer/memleakTest/test-0/sbatch.out   # Standard output and error log
cd /home/parkinsonjl/convection-in-sea-ice/MushyLayer/memleakTest/test-0;  valgrind --leak-check=yes /home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/mushyLayer2d.Linux.64.mpiCC.gfortran.OPT.MPI.ex /home/parkinsonjl/convection-in-sea-ice/MushyLayer/memleakTest/test-0/inputs