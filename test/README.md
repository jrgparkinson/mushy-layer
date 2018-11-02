# Testing

# Introduction
This directory contains python scripts for running various test problems which demonstrate the accuracy of the code. 

# Setup
Before running the test problems, you must first specify a directory to save the output to. This is done by the `getBaseOutputDir()` function in `mushyLayerRunUtils.py`. 
You should also add this directory (`/path/to/mushy-layer/test/`) to the PYTHONPATH.

Analysis of the output is done via matlab scripts located in `/mushy-layer/matlab/`. You should add this directory to your matlab path via your `startup.m` file in order to run the analysis.

Finally, the python scripts assume that the slurm job queuing system is accessible. If it is not, you will have to run the batch files which it creates manually.

# Running the problems
There are five test problems, which can each be run individually by calling `python filename.py (options)`:

* testDiffusiveSolidification.py
* testUniformPorousConvection.py
* testFixedPorousHole.py
* testPorousMushyHole.py
* testFixedChillHeleShaw.py

See the individual files for the options available.

Alternatively you may run all five test problems via `python runMethodPaperTests.py`.

The python scripts will create the relevant folders to hold the output, setup the parameter input files, submit jobs to run the simulations to the slurm queueing system, and also submit a job to run the relevant analysis code once all the jobs have finished. The output from analysis can be found in the sbatch.out file, so a typical directory structure might end up looking like

/baseOutputDir (specified by the user)
  /NoFlow
    runAnalysis.sh
    sbatch.out - output from the analysis
    Various folders containing different simulations, each containing:
      run.sh
      pout.0
      time.table.0
      pltfiles.hdf5
      chkfile.hdf5
  /FixedPorousHole
    ... same as above