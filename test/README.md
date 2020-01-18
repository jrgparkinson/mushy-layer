# Testing

## Introduction
This directory contains python scripts for running various test problems which demonstrate the accuracy of the code.

There are also regression tests in `regression` and unit tests for the python code in `unit`. Both of these are run automatically following each push to the repository, with unit test coveraged published `https://jrgparkinson.github.io/mushy-layer/test/unit/coverage/BRANCH-NAME` e.g. [https://jrgparkinson.github.io/mushy-layer/test/unit/coverage/development].

## Setup
These tests require the code in `mushy-layer/execSubcycle/` and `mushy-layer/setupNewRun/` to be built.

Before running the test problems, you must first specify a directory to save the output to. This is done by the `get_base_output_dir()` function in `mushyLayerRunUtils.py`. 
You should also add the directory (`/path/to/mushy-layer/test/`) to your PYTHONPATH.

Some of the error calculation is done using the Chombo compare utility located in `/path/to/chombo/lib/util/chomboCompare`. Make sure you have compiled this program, and make sure you have defined the environmental variable `CHOMBO_HOME=/path/to/chombo/`.

Analysis of the output is done via matlab scripts located in `/mushy-layer/matlab/`. You should add this directory to your matlab path via your `startup.m` file in order to run the analysis, e.g.
```matlab
addpath(genpath(['~/mushy-layer/matlab/']))
```

The python scripts assume that the slurm job queuing system is accessible. If it is not, you will have to run the batch files which it creates manually. Alternatively, open up the `test/BatchJob.py` file and edit the `BatchJob.write_slurm_file(...)` and `BatchJob.run_task(...)` methods so that they work with your setup.

## Running the problems
There are five test problems, which can each be run individually by 
```console
$ python filename.py [options]
```

*   testDiffusiveSolidification.py
*   testUniformPorousConvection.py
*   testFixedPorousHole.py
*   testPorousMushyHole.py
*   testFixedChillHeleShaw.py

See the individual files for the options available.

Alternatively you may run all five test problems via 

```console
$ python runMethodPaperTests.py
```

There are a lot of individual simulations that need to be run. It is highly recommended that you do this on a machine with lots of cores available so that you can the simulations on different cores in parallel, in which case most jobs should finish in a few hours. The Nusselt number tests with large Rayleigh number on high resolution grids require longer (> 1 day).

The python scripts will create the relevant folders to hold the output, setup the parameter input files, submit jobs to run the simulations to the slurm queueing system, and also submit a job to run the relevant analysis code once all the jobs have finished. The output from analysis can be found in the sbatch.out file, so a typical directory structure might end up looking like

```text
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
```
 
## Figures and analysis
Having run the test problems, most of the figures should already be created, usually in the same directory as the data for each test problem resides in. You can (re)generate the figures used in the Methods paper using
```console
$ python makeFigures.py
```
Use the matlab script `compileNuV2(base_dir)` to compile all the nusselt numbers, e.g.
```matlab
compileNuV2('/path/to/TestDir/ConvectionDB-CFL0.1')
```
