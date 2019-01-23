# AMR Mushy Layer

# Introduction
This is code to simulate flow in a reactive porous media, a "mushy layer", within a square box.

For more information on the equations, numerical method, and motivation, see the [Paper](#) (link to come) 

# Prerequisites
## Required
* [Chombo](https://commons.lbl.gov/display/chombo/Chombo+-+Software+for+Adaptive+Solutions+of+Partial+Differential+Equations) installation ([how to download Chombo](https://anag-repo.lbl.gov/chombo-3.2/access.html)). You will need to configure Chombo with HDF5 support in order to write out data files.

## Optional
* git. For downloading and updating code.
* python. For setting up the test problems.
* MATLAB. For analysis of the test problems.
* doxygen. For generating documentation.
* slurm. The code for running the test problems will create a set of batch scripts and try to submit these to the slurm queuing system. If you do not have slurm installed, you'll need to run the batch scripts manually.
* [visit](http://www.nersc.gov/users/data-analytics/data-visualization/visit-2/). Visualisation software for Chombo AMR files.

# Accessing the code
Cloning the repository from GitHub is straightforward:

```console
$ git clone https://github.com/jrgparkinson/mushy-layer.git
```

To fetch and merge any changes to the library, use `git pull`.


# Installation
To compile the code, you must first have a working Chombo installation. Define the path to your chombo library folder in the GNUMakefile (`/execSubcyle/GNUMakefile/`) then you should be able to compile by

```console
$ cd /execSubcycle/
$ make all
```

You should also add an environment variable `MUSHY_LAYER_DIR` which points to `/path/to/mushy-layer/` for some of the python scripts to work. E.g.

```console
EXPORT MUSHY_LAYER_DIR=/path/to/mushy-layer/
```

# Testing
`/test/` contains python scripts for running various test problems, which demonstrate the accuracy and efficiency of the code. You can run all the problems via the script `runMethodsPaperTests.py`. 

You must ensure that the `/test/` directory is in your `PYTHONPATH`, i.e. place the following in your ~/.bashrc or similar

```console
EXPORT PYTHONPATH=PYTHONPATH:/path/to/mushy-layer/test
```

The analysis is performed by matlab scripts located in `/matlab/*`, which must be added to your matlab path to be executed:

```matlab
addpath(genpath('/path/to/matlab/'))
```

Note that running all the tests will take some time unless you have lots of processors available (i.e. multiple days) and will take up a reasonable amount of disk space (~50GB).

# Running
The executable in `/execSubcycle/` requires inputs in order to run. These can either be supplied at the command line or, more sensibly, by specifying the location of a file which contains a list of inputs, e.g.

```console
$ cd execSubcycle
$ ./mushyLayer2d.Linux.64.mpiCC.gfortran.OPT.MPI.ex inputs
```

where inputs a file in /execSubcycle/

# Documentation
Documentation can be generated using Doyxgen if it is installed. A default configuration file can be found in the `/docs/` directory, which can be run via:

```console
$ cd /docs/
$ doxygen
```

# Extra utilities
This repository contains various other pieces of code for running simulations. 

`/setupNewRun/` takes a checkpoint file and creates a new file with the same data on a domain with a different width, which is useful for computing optimal states.
`/postProcess/` loads checkpoint files and computes various diagnostics that were not run during the simulations.
`/VisitPatch/` describes a small patch for the VISIT software to allow you to open checkpoint files as well as plot files.
`/params/` contains input files for different types of simulations. They may not all work, sorry.
`/grids/` contains gridfiles which can be loaded via `main.gridfile=/path/to/gridfile` in an inputs file, and sets a fixed variable mesh (i.e. not adaptive) for simulations.
`/mk/` contains some custom Makefile options for compiling on some of the machines in AOPP at the University of Oxford.




