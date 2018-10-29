# AMR Mushy Layer

# Introduction
This is code to simulate flow in a reactive porous media, a "mushy layer", within a square box.

For more information on the equations, numerical method, and motivation, see the [Paper](#) (link to come) 

# Prerequisites
## Required
* [Chombo](https://commons.lbl.gov/display/chombo/Chombo+-+Software+for+Adaptive+Solutions+of+Partial+Differential+Equations) installation ([how to download Chombo](https://anag-repo.lbl.gov/chombo-3.2/access.html))
* git (for downloading and updating code)

## Optional
* python (for setting up the test problems and running suites of simulations)
* MATLAB (for analysis of the test problems)
* doxygen (for generating documentation)

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

# Testing
`/test/` contains python scripts for running various test problems, which demonstrate the accuracy and efficiency of the code. You can run all the problems via the script `runMethodsPaperTests.py`. 

The analysis is performed by matlab scripts located in `/matlab/*`, which must be added to your matlab path to be executed:

```matlab
addpath(genpath('/path/to/matlab/'))
```

# Running
The executable in `/execSubcycle/` requires inputs in order to run. These can either be supplied at the command line or, more sensibly, by specifying the location of a file which contains a list of inputs, e.g.

```console
$ cd execSubcycle
$ ./mushyLayer2d.Linux.64.mpiCC.gfortran.OPT.MPI.ex inputs
```

where inputs a file in /execSubcycle/

#Documentation
Documentation can be generated using Doyxgen if it is installed. A default configuration file can be found in the `/docs/` directory, which can be run via:

```console
$ cd /docs/
$ doxygen
```


