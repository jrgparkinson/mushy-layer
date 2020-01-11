# SOFTBALL: SOLidification, Flow and Thermodynamics in Binary ALLoys

This is code to simulate flow in a reactive porous media, a "mushy layer". This if a finite volume/difference code, built on the [Chombo library](https://commons.lbl.gov/display/chombo/Chombo+-+Software+for+Adaptive+Solutions+of+Partial+Differential+Equations), which solves the governing conservation equations on a cuboidal domains with Adaptive Mesh Refinement (AMR).
 
In this readme we first present some instructions on downloading and installing the code, then describe how to recreate the test problems and figures found in [Parkinson et. al. (2020)](https://www.sciencedirect.com/science/article/pii/S2590055219300599). The paper is a good reference for more information on the physics of this model.

If you have any issues, please contact [james.parkinson@physics.ox.ac.uk][james.parkinson@physics.ox.ac.uk].

## Prerequisites
### Required
*   [Chombo](https://commons.lbl.gov/display/chombo/Chombo+-+Software+for+Adaptive+Solutions+of+Partial+Differential+Equations) installation ([how to download Chombo](https://anag-repo.lbl.gov/chombo-3.2/access.html)). You will need to configure Chombo with HDF5 support in order to write out data files. This code requires Chombo version 3.2, patch 6 or above. Earlier versions of Chombo do not include the `AMRFASMultigrid` functionality, required for solving nonlinear elliptic equations.

### Optional
*   git. For downloading and updating code.
*   python. For setting up the test problems.
*   MATLAB. For analysis of the test problems.
*   doxygen. For generating documentation.
*   slurm. The code for running the test problems will create a set of batch scripts and try to submit these to the slurm queuing system. If you do not have slurm installed, you'll need to run the batch scripts manually.
*   [visit](http://www.nersc.gov/users/data-analytics/data-visualization/visit-2/). Visualisation software for Chombo AMR files.

## Accessing the code
Cloning the repository from GitHub is straightforward:

```console
$ git clone https://github.com/jrgparkinson/mushy-layer.git mushy-layer
```

To fetch and merge any changes to the library, use `git pull`.

## Installation
To compile the code, you must first have a working Chombo installation (see the ChomboInstallationGuide.md file in the docs folder). Define the path to your chombo library folder in the GNUMakefile (`/execSubcyle/GNUMakefile`) then you should be able to compile by

```console
$ cd ~/mushy-layer/execSubcycle/
$ make all
```

You should also add an environment variable `MUSHY_LAYER_DIR` which points to `/path/to/mushy-layer/` for some of the python scripts to work. E.g.

```console
export MUSHY_LAYER_DIR=/path/to/mushy-layer/
```

## Running
The executable in `/execSubcycle/` requires inputs in order to run. These can either be supplied at the command line or, more sensibly, by specifying the location of a file which contains a list of inputs, e.g.

```console
$ cd execSubcycle
$ ./mushyLayer2d.Linux.64.mpiCC.gfortran.OPT.MPI.ex inputs
```

where inputs is a file in `/execSubcycle/`.

For more details on how to setup your inputs file for the problem you're interested in, see `/docs/RunningTheMushyLayerCode.md/`.

There are example problems, with input files, in the `/examples/` directory.

For running in parallel (with 2 processors in this example), run the code like

```console
$ cd execSubcycle
$ mpirun -np 2 ./mushyLayer2d.Linux.64.mpiCC.gfortran.OPT.MPI.ex inputs
```

## Documentation
Documentation can be generated using [Doyxgen](http://www.doxygen.nl/) if it is installed. A default configuration file can be found in the `/docs/` directory, which can be run via:

```console
$ cd /docs/
$ doxygen
```

A reasonably up to date version of the documentation can be found [online][https://amr-softball.github.io/doc/html/index.html], if you do not want to compile it yourself. The Classes section is particularly useful.

## Extra utilities
This repository contains various other pieces of code for running simulations. 

`/setupNewRun/` takes a checkpoint file and creates a new file with the same data on a domain with a different width, which is useful for computing optimal states.
`/postProcess/` loads checkpoint files and computes various diagnostics that were not run during the simulations.
`/VisitPatch/` describes a small patch for the VISIT software to allow you to open checkpoint files as well as plot files.
`/params/` contains input files for different types of simulations. They may not all work, sorry.
`/grids/` contains gridfiles which can be loaded via `main.gridfile=/path/to/gridfile` in an inputs file, and sets a fixed variable mesh (i.e. not adaptive) for simulations.
`/mk/` contains some custom Makefile options for compiling on some of the machines in AOPP at the University of Oxford.

The code in `/setupNewRun/` and `/postProcess/` needs to be compiled like the code in `/execSubcycle/`. First update the `GNUMakefile` files in each subdirectory, then run `make all` as before.

## Source code
The source code is spread across a number of directories, which are briefly summarised here.

`/BCutil/` contains code for implementing boundary conditions.
`/util/` contains code for various non-mushy layer specific operations, e.g. computing gradients, divergences etc.
`/src/` contains mushy layer code
`/srcSubcycle/` contains mushy layer code for the subcycled algorithm
`/execSubcyle/` contains the driver code the subcycled application
`/srcNonSubcycle/` contains mushy layer code for the non-subcycled algorithm (currently broken)
`/execNonSubcycle/` contains the driver code for the non-subcycle application (currently broken)

### Methods paper
We have written a paper describing the physics of this model. A draft can be found in the `/docs/` directory. In this paper, we demonstrate the accuracy and convergence properties of our code. Everything you need to run these tests yourself should be found in the `/tests/` directory, as described below. 

## Testing
`/test/` contains python scripts for running various test problems, which demonstrate the accuracy and efficiency of the code. You can run all the problems via the script `runMethodsPaperTests.py`. For more information, see the README located in the `/test/` directory.

You must ensure that the `/test/` directory is in your `PYTHONPATH`, i.e. place the following in your ~/.bashrc or similar

```console
EXPORT PYTHONPATH=PYTHONPATH:/path/to/mushy-layer/test
```

The analysis is performed by matlab scripts located in `/matlab/*`, which must be added to your matlab path to be executed:

```matlab
addpath(genpath('/path/to/matlab/'))
```

Note that running all the tests will take some time unless you have lots of processors available (i.e. multiple days) and will take up a reasonable amount of disk space (~50GB).

## Figures
Having run the test problems, some of the figures found in the paper will have been created automatically in the subdirectory for each tests problem (e.g. `/path/to/test_output/FixedPorousHole-1proc/Fig6Error-xDarcy_velocity-L2.eps`). The python script located at `mushy-layer/test/makeFigures.py` should then make the remaining figures for you, mainly by calling various matlab scripts. 

## The mushy layer code
In this section we briefly describe how the code works with the Chombo library so solve equations on a hierarchy of adaptive meshes. If you are planning on making changes to the code, or want to see which bits solve which equations, this might be a useful starting point. If you are going to look at the code regularly, you are highly recommended to use an Integrated Development Environment (IDE) - e.g. [Eclipse](https://www.eclipse.org/downloads/packages/release/2018-12/r/eclipse-ide-cc-developers). Brief setup instructions for Eclipse can be found below.

The mushy layer code starts in `execSubcycle/mushyLayer.cpp`. The function `mushyLayer(...)` calls various other methods from `src/MushyLayerUtils.cpp` and `srcSubcycle/MushyLayerSubcycleUtils.cpp` to 
a) Create the problem domain
b) Create a 'Mushy Layer Factory' object which tells Chombo how to make a new Mushy Layer object for handling integration of the governing equations on one level of refinement
c) Create an `AMR` object which handles solving the equations over a hierarchy of levels.

The code then tells the AMR object to run until a certain time or number of steps. The AMR object is the heart of Chombo. It continually evolves the solution, periodically calling into functions from the Mushy Layer code to do things, e.g.

*   Initialisation - most of the Mushy layer specific stuff is contained in `srcSubcycle/AMRLevelMushyLayerInit.cpp`.
*   Advancing the solution - this is done `srcSubcycle/AMRLevelMushyLayer.cpp` in the `advance(...)` function (which subsequently calls lots of other mushy layer functions).
*   Finding the current maximum allowed timestep - implemented in `srcSubcycle/AMRLevelMushyLayer.cpp :: computeDt()`.
*   Regridding, if necessary - see `srcSubcycle/AMRLevelMushyLayerRegrid.cpp :: regrid(..)`.
*   Post timestep synchronisation - mainly contained in `srcSubcycle/AMRLevelMushyLayerSync.cpp :: postTimestep(..)`

Within the `advance(..)` method, the key pieces are:

*   `computeAdvectionVelocities(..)` - computes face centered velocities for doing advection with.
*   `exitStatus = multiCompAdvectDiffuse(HC_old, HC_new, srcMultiComp, doFRUpdates, doAdvectiveSrc);` - updates enthalpy and salinity fields.
*   `computeCCvelocity(advectionSourceTerm, m_time-m_dt, m_dt, doFRupdates, doProjection, compute_uDelU);` - compute new cell centred velocity.

## Notes on adding new code
For those new to this code/Chombo, here are some tips:

*   An AMRLevelMushyLayer object exists on a particular level of refinement. You can access the coarser and finer levels via `getCoarserLevel()`, `getFinerLevel()`.

*   During a timestep, `m_time` is the new time ($t^{n+1} = t^n + \Delta t$).

*   Scalar and vector fields are stored in arrays `m_scalarNew`, `m_vectorNew`, indexed by the field name, i.e. enthalpy on a level of refinement is `m_scalarNew[ScalarVars::m_enthalpy]`

*   When using scalar and vector fields, you can ensure all the boundary conditions and ghost cells are filled by calling `fillScalars(..)`, e.g.
` fillScalars(*m_scalarNew[ScalarVars::m_enthalpy], m_time, ScalarVars::m_enthalpy);`

## Eclipse + Chombo + Mushy Layer
[Eclipse](https://www.eclipse.org/downloads/packages/release/2018-12/r/eclipse-ide-cc-developers) is an Integrated Development Environment (IDE). In short, it is a powerful code editor (syntax highlighting, advanced find/replace) which also understands the structure of a program. For codes like this one which are spread over multiple files in multiple folders, and use external libraries, it is invaluable. Useful features include:

*   Find where a variable or function is defined (right click on variable/function -> Open Declaration)
*   Finding all usages of a variable or function (right click -> References -> workspace)
*   Debugging. Insert breakpoints via the top menu: Run  -> Toggle Breakpoint.

In short, it's a bit like MATLAB (but for C++).

Setting up Eclipse with the Mushy Layer is not difficult, but also requires quite a few steps (should take 30 mins, but took me many hours to work them all out from scratch). The following steps should do the job. We assume that you are already able to compile and run the Mushy Layer code.

1.  Create a new project for the Chombo code. Go File->New->Makefile project with existing code. Click Browse and navigate to your /chombo/lib directory	then hit ok. Give the project a suitable name and choose the Linux GCC toolchain.
 
2.  We now need to tell eclipse how to build the Chombo code. Eclipse may not inherit your environmental variables, so you'll need to define these explicitly. To do this, right click on your project in the project window (usually on the left of the screen), choose Properties, and then navigate to C/C++ Build->Environment. Add 
`LD_LIBRARY_PATH = /usr/local/hdf5-1.8.21p-v18/lib:`
The first should point to your hdf5 installation.

3.  Create a new project for your mushy layer code . As before, go File->New->Makefile project with existing code. Set the code location to `/path/to/mushy-layer`, and choose a sensible name (we went for 'MushyLayer'). Again choose the Linux GCC toolchain.
 
4.  Set Chombo as a dependency of mushy layer. Open the project properties for the Mushy Layer project (right click in the project window->properties). Choose 'Project references' on the left and select your Chombo project.

5.  Setup build profile for mushy layer. Go to C/C++ build->Environment Variables and check that the same environment variables you set for the Chombo project are also here. Now go to C/C++ build and, under the Builder settings tab, set the build directory to `${workspace_loc:/MushyLayer}/execSubcycle/`. When you use this build profile, you will create an executable using the default options in your Make.defs.local file.

6.  Create a debugging build profile. Also on the C/C++ build page, at the top hit 'Manage Configurations' then 'New'. Set the name to something like 'Mushy Layer Debug', and copy settings from the default configuration. Then select the debug configuration on the C/C++ build page, untick 'Use default build command' and enter `make DEBUG=TRUE OPT=FALSE`. When you use this build profile, you will force the compiler to use the DEBUG flags and not use the OPT flags. This creates and executable suitable for debugging with.

7.  Check everything works so far. Build the project using both the configurations we've created via Project->Build all. Use Project->Build Configurations-Set active to change between build configurations.

8.  Create a debug configuration. Now Eclipse knows how to build your project, but you haven't explicitly told it where the executable file is located, or what options it should be run with. From the top menu, go Run->Debug configurations. Click the little 'New Debug configuration' icon (top left), give it a name, and for C/C++ application navigate to the DEBUG version of your compiled codes e.g. `mushyLayer2d.Linux.64.mpiCC.gfortran.DEBUG.MPI.ex`. For Build configuration, choose the debug configuration you just created. Now go to the Arguments tab and enter the name of the inputs file you'll want to use for your debugging (you'll probably find you change this a lot). Make sure the Working Directory at the bottom is set to `${workspace_loc:MushyLayer}/execSubcycle`. Now go to the Environment tab and set
`CH_TIMER=1` - this makes sure chombo produces timing output from simulations, which can be useful.
`VISITHOME=/path/to/visit2.12.0/src` - this allows Eclipse to find your installation of Visit (if you've installed Visit) so it can be used for debugging. (your visit installation is likely in a different place)
`LD_LIBRARY_PATH=/usr/local/hdf5-1.8.21p-v18/lib:` - make sure Eclipse can find HDF5. (your hdf5 installation is likely in a different place)
Finally, go to the Debugger tab and uncheck 'Stop on startup at: main' as otherwise this will drive you crazy.

You should now be all set up for debugging. Here is an example. Open up `srcSubcycle/AMRLevelMushyLayer.cpp`, find the `advance(...)` method (the outline panel on the right hand side is useful for this) and find the code
```cpp
if (solvingFullDarcyBrinkman())
  {
    // This fills all the ghost cells of advectionSourceTerm
    computeAdvectionVelSourceTerm(advectionSourceTerm);

    //    Real vel_centering = 0.5;
    //    ppMain.query("adv_vel_centering", vel_centering);
    computeAdvectionVelocities(advectionSourceTerm, m_adv_vel_centering);
  }
```
Insert a breakpoint on the `if...` line (Run->Toggle breakpoint) and then go Run->Debug. The code should run until this line then stop. Let's have a look at the source term for the advection velocity solve. In the console at the bottom of the window, enter 
`call viewLevel(&advectionSourceTerm)`
and an instance of Visit should open showing you the source term field on this level of refinement. Now step through the lines (F6) until you get past this `if/else` clause, by which point you will have calculated the advection velocity. We can look at this by entering
`call viewFluxLevel(&m_advVel)`
Other useful options available for viewing data include

*   `viewFAB(&fab)` for FArrayBoxes (LevelData's are composed of FArrayBoxes). If you have a FluxBox `fluxbox` instead, you can look at each component in turn with `viewFAB(&fluxbox[0])`, `viewFAB(&fluxbox[1])` etc.
*   `p m_time` return the value of any variable (you can also hover over a variable in the editor window).

There are probably more, and these are documented in the ChomboDesign.pdf document.
