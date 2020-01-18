# Running the Mushy Layer code
## Introduction
This document explains how to use the mushy layer code for various physical applications. 

It assumes you already have a compiled version of the code which you can run.

## Basics
There are over a hundred options used in the code, which are stored on the `AMRLevelMushyLayer` class in one of two places:

1. The member variable `AMRLevelMushyLayer.m_parameters` which is a `MushyLayerParams` object. This is where most physical parameters (e.g. Rayleigh number, Lewis number) are stored.
2. The member variable `AMRLevelMushyLayer.m_opt`, which is a struct of type `MushyLayerOptions` and contains options on how to solve the equations (e.g. solver thresholds, how to do projection, and lots of debugging options).

Inputs are loaded in the code via the `ParmParse` class, e.g.

```cpp
Real multigrid_threshold = 1e-10;
ParmParse pp("main");
pp.query("mg_thresh", multigrid_threshold)
```
which would load some inputs like
```
main.mg_thresh=1e-15
```

We do most of this in either `/mushy-layer/srcSubcycle/MushyLayerSubcycleUtils.cpp` or `/mushy-layer/src/MushyLayerParams.cpp`.

Most of these options/parameters can just be left at their default values. The inputs file at `/mushy-layer/examples/meltponds/inputsGrowSeaIce` illustrates the main options you may wish to change, which are described in more detail below under some roughly suitable headings.

## General options
`main.verbosity=2` larger verbosity means there will be more text output produced

`main.output_folder=.` folder to save plot/checkpoint files to a period `.` means the directory from which the code was executed

`main.plot_prefix=plt` prefix for plot files, which will then be like `plt001234.2d.hdf5`

`main.chk_prefix=chk` prefix for checkpoint files

`main.plot_interval=n`  produce plot files every `n` steps (-1 = don't use this option)

`main.checkpoint_interval=n` produce checkpoints files every `n` steps

`main.plot_period=0.005` time interval at which to write out plot files

`main.debug=false`  set to true to write more fields to the plot files

## Timestepping
`main.cfl=0.1` max allowed CFL number

`main.initial_cfl=2e-05` start simulations with a smaller CFL number. The timestep  increases thereafter as determined by `main.max_dt_growth`

`main.max_dt=0.1` max allowed timestep

`main.max_dt_growth=1.05` max fractional increase in $\Delta t$ allowed from one timestep to the next

`main.max_step=10000` stop simulations after this many timesteps

`main.max_time=10.0` stop simulations once this time has been reached

`main.min_time=0.1` ensure simulations run until at least this time (even if a steady state condition is reached)

`main.dt_tolerance_factor=1.01` Set the factor by which the current dt must exceed the new (max) dt for time subcycling to occur (i.e., reduction of the current dt by powers of 2).

## Domain setup
`main.num_cells=Nx Ny Nz` e.g.  `main.num_cells=64 128` size of the grid on the coarsest mesh. Make sure each dimensions is a power of 2 unless you want the multigrid to be slow. Only need to specify 2 dimensions if the code is compiled for 2 dimensions. The final dimension is always the 'vertical' with respect to gravity.

`main.domain_height=1.0` domain height (width is then calculated from the number of grid cells specified)

`main.periodic_bc=1 0` whether to enforce periodic boundaries (1) or not (0) in each dimension

`main.max_grid_size=256` split domain into boxes (for sending to different processors) of max size in any dimension given by this parameter.

## Initial and boundary conditions
See the [separate document](InitialAndBoundaryConditions.md).


## Dimensionless parameters
`main.nondimensionalisation=0` choose nondimensionalisation, different options are:
0. Diffusive timescale, advective velocity scale
1. Darcy timescales, advective velocity scale
2. Darcy timescale, darcy velocity scale
3. Advective timescale, darcy velocity scale
4. Buoyancy timescale, advective velocity scale

`parameters.problem_type=0` can specify different problem types, which can tell the code to solve some of the equations differently. Lots of the options aren't really used anymore (and maybe don't work). Ones which definitely should work are in bold.
0. *mushyLayer*
1. burgers equation
2. burgers equation on a periodic domain
3. poiseuille flow
4. diffusion
5. *solidificationNoFlow*
6. corner flow
7. sidewall heating
8. vertical heating (with Darcy's equation)
9. rayleigh-benard (navier-stokes)
10. a solute flux test
11. a reflux test
12. a test for solving with zero porosity
13. simulating a melting ice block
14. *convectionMixedPorous*
15. vortex pair (navier-stokes benchmark problem)

`parameters.compositionRatio=1.18` composition ratio, should be >= 1.0.


`parameters.rayleighComp=500` rayleigh number for compositional differences. When just solving Darcy's equation, this is the mushy layer number. Otherwise, it's the fluid rayleigh number.

`parameters.rayleighTemp=25.0` rayleigh number for temperature.

`parameters.stefan=5.0` stefan number

`parameters.K=1` ratio of solid/liquid heat diffusivity

`parameters.heatConductivityRatio=1` ratio of solid/liquid heat conductivity

`parameters.specificHeatRatio=1` ratio of solid.liquid specific heat

`parameters.lewis=100` lewis number (liquid heat diffusivity/liquid salt diffusivity)

`parameters.darcy=0.0` darcy number. Set to 0 to solve Darcy's equation rather than Darcy-Brinkman

`parameters.prandtl=0` prandtl number 

`parameters.heleShaw=true` constrain permeabilities by a Hele-Shaw cell

`parameters.nonDimReluctance=0.1` inverse of the dimensionless Hele-Shaw permeabilty. Smaller means wider hele-shaw gap width. Alternatively, you can now specify the dimensionless Hele-Shaw cell permeability e.g. `parameters.heleShawPermeability=10.0`
 

`parameters.nonDimVel=0.0` dimensionless frame advection velocity
`parameters.permeabilityFunction=2` mathematical form of the permeability function. Options include:
0. pureFluid (permeability=1)
1. cubic (permeability=porosity^3)
2.   kozeny-carman (permeability = porosity^3 / (1-porosity)^2
3.   log function (permeability = - porosity^2 * log(1-porosity) )
4.   spatially dependent permeability with a channel
5.   permeability = porosity


`parameters.waterDistributionCoeff=1e-5` water distribution coefficient

## Phase diagram
`parameters.eutecticComposition=230`

`parameters.eutecticTemp=-23`

`parameters.initialComposition=30`

`parameters.liquidusSlope=-0.1`

## AMR options
`main.max_level=0` max allowed level for AMR simulations.

`main.block_factor=8` all boxes must be divisible by at least this length in all dimensions. 

`main.use_subcycling=1` whether or not to use subcycling in time

`main.ref_ratio=2 2 2` refinement ratios between successive levels 

`main.refine_thresh=0.1` refinement threshold

`main.regrid_interval=8 8 8`  regrid interval on each level 

`main.tag_buffer_size=4` number of cells by which to grow the set of cells tagged for refinement, before creating the new grids

`main.fill_ratio=0.75` should be between 0 and 1. A value of 1 means the grids generated on refined levels will be as tightly grouped as possible to the cells tagged for refinement. A smaller value means that Chombo will cover a larger area on refined levels with grid cells. 

`main.grid_buffer_size=0 0 0` this is the 'padding' between grids on different levels of refinement

`projection.eta=0.0` Freestream correction coefficient. Should be less than 1 for stability, but close to 1 for accuracy.


## Projection
The projection solver has some tolerance specified by `projection.solverTol`. If the unprojected velocity has some divergence $d$, then the projected velocity will have some divergence of order $d \times $`projection.solverTol`. In reality, we care more about the absolute value of the final divergence than how much it is has been reduced. Therefore we introduce an option to adapt the solver tolerance to achieve a particular final divergence:

`projection.adaptSolverParamsDivU = 1`  turns on this option (set 0, or don't set at all, to not use this option)
`projection.divergence_tolerance = 1e-9`  absolute divergence we wish to aim for

If you are struggling to achieve this final divergence, you can try increasing either the number of multigrid smooths or number of multigrid iterations:

`projection.numSmoothUp=64`
`projection.maxIter=40`

Additionally, you can try using the pressure from the previous timestep to remove a significant ammount of the divergence before projection:

`projection.useIncrementalPressure=true`

i.e. after computing the unprojected velocity $\mathbf{U}^*$, then subtract off the old pressure gradient to find $\mathbf{U}^* - \chi \nabla p$ before projecting this to find a (hopefully) small extra pressure correction. This can significantly speed up the solve.

