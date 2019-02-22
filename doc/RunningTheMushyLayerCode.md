## Running the Mushy Layer code
# Introduction
This document explains how to use the mushy layer code for various physical applications. 

It assumes you already have a compiled version of the code which you can run.

# Basics
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

# General options
`main.verbosity=2` larger verbosity means there will be more text output produced

`main.output_folder=.` folder to save plot/checkpoint files to a period `.` means the directory from which the code was executed

`main.plot_prefix=plt` prefix for plot files, which will then be like `plt001234.2d.hdf5`

`main.chk_prefix=chk` prefix for checkpoint files

`main.plot_interval=n`  produce plot files every `n` steps (-1 = don't use this option)

`main.checkpoint_interval=n` produce checkpoints files every `n` steps

`main.plot_period=0.005` time interval at which to write out plot files

`main.debug=false`  set to true to write more fields to the plot files

# Timestepping
`main.cfl=0.1` max allowed CFL number

`main.initial_cfl=2e-05` start simulations with a smaller CFL number. The timestep  increases thereafter as determined by `main.max_dt_growth`

`main.max_dt=0.1` max allowed timestep

`main.max_dt_growth=1.05` max fractional increase in $\Delta t$ allowed from one timestep to the next

`main.max_step=10000` stop simulations after this many timesteps

`main.max_time=10.0` stop simulations once this time has been reached

`main.min_time=0.1` ensure simulations run until at least this time (even if a steady state condition is reached)

`main.dt_tolerance_factor=1.01` Set the factor by which the current dt must exceed the new (max) dt for time subcycling to occur (i.e., reduction of the current dt by powers of 2).

# Domain setup
`main.num_cells=Nx Ny Nz` e.g.  `main.num_cells=64 128` size of the grid on the coarsest mesh. Make sure each dimensions is a power of 2 unless you want the multigrid to be slow. Only need to specify 2 dimensions if the code is compiled for 2 dimensions. The final dimension is always the 'vertical' with respect to gravity.

`main.domain_height=1.0` domain height (width is then calculated from the number of grid cells specified)

`main.periodic_bc=1 0` whether to enforce periodic boundaries (1) or not (0) in each dimension

`main.max_grid_size=256` split domain into boxes (for sending to different processors) of max size in any dimension given by this parameter.

# Boundary conditions
Boundary conditions are defined in two ways. Firstly we define the type of boundary condition to apply to each side, then we define the value to use (where needed) at each side. For scalars, options are:

0. Dirichlet
1. Neumann
2. Inflow/outflow

For vectors, there are lots of options, some of the most important being

0. Solid wall
3. Inflow/outflow
6. Reflection
9. Pressure head

See BCUtil/PhysBCUtil.h for a full list of possible options.

We set boundary conditions on the 'high' side of the domain in each dimension, then the 'low' side of the domain in each dimension e.g.

`bc.scalarHi=1 0`  scalar BCs are neumann on the right boundary, dirichlet on the top

`bc.scalarLo=1 2` scalar BCs are neumann on the left boundary, inflow/outflow at the bottom

`bc.velHi=6 0` 

`bc.velLo=6 3`

We then specify the values of the bulk concentration and enthalpy on each boundary in a similar way

`bc.bulkConcentrationHiVal=-1 -1`

`bc.bulkConcentrationLoVal=-1 -1`

`bc.enthalpyHiVal=0 0.9`

`bc.enthalpyLoVal=0.0 6.03`

Boundary values for temperature, porosity, liquid salinity etc. are computed from the specified bulk concentration and enthalpy. If needed, they can be provided explicitly using e.g
```
bc.temperatureLo
bc.porosityLo
bc.liquidConcentrationLo
permeabilityLo
```

# Dimensionless parameters
main.nondimensionalisation=0

parameters.problem_type=0 # 0 = mushy layer

parameters.compositionRatio=1.18
parameters.nonDimReluctance=0.1 # smaller means wider hele-shaw gap width
parameters.rayleighComp=500
parameters.rayleighTemp=25.0
parameters.stefan=5.0

parameters.K=1  # solid/liquid heat diffusivity ratio
parameters.heatConductivityRatio=1

parameters.lewis=100
parameters.darcy=0.0
parameters.prandtl=0

parameters.heleShaw=true
parameters.nonDimVel=0.0
parameters.permeabilityFunction=2  # 2 = carman-kozeny

parameters.specificHeatRatio=1
parameters.waterDistributionCoeff=1e-5

Specify phase diagram:
parameters.eutecticComposition=230
parameters.eutecticTemp=-23
parameters.initialComposition=30
parameters.liquidusSlope=-0.1

# AMR options
main.max_level=0
main.block_factor=8
main.use_subcycling=1
main.ref_ratio=2 2 2
main.refine_thresh=0.1
main.regrid_interval=8 8 8
main.tag_buffer_size=4
main.fill_ratio=0.75
main.grid_buffer_size=0 0 0
projection.eta=0.0 # No freestream correction, as only one level


