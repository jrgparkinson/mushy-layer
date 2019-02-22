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
*main.verbosity=2
*main.output_folder=.
*main.plot_prefix=plt
*main.plot_interval=-1  # step interval to produce plot files (-1 = don't use this option)
*main.checkpoint_interval=1000
*main.debug=false  # set to true to produce more output
*main.plot_period=0.005  # time intervals to produce plot files
*main.chk_prefix=chk

# Timestepping
*main.cfl=0.1
*main.initial_cfl=2e-05

*main.max_dt=0.1
*main.max_dt_growth=1.05

*main.max_step=10000
*main.max_time=10.0
*main.min_time=0.1
*main.dt_tolerance_factor=1.01

# Domain setup
*main.domain_height=1.0
*main.num_cells=64 128
*main.periodic_bc=1 0 0
*main.max_grid_size=256

# Boundary conditions
Boundary condition values:
*bc.bulkConcentrationHiVal=-1 -1
*bc.bulkConcentrationLoVal=-1 -1
*bc.enthalpyHiVal=0 0.9
*bc.enthalpyLoVal=0.0 6.03

# Boundary conditions types. 
For scalars: 0 = dirichlet, 1 = neumann, 2 = reflection
For velocity: 0 = solid wall, 3 = inflow/outflow, 6 = reflection,  9 = pressure head
See BCUtil/PhysBCUtil.(h/cpp) for more info/options
bc.scalarHi=1 0 0
bc.scalarLo=1 2 2
bc.velHi=6 0 0
bc.velLo=6 3 0
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


