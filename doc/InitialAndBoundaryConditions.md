# Initial and boundary conditions
This document describes how to implement the initial conditions and boundary conditions.

## Initial conditions


## Boundary conditions

Boundary conditions are defined in two ways. Firstly we define the type of boundary condition to apply to each side, then we define the value to use (where needed) at each side. For scalars, options are:

0. Dirichlet
1. Neumann
2. Inflow/outflow

For velocity, there are lots of options, some of the most important being

0. Solid wall
2. Inflow/outflow - pressure = 0, velocity determined by the projection
3. Inflow/outflow with purely normal flow enforced (a bit dodgy)
6. Reflection
9. Pressure head

See BCUtil/PhysBCUtil.h for a full list of possible options.

We set boundary conditions on the 'high' side of the domain in each dimension, then the 'low' side of the domain in each dimension e.g.

`bc.scalarHi=1 0`  scalar BCs are neumann on the right boundary, dirichlet on the top

`bc.scalarLo=1 2` scalar BCs are neumann on the left boundary, inflow/outflow at the bottom

`bc.velHi=6 0` 

`bc.velLo=6 2`

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

