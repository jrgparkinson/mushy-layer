# Initial and boundary conditions
This document describes how to implement the initial conditions and boundary conditions.

## Initial conditions

There are two options for setting the initial data in a simulation. 

### Restart file
Firstly, you may specify a previous data file from which the simulation will restart, via the input option `main.restart_file`, e.g.

```
main.restart_file = /full/path/to/chk99999.2d.hdf5
```

Note that the restart file must be a checkpoint file (as opposed to a plot file), and must have the same data structure as your simulation (i.e. the same number of potential levels of refinement).

### Specify in code
The second option is to declare the initial conditions analytically within the code. This is all controlled by  the function `AMRLevelMushyLayer::initialData()`, which can be found in `/srcSubcycle/AMRLevelMushyLayerInit.cpp`. For some types of problems (declared by `parameters.problem_type`) there are specific initial conditions defined. For the standard problem type, `parameters.problem_type=0` corresponding to mushy layer simulations, the initial data is defined in `AMRLevelMushyLayer::initialDataMushyLayer()` where there are lots of possible variations but typically we just set the enthalpy and bulk concentration to their boundary values on the bottom boundary (i.e. the vertical direction, lower face). This is designed for simulations where ice is grown from a cold upper boundary, and hence the bottom boundary corresponds to the warm ocean. 

You may overwrite the default initial conditions using the `main.initData` input option.  Currently there are just two options:

`main.initData=0`: default, set H and C to their values at the bottom boundary (within `AMRLevelMushyLayer::initialDataDefault()`)
`main.initData=1`: create a porous mushy hole in the center of the domain, for one of the benchmark problems (within `AMRLevelMushyLayer::initialDataPorousHole()`) 

Here is an example of how you would add a new initial condition that establishes a linear enthalpy gradient and constant bulk concentration across the domain:

1. Create a new initial data function `AMRLevelMushyLayer::initialDataLinearGradient()`, e.g.

```c++
void AMRLevelMushyLayer::initialDataLinearGradient()
{

  // Iterate over the different boxes on this level of refinement
  DataIterator dit = m_grids.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
  {
    // Iterate over the cells within each box
    BoxIterator bit((*m_scalarNew[0])[dit].box());
    for (bit.reset(); bit.ok(); ++bit)
    {
      // This is the location of the grid cell in indexing space
      IntVect iv = bit();
      
      // Get the location of this grid cell in x-y space
      RealVect loc;
      getLocation(iv, loc, m_dx);

      // Set the enthalpy = stefan number - 0.1*z, where z is the vertical co-ordinate
      // note that loc[1] is the vertical co-ordinate
      (*m_scalarNew[ScalarVars::m_enthalpy])[dit](iv) = m_parameters.stefan - 0.1*loc[1];
      
      // Set the bulk concentration
      (*m_scalarNew[ScalarVars::m_bulkConcentration])[dit](iv) = -1.0;
    }
  }


}
```

2. Also add this function to the header file `AMRLevelMushyLayer.h`, i.e.
```c++
/// Mushy layer with a porous hole in the middle
  void initialDataDefault();
  
  /// New function
  void initialDataLinearGradient();
```

3. Add the new function to the list of possible options in `AMRLevelMushyLayerInit.cpp`,

```c++
switch(m_opt.customInitData)
    {
      case 1:
        initialDataPorousHole();
        break;
      case 2:   // if main.initData = 2, use this function:
         initialDataLinearGradient();
         break; 
      default:
        initialDataDefault();
    }
```



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

