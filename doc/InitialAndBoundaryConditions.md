# Initial and boundary conditions
This document describes how to implement the initial conditions and boundary conditions.

The python script in `/plotting/boundaryConditionsFigure.py` is designed to create a visualisation of the boundary conditions specified in an inputs file, which you may find useful. 

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

Boundary conditions are defined in two ways. Firstly we define the type of boundary condition to apply to each side, then we define the value to use (where needed) at each side. For scalars, the simple options are:

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

We must now specify the values of the various thermodynamic fields on each boundary. In older versions of the code, the only way to do this was to specify the enthalpy and bulk salinity, e.g. 

```
bc.bulkConcentrationHiVal=-1 -1
bc.bulkConcentrationLoVal=-1 -1
bc.enthalpyHiVal=0 0.9
bc.enthalpyLoVal=0.0 6.03
```

From which the temperature, porosity, liquid salinity and solid salinity are computed using the phase diagram.

However, this is not always particularly convenient as it doesn't give precise control over the temperature at a boundary. If you instead want to fix other fields, this can be done directly e.g. 

```
bc.temperatureHi = 0 1
bc.temperatureHiVal = 0.5 0.0
```

would set the temperature at the top boundary to be 0.5. Other relevant thermodynamic variables can similarly be controlled e.g.

```
bc.porosityLo
bc.liquidConcentrationLo
permeabilityLo
```

### Mixed temperature BC
A notable extra boundary condition available is designed to enforce the following condition:

```
a \frac{d \theta}{dz} - F - b(\theta - \theta_{ref}) = 0
``` 

where `\theta` is the dimensionless temperature. The correct implementation of this BC can be confirmed using the test problem in `/test/diffusiveGrowth/`.

To use this BC, first set the temperature to use BC type `9` on the appropriate face(s). E.g. to use this BC on the top face,

```
bc.temperatureHi=1 9   # use no flux in x direction, and the mixed BC above in the vertical direction
```

Then define the values of a, b, F, \theta_{ref} via:

```
bc.TRefHi=0 -0.8  # \theta_{ref}
bc.aHi=0 0.5  # a
bc.bHi=0 -1.0 # b
bc.temperatureHiVal=0 -1.0 # F
```

Note that we have also defined values for these variables in the x-direction, but these are irrelevant and not used.

### Time dependent BCs
Time dependent BCs are possible. At the moment the only option implemented is a sinusoidally varying temperature. An example of this can be found in `/examples/sinusoialcooling/` like

```
bc.timeDependent=1
bc.temperatureHi=1 0
bc.temperatureHiVal=0 -999 # last value here shouldn't matter as we do time dependent bcs
main.plot_prefix=Sinusoidal

bc.sinusoidal_temperature_bc_timescale = 0.2
bc.sinusoidal_temperature_bc_amplitude = 0.2
bc.sinusoidal_temperature_bc_av = -0.8
bc.sinusoidal_temperature_bc_phase_diff = 0.0
```

These BCs are implemented in the function `PhysBCUtil::updateTimeDependentBCs()`, where should serve as a model for adding your own other time dependent BCs.

### 3D

In 3D, you must specify three values as there are three dimensions. The vertical dimension becomes the final item in the list, i.e.

```
bc.temperatureHi = 0.0 0.0 -1.0  # x, y, and z values
```

### Allowed text keywords
The latest version of the code will let you define BC types using text rather than numbers, which can make the inputs file read a little easier. This is handled by the function `MushyLayerParams::parseBCs()`. Currently supported options are:

Scalars:
* 'noflux'
* 'fixed'
* 'open'
* 'inflow'
* 'mixed' 

Velocity:
* 'solidwall' or 'noflow'
* 'inflow'
* 'outflow' or 'open'
* 'outflownormal'
* 'inflowoutflow'
* 'noshear'
* 'symmetry'
* 'pressureHead'


Example usage:

```
bc.scalarHi = noflux fixed
bc.velLo = symmetry open
```

which is equivalent to

```
bc.scalarHi = 1 0
bc.velLo = 6 2
```

Note that in one particular line, you must use either all integers or all text to define the BCs.