# Example problems

This directory contains a range of example problems. In each folder there should be, as a bare minimum, an inputs file which you can run like

```console
$ cd /mushy-layer/examples/sinusoidalcooling/
$ ../../execSubcycle/mushyLayer2d.Linux.64.mpiCC.gfortran.OPT.MPI.ex  inputs 
```

There is not a lot of other documentation yet - sorry. You may find the python script at `/plotting/boundaryConditionsFigure.py` useful for visualising the boundary conditions and physical parameters in each of these examples.

## 3D
Most of these inputs are designed for 2D simulations. To use them for 3D simulations is fairly trivial, just remember to updated the boundary conditions, `bc.XXX` and `main.periodic_bc` so that the final element is always the vertical direction. 
