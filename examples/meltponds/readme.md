# Melt pond simulations
The inputs files in this directory are designed to 

1.  simulate sea ice growth
2.  add a melt pond to some sea ice
3.  simulate the melt pond evolution

## Prerequisites
To run the code you'll need to have compiled the code `/mushy-layer/execSubcycle/`, and also the code in `/mushy-layer/setupNewRun/`

```console
$ cd ../../execSubcycle/; make all
$ cd ../../setupNewRun/; make all
```

## Simulate sea ice growth
Assuming your compiled mushy layer executable is called `mushyLayer2d.Linux.64.mpiCC.gfortran.OPT.MPI.ex` (it may not be), you can start growing sea ice via

```console
$ cd ~/mushy-layer/examples/meltponds/
$ ../../execSubcycle/mushyLayer2d.Linux.64.mpiCC.gfortran.OPT.MPI.ex  inputsGrowSeaIce
```

This is a very coarse simulation, which should grow a sensible ammount of ice in a few minutes. Once you're happy that you have enough ice, cancel the execution and move on to the next step.

## Adding a melt pond
Adding a melt pond is done via a seperate program. You'll need to edit the `setupNewRunAddPond` file to ensure that the `inFile`, `run_inputs`, and `outfile` paths are correct. Then run the program via

```console
../../setupNewRun/setupnewrun2d.Linux.64.mpiCC.gfortran.OPT.MPI.ex  setupNewRunAddPond
```

This should create a `restart.2d.hdf5` file at the location specified in the `setupNewRunAddPond` inputs file.

## Simulating melt pond evolution
Make sure the parameter `main.restart_file` in `inputsEvolveMeltPond` points to the newly created file, then start simulations like before but with the new inputs file
```console
../../execSubcycle/mushyLayer2d.Linux.64.mpiCC.gfortran.OPT.MPI.ex  inputsEvolveMeltPond
```

The key difference between growing ice and evolving the melt pond simulations are:
```make
bc.enthalpyHiVal=0 6.02
parameters.pressureHead=5
bc.velHi=6 9 0 # pressure head top bc
main.plot_period=0.00002
```
We increase the enthalpy boundary at the top of the domain in the vertical direction, specify the pressure head, change velocity boundary condition at the top of the domain to use this pressure head, and make the plotting period more frequent (as we now have much stronger flow, so smaller timesteps).
