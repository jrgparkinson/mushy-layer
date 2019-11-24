# Regression tests
This directory contains regression tests. The idea is to run a simulation and check that the output hasn't changed (hasn't 'regressed'). Assuming the original output was correct, the program still works.

We can compare both the text output (e.g `pout.0` or `diagnostics.csv`) and the data output (`.hdf5` files). Each subdirectory contains a test, which at a bare minimum consist of the following files:

1. `inputs`: the input parameters for the simulation
2. `regression.json`: the options for running the test (number of dimensions, processors etc.). See an existing file for a template.
3. At least one 'expected' file e.g. `pout.0.expected` or `plt000010.2d.hdf5.expected` which will be compared with the computed output files (less the `.expected` suffix) to check for any deviations.

# Prerequisites
1. Compile the mushy layer program. Needs to be done separately in 2D and 3D for the relevant tests to work.
2. Compile the Chombo Compare utility (at `$CHOMBO_HOME/lib/utils/ChomboCompare`) for comparing hdf5 files, and ensure the `CHOMBO_HOME` environmental variable is set.

# Usage
Run all the tests using the python file `run_regression_tests.py` e.g.

```bash
~/mushy-layer/test/regression$ python run_regression_tests.py
```

Check the output to the terminal, which will report which tests fail and which succeed, along with further information on test failures.
