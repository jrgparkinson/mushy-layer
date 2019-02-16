## Installing everything
This is an attempt at listing the steps required to install Chombo from a blank slate of a unix system. It likely contains various holes where I forgot I had done something to get it working - please send any additions or corrections to [jamie.parkinson@gmail.com], thanks.

# Fortran
gfortran works well on ubuntu, my installation looks like

```
parkinsonjl@atmlxlap005:~$ gfortran --version
GNU Fortran (Ubuntu 4.8.4-2ubuntu1~14.04.4) 4.8.4
```

# C++ compiler
Make sure you have a sufficient compiler version, as detailed in the ChomboDesign doc. Mine looks like

```
parkinsonjl@atmlxlap005:~$ mpicc --version
gcc (Ubuntu 4.8.4-2ubuntu1~14.04.4) 4.8.4
```

# HDF5
The latest version of hdf5 (1.10.xx) doesn't seemt to work with Chombo at the moment. You can download the source code for version 1.8.21 from [https://portal.hdfgroup.org/display/support/HDF5+1.8.21#files], which does work. There are decent install instructions in release_docs/INSTALL_parallel. Make sure you use the `--enable-parallel` option. I used the command

`./configure --enable-parallel --prefix=/home/parkinsonj/soft/hdf5/`

and it seems to work fine with Chombo.  Once installed make sure you set the evironment variable 

`LD_LIBRARY_PATH=/path/to/hdf5` e.g. `LD_LIBRARY_PATH=/home/parkinsonj/soft/hdf5/`.

Do this in your .bashrc or whatever gets called whenever you open a terminal so it's always there.


# Chombo
Chombo download instructions are available at [https://anag-repo.lbl.gov/chombo-3.2/access.html]. You need to register for an account, then you can get the code using svn:

```
svn --username username co https://anag-repo.lbl.gov/svn/Chombo/release/3.2 /your/path/to/chombo/
```

I found the following options in Make.defs.local worked:

```
DIM           = 2            # 3 also works
DEBUG         = FALSE        # true also works
OPT           = TRUE         # false also works
FC            = gfortran
MPI           = TRUE
MPICXX        = mpiCC
USE_HDF       = TRUE
HDFMPIINCFLAGS= -I/usr/local/hdf5-1.8.21p-v18/include
HDFMPILIBFLAGS= -L/usr/local/hdf5-1.8.21p-v18/lib -lhdf5 -lz 
flibflags     = -lblas -llapack -lgfortran 
LAPACKLIBS = -L/usr/lib -L/usr/lib/lapack -L/usr/lib/libblas -llapack  -lblas
```

If the blas or lapack packages don't exist, you may need to run something like `sudo apt-get install libblas-dev liblapack-dev`. These are tools for mainly doing linear algebra (inverting matrices etc.) very quickly.

# Visit
Visit is very useful for quickly visualising data. Download Visit from [https://wci.llnl.gov/simulation/computer-codes/visit/source]. There are good install instructions supplied with the source code/on their website. Add your visit install location to your PATH.

# doxygen
Doxygen let's you build documentation from the source code (and the comments within it). It can be downloaded from [http://www.doxygen.nl/download.html]. Then go to the `/path/to/mushy-layer/doc` directory and run `doxygen` in a terminal - it should pick up the `Doxyfile` already there, which specifies how to build the documentation.

