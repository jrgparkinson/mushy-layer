## Guide to installing Chombo on AOPP servers

1. Login into the server e.g. `$ ssh -X gyre1`
2. **10 minutes**. Get the latest chombo version. 
Whilst the [Chombo website](https://anag-repo.lbl.gov/chombo-3.2/access.html) will tell you to get the repository located at 
`https://anag-repo.lbl.gov/svn/Chombo/release/3.2`,
you should actually get the code from 'Chombo Trunk' which is like a development branch and contains code needed for the Mushy Layer program. This is located at
`https://anag-repo.lbl.gov/svn/Chombo/trunk/`, 
so, using `svn`, you should run something like (replace 'jparkinson' with your username)
```console
$ svn --username jparkinson co https://anag-repo.lbl.gov/svn/Chombo/trunk/ ~/Chombo
```
3. Copy the makefile from `/mushy-layer/doc/Make.defs.AOPP` (in the same directory as this readme) to `Chombo/lib/mk/Make.defs.local` e.g.
```console
$ cp ~/mushy-layer/doc/Make.defs.AOPP ~/Chombo/lib/mk/Make.defs.local
```
4. Add the following to your ~/.login script (note this is for csh/tcsh shells, if you're using a bash shell then you'll need to use [different commands](https://web.fe.up.pt/~jmcruz/etc/unix/sh-vs-csh.html) to set the environmental variables).
```bash
module load chombo/3.2__hdf5v18
setenv CH_TIMER "true"
```
And to your .bashrc
```bash
setenv CHOMBO_HOME /path/to/Chombo/
```
5. Logout and log back in.
6. **10 minutes**. Go to the Chombo code and start compiling
```console
$ cd $CHOMBO_HOME/lib/
$ make lib
```
Technically you should now build and run all the test problems, but if you're feeling lucky you can skip straight to step 7.
```console
$ make test
$ make run | grep 'finished with'
```
Everything test should report 'finished with status 0'
7. **5 minutes**. Build and run one of the released examples as a test
```console
$ cd $CHOMBO_HOME/releasedExamples/AMRGodunov/execAdvectDiffuse/
$ make all
$ ./amrGodunov2d.Linux.64.mpiCC.gfortran.OPT.MPI.ex diffuse.inputs 
```
The output files produced are most easily opened using the [Visit software](https://wci.llnl.gov/simulation/computer-codes/visit), which is not currently installed on the AOPP servers (though it may be soon - I have asked for it).


