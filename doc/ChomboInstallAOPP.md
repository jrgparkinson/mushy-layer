## Guide to installing Chombo on AOPP servers

1. Login into the server e.g. `$ ssh -X gyre1`
2. Pull down the latest chombo version. 
Whilst the [Chombo website](https://anag-repo.lbl.gov/chombo-3.2/access.html) will tell you to get the repository located at 
`https://anag-repo.lbl.gov/svn/Chombo/release/3.2`,
you should actually get the code from 'Chombo Trunk' which is like a development branch and contains code needed for the Mushy Layer program. This is located at
`https://anag-repo.lbl.gov/svn/Chombo/trunk/`, 
so, using `svn`, you should run something like (replace 'jparkinson' with your username)
```console
$ svn --username jparkinson co https://anag-repo.lbl.gov/svn/Chombo/trunk/ Chombo
```
3. Copy the make file from /mushy-layer/doc/Make.defs.AOPP (in the same directory as this readme) to Chombo/lib/mk/Make.defs.local
4. Add the following to your ~/.login script 
```bash
module load chombo/3.2__hdf5v18
setenv CH_TIMER "true"
```
And to your .bashrc
```bash
setenv CHOMBO_HOME /path/to/Chombo/
```
5. Logout and log back in.
6. Go to the Chombo code and start compiling
```console
$ cd $CHOMBO_HOME/lib/
$ make lib
```
7. Build and run one of the released examples as a test
```console
$ cd $CHOMBO_HOME/releasedExamples/AMRPoisson/
```



