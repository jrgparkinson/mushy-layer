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


## Running code on AOPP servers
There is an example SLURM batch script in `execSubcycle/slurmExample`.
There is a batch script for building all the mushy layer code (with the relevant module dependencies) in `execSubcycle/buildMushyLayer.sh`, which can be modified for your needs.


## Time saving tips for AOPP
1.  **Automating tasks on login**. After logging in to the machine where you run code from (e.g. atmlxmaster), add commands to your ~/.login file so you don't have to run them manually each time you login, e.g.
```bash
module load visit/2.13.2 # gives you access to the visit command, which will let you visualise hdf5 files
module load hdf5/1.8.18v18__parallel # load hdf5 for data I/O
module load chombo/3.2__hdf5v18 # load chombo (though you may find it easier to just compile your own copy locally)

# Set default squeue format to display longer job names
setenv SQUEUE_FORMAT "%.18i %.9P %.25j %.8u %.2t %.10M %.6D %.20R"

# Make sure we produce time.table files in chombo
setenv CH_TIMER "true"
```

2. **Mounting remote drives**. (Currently only tested on linux machines). You can mount remote drives like shared storage or your home directory over SSH, so that you can access them from your default file browser. On your local computer, first create a blank directory into which you will mount the drive, e.g.
```console
$ cd ~
$ mkdir mnt
$ cd mnt
$ mkdir shared_storage
```
Now setup an SSH alias to tunnel through the gateway server, by adding the following (change your username!) to `~/.ssh/config` (create the file if it doesn't already exist):
```
Host atmlxmaster
User          parkinsonj
HostName      atmlxmaster
ProxyCommand  ssh parkinsonj@ssh-gateway.physics.ox.ac.uk nc %h %p %r 2> /dev/null
```
then run the command (change your shared storage folder unless you want to mount mine!)
```console
$ sshfs -o idmap=user atmlxmaster:/network/group/aopp/oceans/AW002_PARKINSON_MUSH ~/mnt/shared_storage
```
you will be prompted for your SSH password (possibly twice). Thereafter you can access all your files from your local machine by going to
```
cd ~/mnt/shared_storage
```
You can also mount your home drive in the same way i.e.
```console
$ sshfs -o idmap=user atmlxmaster:/home/parkinsonj/  ~/mnt/home
```
