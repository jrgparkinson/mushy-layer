#!/bin/bash -l

# Batch script to build all mushy layer applications on AOPP servers
# Makes sure you use gcc 4.8.3 to build the code 
# (I load gcc 5.2.0 by default else the latest MATLAB crashes)
# Assumes you have set the MUSHY_LAYER_DIR environmental variable
# Add an alias to your bashrc (or equivalent): 
#    alias buildmush 'MUSHY_LAYER_DIR/execSubcycle/buildMushyLayer.sh'
# in order to access this script by just running 
#    $ buildmush 
# from the command line

module unload gcc/5.2.0
module load gcc/4.8.3

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

echo "Building $DIR"
#cd ~/mushy-layer/execSubcycle/
cd DIR
make all

cd ../setupNewRun/
# make clean
make all

cd ../postProcess/
make all

module unload gcc/4.8.3
module load gcc/5.2.0

