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

# Dimensions to build in (usually 2)
D=2
echo "Build with DIM=$D"

#cd ~/mushy-layer/execSubcycle/
cd $MUSHY_LAYER_DIR/execSubcycle/
make all DIM=$D | grep --color -E 'is up to date'

cd $MUSHY_LAYER_DIR/setupNewRun/
# make clean
make all DIM=$D | grep --color -E 'is up to date'

cd $MUSHY_LAYER_DIR/postProcess/
make all DIM=$D | grep --color -E 'is up to date'

module unload gcc/4.8.3
module load gcc/5.2.0

