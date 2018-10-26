# Script to find optimal states, varying the domain
# width to find the optimal flux for fixed Ra, then increasing Ra and 
# varying the width again

import os
import sys
parDir = os.path.abspath(os.pardir)
pythonDir = os.path.join(parDir, 'python')
sys.path.append(pythonDir)
import numpy

import csv

from MushyLayerRun import MushyLayerRun

import math
from subprocess import Popen
from mushyLayerRunUtils import constructRunName, readInputs, writeParamsFile
from ChkFile import ChkFile


import getopt

def main(argv):
    
    folder = 'middleton/CR1.179RaC80Le200KozenyPermeabilitypts256-S3.5-TB-10.0-R0.13-restart/'
    prefix = 'plt'
    start = 24150
    interval = 100
    end = 200750
    
    try:
       opts, args = getopt.getopt(argv,"f:s:i:e:p:")
    except getopt.GetoptError:
       print 'deletePltFiles.py -n <num processors> -f <folder> -H <domain height>'
       sys.exit(2)
    for opt, arg in opts:
        if opt in ("-s"):
            start = int(arg)
        elif opt in ("-f"):
            base_folder = str(arg)
        elif opt in ("-i"):
            interval = int(arg)
        elif opt in ("-e"):
            end = int(arg)
        elif opt in ("-p"):
            prefix = str(arg)
          
          
   
    
    base_folder = '/home/parkinsonjl/mnt/sharedStorage/'
    #base_folder = '/network/group/aopp/oceans/AW002_PARKINSON_MUSH/'
    base_dir = os.path.join(base_folder, folder)
    
    frames = numpy.arange(start,end,interval)
    
    for frame in frames:
        name = "%06d" % frame
        filename = base_dir + prefix + name + ".2d.hdf5"
        if os.path.exists(filename):
            
            print("Delete " + filename)
            try:
                os.remove(filename)
            except e:
                print("Delete failed")
                #pass
    
           

if __name__ == "__main__":
    print('Running script')
    main(sys.argv[1:])
                
