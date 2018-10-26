from MushyLayerRun import MushyLayerRun
import os
import math
from subprocess import Popen
from mushyLayerRunUtils import constructRunName, readInputs, writeParamsFile
from ChkFile import ChkFile
import sys, getopt
from os import listdir
from os.path import isfile, join


def runPostProcessFolder(folder, extraParams, execDir, outfile):
    
    # 1) Get final checkpoint file/checkpoint prefix
    inputs = readInputs(os.path.join(folder, 'inputs'))
    chkPrefix = inputs['main.chk_prefix']
    
    # 2) Make inputs for post processing  
    if not folder[-1] == '/':
        folder = folder + '/'
         
    postProcessInputs = dict()
    
    # Where are the files to analyse located?
    postProcessInputs['in_folder'] = folder
    postProcessInputs['prefix']  = chkPrefix # e.g. chkCR1.1RaC125Le200ChiCubedPermeabilitypts64-276000.2d.hdf5
    postProcessInputs['out_file'] = os.path.join(folder, outfile)
    
    for key in extraParams.keys():
        postProcessInputs[key] = extraParams[key]
        
    # Write this file out somewhere
    inputsLoc = os.path.join(folder, 'postProcess.inputs')
    writeParamsFile(inputsLoc, postProcessInputs)
    
    # 3) Execute post process code
    execLoc = join(execDir, 'process2d.Linux.64.mpiCC.gfortran.OPTHIGH.MPI.ex')
    cmd = 'cd ' + execDir + '; mpirun -np 1 '  + execLoc + ' ' + inputsLoc 
    os.system(cmd)

