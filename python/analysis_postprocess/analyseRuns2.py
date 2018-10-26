from MushyLayerRun import MushyLayerRun
import os
#import numpy
import math
from subprocess import Popen
from mushyLayerRunUtils import constructRunName, readInputs, writeParamsFile, timeSinceFolderUpdated
from ChkFile import ChkFile
import sys, getopt
from os import listdir
from os.path import isfile, join
import shutil
import time
import re
from runPostProcessFolder import runPostProcessFolder
import csv

#thisDir = os.getcwd()
#sys.path.append('/path/to/application/app/folder')
currentDir = os.getcwd()
MushyLayerDir = currentDir.replace('python', '')

katzUnits = True

# Code to loop through all runs in specified folder and do something with them

# Remote:
#baseDir = '/network/group/aopp/oceans/AW002_PARKINSON_MUSH/optimalStates-highRes-newSoluteFlux/' 
baseDir = '/network/group/aopp/oceans/AW002_PARKINSON_MUSH/optimalStates-highRes-new/' 


# Local
baseDir = '/home/parkinsonjl/mnt/sharedStorage/optimalStates-highRes-new/'


# 4 processors should speed things up a little
num_proc = 1

foldersToProcess = []
emptyFolders = []
allFolders = []
# Get all folders that aren't 'steady'

for f in next(os.walk(baseDir))[1]:
    #print(f[-2:])
    allFolders.append(f)
    #if f[-2:] == '-0':
    identifier = 'steady'
    if not f[-len(identifier):] == identifier:    
        foldersToProcess.append(f)
    else:
        emptyFolders.append(f)

     
        
        
for folder in foldersToProcess:
    fullFolderLoc = os.path.join(baseDir, folder)
        
    inputsPath = os.path.join(fullFolderLoc, 'inputs')
    if not os.path.exists(inputsPath):
        print('No inputs file found')
        continue
    
    inputs = readInputs(inputsPath)
    
    print(folder)
    
    # If run is still going, do nothing
    stillRunning = True
    
    if os.path.exists(os.path.join(fullFolderLoc, 'time.table.0')):
        stillRunning = False
        
    
    # Check last time a file was updated    
    minutesSinceFolderUpdated = timeSinceFolderUpdated(fullFolderLoc)
    if minutesSinceFolderUpdated > 120:
        stillRunning = False
    
    
    if stillRunning:
        print('  Still running, skip')
        continue
    
    # Just get the diagnostics to decide if this is actually steady
    steadyState = False
    
    # Load runs for the csv file
    with open(os.path.join(fullFolderLoc, 'diagnostics.out'), 'rb') as f:
        reader = csv.DictReader(f)
        
        Fs = []
        dts = []
        for row in reader:
            Fs.append(float(row['Fs_vertical_av']))
            dts.append(float(row['dt']))
            
        if len(Fs) > 2:
            dt = dts[-1] # guess
            dFdt = abs((Fs[-1] - Fs[-2])/Fs[-1])/dt
            print('  (delta flux) / flux = ' + str(dFdt))
            if dFdt < 1e-3:
                steadyState = True
    
    # Now rename appropriately
    parts = folder.split('-')
    oldEnding = parts[1]
    
    if steadyState:
        newEnding = 'steady'
    else:
        newEnding = 'unsteady'
    
    newName = folder[:-len(oldEnding)] + newEnding
        
    print('  Rename to  ' + newName)
        
    try:
        os.rename(fullFolderLoc, join(baseDir, newName))
    except OSError as (errno, strerror):
        print "OS error({0}): {1}".format(errno, strerror)
   
    
print('Finished analysing all folders')
