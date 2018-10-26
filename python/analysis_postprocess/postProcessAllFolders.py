from MushyLayerRun import MushyLayerRun
import os
#import numpy
import math
from subprocess import Popen
from mushyLayerRunUtils import constructRunName, readInputs, writeParamsFile, timeSinceFolderUpdated, timeSinceFileUpdated
from ChkFile import ChkFile
import sys, getopt
from os import listdir
from os.path import isfile, join
import shutil
import time
import re
from runPostProcessFolder import runPostProcessFolder
import csv

import getopt


def getFolders(baseDir):

	folders = []

	for root, dirs, files in os.walk(baseDir):
	    path = root.split(os.sep)
	    #print((len(path) - 1) * '---', os.path.basename(root))

	    for directory in dirs:
	    	fullDir = root + directory
	    	print(fullDir)
	    	if directory[1:3] == 'pts':

	    		print(fullDir)
	    		folders.append(fullDir)
	    #for file in files:
	    #    print(len(path) * '---', file)



	return folders


def main(argv):
    folder = 'optimalStates-highRes-new'
    forceRedo = False
    
    try:
       opts, args = getopt.getopt(argv,"rf:")
    except getopt.GetoptError:
       print 'postProcessAllFolders.py -f <base folder> (-r)'
       sys.exit(2)
    for opt, arg in opts:
        if opt in ("-f"):
            folder = str(arg)
        elif opt in ("-r"):
            forceRedo = True
        

    #local 
    #baseDir = '/home/parkinsonjl/mnt/sharedStorage/optimalStates-highRes-new/'
    # remote
    baseDir = os.path.join('/network/group/aopp/oceans/AW002_PARKINSON_MUSH/', folder)
    
    print 'Base folder: ', str(baseDir)

    #thisDir = os.getcwd()
    #sys.path.append('/path/to/application/app/folder')
    currentDir = os.getcwd()
    MushyLayerDir = currentDir.replace('python', '')

    katzUnits = True

    # Code to loop through all runs in specified folder and do something with them

    # Remote:
    #baseDir = '/network/group/aopp/oceans/AW002_PARKINSON_MUSH/optimalStates-highRes-newSoluteFlux/' 


    # Local


    # 4 processors should speed things up a little
    num_proc = 1

    foldersToProcess = []
    emptyFolders = []
    allFolders = []

    # Get all folders

    for f in next(os.walk(baseDir))[1]:
        #print(f[-2:])
        allFolders.append(f)
        #if f[-2:] == '-0':
        #identifier = 'steady'
        #if not f[-len(identifier):] == identifier:    
        #    foldersToProcess.append(f)
        #else:
        #    emptyFolders.append(f)


    allFolders = getFolders(baseDir)
         
            
            
    for folder in allFolders:
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
        
        
        if stillRunning and not forceRedo:
            print('  Still running, skip')
            continue
            
        # Check if we need to run post process
        # Either a) we haven't done this yet, or
        # b) files have changed this we last did post processing
        outfile = 'cpp_diagnostics.out'
        cpp_diags_file = os.path.join(fullFolderLoc, outfile)
        if not os.path.exists(cpp_diags_file) or forceRedo:
            doPostProcess = True
        else:
            doPostProcess = False
        
        if doPostProcess:
            print('cpp_diagnostics.out does not exist')
        else:
            timeSincePostProcess = timeSinceFileUpdated(cpp_diags_file)
            if (timeSincePostProcess-1) > minutesSinceFolderUpdated:
                doPostProcess = True
                print('Need to redo cpp_diagnostics. Folder last updated: ' + str(minutesSinceFolderUpdated) + ' mins ago, postProcess last done ' + str(timeSincePostProcess) + ' mins ago')
                
        #doPostProcess = False
        if doPostProcess:

            extraParams = dict()
            
            if katzUnits:
                # add some things to extra params in order to convert to Wells units
                ##################################
                # To convert Wells to Katz
                ##################################
                extraParams['parameters.referenceTemperature'] = -3.5
                extraParams['parameters.referenceSalinity'] = 35
                    
                # Katz ref:
                extraParams['parameters.prevReferenceTemperature']= -23
                extraParams['parameters.prevReferenceSalinity'] = 230
                    
                # Convert boundary conditions:
                extraParams['parameters.bottomEnthalpy'] = float(inputs['parameters.bottomEnthalpy']) - 1.0
                extraParams['parameters.topEnthalpy'] =  float(inputs['parameters.topEnthalpy']) - 1.0
                extraParams['parameters.topBulkConc'] = float(inputs['parameters.topBulkConc']) + 1.0
                extraParams['parameters.bottomBulkConc'] = float(inputs['parameters.bottomBulkConc']) + 1.0
                    
                # Convert non dim vals:
                extraParams['parameters.compositionRatio'] = float(inputs['parameters.compositionRatio']) - 1.0

                
            # MAKE THIS RIGHT!
            ppexecDir = join(MushyLayerDir, 'postProcess')
            runPostProcessFolder(fullFolderLoc, extraParams, ppexecDir, outfile)
       
        
    print('Finished analysing all folders')
    
    
if __name__ == "__main__":
    print('Running script')
    main(sys.argv[1:])
     
