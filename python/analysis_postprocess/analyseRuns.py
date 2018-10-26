from MushyLayerRun import MushyLayerRun
import os
#import numpy
import math
from subprocess import Popen
from mushyLayerRunUtils import constructRunName, readInputs, writeParamsFile
from ChkFile import ChkFile
import sys, getopt
from os import listdir
from os.path import isfile, join
import shutil
import time
import re
from runPostProcessFolder import runPostProcessFolder

import getopt

def main(argv):
    folder = 'optimalStates-highRes-new'
    
    try:
       opts, args = getopt.getopt(argv,"f:",["f="])
    except getopt.GetoptError:
       print 'analyseRuns.py -f <base folder>'
       sys.exit(2)
    for opt, arg in opts:
        if opt in ("-f"):
            folder = str(arg)
            
    # Local
    #baseDir = '/home/parkinsonjl/mnt/sharedStorage/AW002_PARKINSON_MUSH/
    
    # Remote
    baseDir = os.path.join('/network/group/aopp/oceans/AW002_PARKINSON_MUSH/', folder)

    currentDir = os.getcwd()
    MushyLayerDir = currentDir.replace('python', '')

    katzUnits = True
    doRenameFinishedEarly = False

    # 4 processors should speed things up a little
    num_proc = 1

    foldersToProcess = []
    emptyFolders = []
    allFolders = []
    # Get list of folders that end in -0

    for f in next(os.walk(baseDir))[1]:
        #print(f[-2:])
        allFolders.append(f)
        if f[-2:] == '-0':
            foldersToProcess.append(f)
        else:
            emptyFolders.append(f)

    # Temporarily just process a few folders (for testing)
    #foldersToProcess = foldersToProcess[:5]
        
    # Remove empty folders    
    for emptyFolder in allFolders:
        fullFolderLoc = os.path.join(baseDir, emptyFolder)
        
        #filesInFolder = [f for f in listdir(fullFolderLoc) if (isfile(join(fullFolderLoc, f)) and not  f[:12] == '.fuse_hidden')]
        filesInFolder = [f for f in listdir(fullFolderLoc) if isfile(join(fullFolderLoc, f))]
        fuseFiles = [f for f in listdir(fullFolderLoc) if (isfile(join(fullFolderLoc, f)) and f[:12] == '.fuse_hidden')]
        
            
        
        
        if len(filesInFolder) == len(fuseFiles):
            print('Folder: ' + emptyFolder)
            #print('    Num files: ' + str(len(filesInFolder)))
            print('    ' + str(filesInFolder))
            
        #    delFile = raw_input('Remove folder? (y/n) ')
        #    if delFile == 'y':
        #        cmd  = 'rm -rf ' + fullFolderLoc
        #        #print(cmd)
        #        shutil.rmtree(fullFolderLoc, ignore_errors=True, onerror=None)
        #        print('    deleted')
            #else:
            #    print('    not deleted')
                
            # Just delete by default
            shutil.rmtree(fullFolderLoc, ignore_errors=True, onerror=None)
            print('    deleted')
            
        #else:
        #    print('    Folder is not empty')
            
         
    for folder in foldersToProcess:
        fullFolderLoc = os.path.join(baseDir, folder)
            
        inputsPath = os.path.join(fullFolderLoc, 'inputs')
        if not os.path.exists(inputsPath):
            print('No inputs file found')
            continue
        
        inputs = readInputs(inputsPath)
        
        
        print(folder)
        
        # Get the final file
        files = [f for f in listdir(fullFolderLoc) if (isfile(join(fullFolderLoc, f)) and f[:3] == 'chk')]
        allFiles = [f for f in listdir(fullFolderLoc) if isfile(join(fullFolderLoc, f))]
        allFiles = sorted(allFiles)
        files = sorted(files)
        restart_file = ''
        if len(files) > 0:
            restart_file = files[-1]
        #else:
        #    print('no restart files!')
        #    print(files)
        #    continue


        if not os.path.exists(join(fullFolderLoc, 'pout.0')):
            print('No pout file - finished early')
            newName = folder[:-1] + 'finishedEarly'
            
            os.rename(fullFolderLoc, join(baseDir, newName))
            continue

        # Check we haven't finished early
        # First find how recently pout.0 was edited
        st=os.stat(join(fullFolderLoc, 'pout.0'))
        age=float(time.time()-st.st_mtime) # Should be in seconds
        
        
        steadyState = os.path.isfile(join(fullFolderLoc, 'time.table.0'))
        
        if not steadyState:
            # Try slower test
            with open(join(fullFolderLoc, 'pout.0'), 'r') as f:
                poutContents = f.read().replace('\n', '')
                if poutContents.find('all amr levels said steady state was reached') > -1:
                    steadyState = True

        # 5 minutes
        if steadyState:
            # Check the final dt wasn't tiny
            small_dt = False
            with open(join(fullFolderLoc, 'pout.0'), 'r') as f:
                poutContents = f.read()
                #poutContents = 'euikwfuwibg dt feoghow'
                #print(poutContents)
                matches = re.findall(r'dt = ([^\s]+)  wallclocktime', poutContents)
                if matches:
                    final_dt = float(matches[-1])
                    #for m in matches:
                        #print(m)
                    if final_dt < 1e-10:
                        small_dt = True
                        break
                            
                else:
                    print('    no matches: ')
                    print(matches)
                    
            if small_dt:
                print('    dt < 1e-10, deleting folder')
                #newName = folder[:-1] + 'small_dt'
                #os.rename(fullFolderLoc, join(baseDir, newName))
                #shutil.rmtree(fullFolderLoc, ignore_errors=True, onerror=None)
            else:
                print('    Reached steady state')
                
                
                
                
                # At this point, execute some analysis code if we haven't already
                outfile = 'cpp_diagnostics.out'
                if not os.path.exists(os.path.join(fullFolderLoc, outfile)):
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
                    
                
                newName = folder[:-1] + 'steady'
                print('    Renaming to + ' + newName)
                try:
                    os.rename(fullFolderLoc, join(baseDir, newName))
                except OSError as (errno, strerror):
                    print "OS error({0}): {1}".format(errno, strerror)
            
        elif age < 300:
            pass
            print('    Still running')
        elif doRenameFinishedEarly:
            
            for aFile in allFiles:
                print('        ' + aFile)
                
            print('    stopped before steady state, pout.0 last edited ' + str(age) + ' seconds ago')
                
            newName = folder[:-1] + 'finishedEarly'
            # stop renaming
            #renameFolder = raw_input('Rename folder to ' + newName + '? (y/n) ')
            #if renameFolder == 'y':
            #    os.rename(fullFolderLoc, join(baseDir, newName))
            #try:
            #    os.rename(fullFolderLoc, join(baseDir, newName))
            #except OSError as (errno, strerror):
            #    print "OS error({0}): {1}".format(errno, strerror)
            
            
            #print('    ' + str(allFiles))
        
    print('Finished analysing all folders')


if __name__ == "__main__":
    print('Not running analyseRuns - be careful and check this script does not delete anything first.')
    #main(sys.argv[1:])
     
