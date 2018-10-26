    # Script to test convergence of uniform and AMR solutions

import os
import getpass
import sys
parDir = os.path.abspath(os.pardir)
pythonDir = os.path.join(parDir, 'python')
sys.path.append(pythonDir)

from MushyLayerRun import MushyLayerRun
import math
from subprocess import Popen
from mushyLayerRunUtils import constructRunName, readInputs, writeParamsFile, getFinalPlotFile
#from AMRConvergenceTest import AMRConvergenceTest
from SlurmTask import SlurmTask
import re

import getopt

def getRes(folder):
    m = re.findall('-(\d+)-', folder)
    if m:
        resolution = int(m[0])
        #print('Resolution of folder ' + folder + ' determined as ' + str(resolution))
        return resolution

def makeCompareScript(baseDir, coarseFolder, fineFolder, execFile, runType):
    #TODO: write this
    coarseRes = getRes(coarseFolder)
    fineRes = getRes(fineFolder)

    print('Run type: ' + runType)


    outputFolderName = runType + str(coarseRes) + '-' + str(fineRes)
    outputFolder = os.path.join(baseDir, 'ChomboCompare', outputFolderName )

    print('Output folder name: ' + outputFolderName)

    # Don't redo analysis. 
    # Create folder if we haven't done analysis yet
    if os.path.exists(outputFolder):
        return
    else:
        os.makedirs(outputFolder)   

    name = 'compare' + outputFolderName

    #runFileName = name + '.sh'
    runFileName = 'run.sh'
    inputsFileName = 'inputs'
    #inputsFileName = name + '.inputs'

    inputsFile = os.path.join(outputFolder, inputsFileName)
    #runFile = os.path.join(outputFolder, runFileName)
    errFile = os.path.join(outputFolder, name + '.hdf5')
    exactFile = getFinalPlotFile(fineFolder)
    coarseFile = getFinalPlotFile(coarseFolder)

    if not exactFile:
        print('Could not find exact file in directory ' + fineFolder)
        return

    print('exactFile: ' + exactFile)
    print('coarseFile: ' + coarseFile)

    # Construct inputs
    params = {}
    params['compare.sameSize'] = 0

    #name of file for "exact" solution (single-level fine-resolution solution)
    params['compare.exactRoot'] = os.path.join(fineFolder, exactFile)

    #name of file containing "computed" solution (which may be AMR )
    params['compare.computedRoot'] = os.path.join(coarseFolder, coarseFile)

    #name of file into which to write error plotfile
    params['compare.errorRoot'] = errFile

    #these are only important if you want to do a time-series of plotfiles
    params['compare.isTimeDep'] = 0
    #compare.numCrseStart = 0
    #compare.crseStep = 1
    #compare.mult = 16

    #dump out a plotfile of the error?
    params['compare.doPlots'] = 1

    #use higher-order averaging? (1 = yes)
    # generally, want to use HO averaging for point values,
    # non-HO-averaging when computing cell-average (conservative) 
    # solutions
    params['compare.HOaverage'] = 0 # try turning this off

    #which variables to compute error for
    #params['compare.error_var'] = '' # blank string - do them all

    #these variables are not compared against a finer plotfile
    #(it's assumed that the "exact" solution is 0)
    #params['compare.no_average_var'] = 'divergence'

    

    writeParamsFile(inputsFile, params)

    s = SlurmTask(outputFolder, name, execFile)
    s.writeSlurmFile(runFileName, inputsFileName)
    s.runTask(runFileName)

    


def main(argv):
    
    # Defaults:
    folder = 'MushyConvectionLiquid'
    highestResolution = 1024
    
    try:
       opts, args = getopt.getopt(argv,"f:r:")
    except getopt.GetoptError:
       print 'computeError.py -f <folder> -r <highest resolution>'
       sys.exit(2)
    for opt, arg in opts:
        if opt in ("-f"):
            folder = str(arg)
        if opt in ("-r"):
            highestResolution = int(arg)
        
    print 'Analysing folder: ', folder
    print 'Highest resolution: ', str(highestResolution)

    #################################
    # These shouldn't change
    #################################
    #cwd = os.getcwd()
    #mushyLayerBaseDir = cwd.replace('/python', '')
    username = getpass.getuser()
    mushyLayerBaseDir = os.path.join('/home',username, 'convection-in-sea-ice','MushyLayer')

    os.environ["CH_TIMER"] = "1"
    dropbox_data_folder = "/home/parkinsonjl/Dropbox/DPhil/Data/"
    subcycled = True # subcycled code

    data_dir = '/network/group/aopp/oceans/AW002_PARKINSON_MUSH/' # for writing to shared storage
    baseDir = os.path.join(data_dir, 'AMRConvergenceTest', folder)

    chomboDir = mushyLayerBaseDir.replace('convection-in-sea-ice/MushyLayer', 'chombofork')
    execFile = os.path.join(chomboDir, 'lib/util/ChomboCompare/compare2d.Linux.64.mpiCC.gfortran.OPT.MPI.ex')

    # Get all folders
    folders = [f for f in os.listdir(baseDir) if os.path.isdir(os.path.join(baseDir, f))]
    
    # Find fine res uniform run
    fineFolder = '' 

    for f in folders:
        m = re.findall('Uniform(.*)-' +str(highestResolution)+ '-*', f)
        if str(highestResolution) in f:

            if m:
                fineFolder = os.path.join(baseDir, f)
                print('Fine uniform run is ' + f)
           
    if not fineFolder:
        print('Unable to find folder with resolution ' + str(highestResolution) + ' in output folder: ' + folder)
        sys.exit()



    # Create run scripts and execute for each folder

    # For comparison to fine res:
    for f in folders:
        if f[-2:] == '-0':
            coarseFolder = os.path.join(baseDir, f)

            timeTableFile = os.path.join(coarseFolder, 'time.table.0')
            if os.path.exists(timeTableFile):

                #prefix = re.sub(folder+'.*', '', f)
                prefix = re.sub('\d+--0$', '', f)


                makeCompareScript(baseDir, coarseFolder, fineFolder, execFile, prefix)


    # For Richardson error:
    prefix = 'Richardson'
    for f in folders:
        if f[-2:] == '-0' and f[:7] == 'Uniform':
            coarseFolder = os.path.join(baseDir, f)

            timeTableFile = os.path.join(coarseFolder, 'time.table.0')
            if not os.path.exists(timeTableFile):
                continue

            coarseRes = getRes(coarseFolder)
            fineRes = int(coarseRes)*2
            fineFolder = coarseFolder.replace(str(coarseRes), str(fineRes))

            if os.path.exists(fineFolder):

                makeCompareScript(baseDir, coarseFolder, fineFolder, execFile, prefix)



    


if __name__ == "__main__":
    #print('Running script')
    main(sys.argv[1:])
                 
