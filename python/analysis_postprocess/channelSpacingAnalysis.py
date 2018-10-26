from PltFile import PltFile
import numpy as np
#from matplotlib import mpl,pyplot
from os import listdir
from os.path import isfile, join
import re
import os
import csv
from MachineSpecific import MachineSpecific


#folders = ['/media/parkinsonjl/FREECOM HDD/channelSpacing-periodic/CR1.25RaC800Le200ChiCubedPermeabilitypts1024-0'

# Which Rayleigh numbers to recalculate
#base_folder = '/media/parkinsonjl/FREECOM HDD/channelSpacing-periodic/'
#base_folder = '/media/parkinsonjl/FREECOM HDD/channelSpacing-1proc-periodic/'

m = MachineSpecific()
if m.machine == m.LAPTOP:
    base_folder = '/home/parkinsonjl/mnt/sharedStorage/channelSpacing/'
else:
    base_folder = '/network/group/aopp/oceans/AW002_PARKINSON_MUSH/channelSpacing/'
    
outputFilename = 'channelSpacing.csv'

# Option 1: Just do all folders within base_folder
folders = []
folders = next(os.walk(base_folder))[1]
for i in range(0, len(folders)):
	folders[i] = os.path.join(base_folder, folders[i])
	
#base_folder + 'CR1.25RaC80Le200ChiCubedPermeabilitypts2048-0',
#folders = [base_folder + 'CR1.25RaC90Le200ChiCubedPermeabilitypts2048-0',
#base_folder + 'CR1.25RaC100Le200ChiCubedPermeabilitypts2048-0',
#base_folder + 'CR1.25RaC125Le200ChiCubedPermeabilitypts2048-0',
#base_folder + 'CR1.25RaC150Le200ChiCubedPermeabilitypts2048-0',
#base_folder + 'CR1.25RaC200Le200ChiCubedPermeabilitypts2048-0',
#base_folder + 'CR1.25RaC250Le200ChiCubedPermeabilitypts2048-0',
#base_folder + 'CR1.25RaC100Le200ChiCubedPermeabilitypts1024-0',
#base_folder + 'CR1.25RaC150Le200ChiCubedPermeabilitypts1024-0',
#base_folder + 'CR1.25RaC200Le200ChiCubedPermeabilitypts1024-0',
#base_folder + 'CR1.25RaC250Le200ChiCubedPermeabilitypts1024-0',
#base_folder + 'CR1.25RaC300Le200ChiCubedPermeabilitypts1024-0'
#]

filesPerFolder = 1000

# Option 2:
#folders = []
#Ras = [200, 300]
#for Ra in Ras:
#    folders.append(base_folder + 'CR1.25RaC' + str(Ra) + 'Le200ChiCubedPermeabilitypts1024-0')
    
for folder in folders:
    print('==========================================================')
    print('Processing ' + folder)
    
    diagFile = join(folder, outputFilename)
    alreadyProcessedFiles = []
    
    oldFilename = join(folder, 'channelSpacing.out')
    if os.path.isfile(oldFilename):
        os.rename(oldFilename, diagFile)
    
    if os.path.isfile(diagFile):
        
        print('Diag file already exists')
        
        # Skip these for now
        #continue
        
        modernFormat = True
       
        with open(diagFile) as csvfile:
            #print('Opened diag file as csv file')
            
            reader = csv.DictReader(csvfile)
              
            
            for row in reader:
                print(row)
                if 'file' in list(row.keys()):      
                    alreadyProcessedFiles.append(row['file'])
                    
                else:
                    modernFormat = False
                
        # Check it's the correct format CSV file
        if modernFormat:        
            F = open(diagFile, 'a') # open to append
        else:
            F = open(diagFile, 'w') # open to overwrite as wrong file format before
            headerStr = 't,spacing,file'
            print(headerStr)
            F.write(headerStr + '\n')
                
                
        print('Already processed: ' + str(alreadyProcessedFiles))
    else:
        F = open(diagFile, 'w') # open to write to
        headerStr = 't,spacing,file'
        print(headerStr)
        F.write(headerStr + '\n')
    
    #folder = '/media/parkinsonjl/FREECOM HDD/channelSpacing-periodic/CR1.25RaC400Le200ChiCubedPermeabilitypts1024-0'
    
    #time_interval = float(100)
    
    # Find all plot files in folder, load them, calculate chimney spacing diagnostics, and write out file with results
    
    # Get all files in folder
    files = [f for f in listdir(folder) if isfile(join(folder, f))]
    files = sorted(files)
    
    #print(files)
    
    # Filter out any that aren't plt files
    pltFiles = []
    nFiles = 0
    for filename in files:
        # Skip files we've already processed
        if filename in alreadyProcessedFiles:
            continue
            
        if nFiles >= filesPerFolder:
            continue
            
        try:
        
            m = re.search('^(?!chk).*\-(\d+)\.2d\.hdf5', filename)
            if m:
                timestep = float(m.group(1))
                #if timestep % time_interval == 0:
                
                full_dir = join(folder, filename)
                p = PltFile(full_dir)
                #print('Loaded ' + file)
                
                [numChannels, relPos, channelSpacing] = p.channelProperties(False)
        
                lev0Dx = p.levels[0][p.DX]
                #print('dx = ' + str(lev0Dx))
                domainSize = p.domainSize
                #print('Domain extent: ' + str(domainSize))
                if numChannels == 0:
                    avChanSpacing = float('NaN')
                else:
                    #avChanSpacing = np.mean(channelSpacing)
                    avChanSpacing = (p.domainSize[2] - p.domainSize[0])/float(numChannels)
                
                #print('Num channels: ' + str(numChannels) + ', relative positions (0 = left, 1 = right): ' + str(relPos))
                #print('Average channel spacing: ' + str(avChanSpacing))     
                #p.plotField('Porosity')
                
                #channelSpacings.append(avChanSpacing)
                #times.append(p.time)
            
                outLine = '%1.10f,%-1.5f,%s' % (p.time, avChanSpacing, filename)
                print(outLine)
                F.write(outLine + '\n')
                nFiles = nFiles + 1
                
        except Exception as e:
            logging.error(traceback.format_exc())
        
            
            
    F.close()
        
