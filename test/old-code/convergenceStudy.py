import numpy as np
import matplotlib.pyplot as plt
from MachineSpecific import MachineSpecific
from os.path import join, isfile
import re
import ModelOutput.modelOutput

#plot two convergence studies on the same graph. 
#Assumes all convergence studies have the same:
#compareFrom, compareTo
def plotError(cStudyArr, var, errType, plotComparison):
    plt.figure()
    
    log_errArrays = []
    compareFromPts_log_arrs = []
    plot_titles = []
    
    
    for cStudy in cStudyArr:
        [log_errArray, firstOrderComp, secondOrderComp] = cStudy.plotVars(var, errType)
         
        log_errArrays.append(log_errArray)
        compareFromPts_log_arrs.append(cStudy.compareFromPts_log_arr)
        
        plt.plot(cStudy.compareFromPts_log_arr, log_errArray, label=cStudy.plotTitle, marker='x')
        
        parts = str.split(cStudy.data_dir,'/')
        plot_titles.append(parts[-1])
        
    if plotComparison:
        plt.plot(compareFromPts_log_arrs[0], firstOrderComp, label="1st order", marker='x')
        plt.plot(compareFromPts_log_arrs[0], secondOrderComp, label="2nd order", marker='x')
    
    plt.legend(loc = 3) 
    plt.title(errType + " Error convergence of " + cStudyArr[0].compareFrom + " solution to " + cStudyArr[1].compareTo + " solution ("+convergenceStudy.VAR_NAMES[var]+")")
    plt.xlabel('log10(grid points)')
    plt.ylabel('log10(' + errType + ' error)')
    plt.grid()
    
    filename = 'comparison-'+'-'.join(plot_titles)+".png" 
    plt.savefig(filename)
         
         
    plt.draw()
    
def plotErrorTime(cStudyArr, var, errType):
    plt.figure()
   
    plot_titles = []
    
    for cStudy in cStudyArr:
        [log_errArray, firstOrderComp, secondOrderComp] = cStudy.plotVars(var, errType)

        plt.plot(np.log10(cStudy.getTimes()), log_errArray, label=cStudy.plotTitle, marker='x')
        
        parts = str.split(cStudy.data_dir,'/')
        plot_titles.append(parts[-1])
        
    
    plt.legend(loc = 3) 
    plt.title(errType + " Error (" + cStudyArr[0].compareFrom + " solution vs " + cStudyArr[1].compareTo + " solution) and run time ("+convergenceStudy.VAR_NAMES[var]+")")
    plt.xlabel('log10(runtime)')
    plt.ylabel('log10(' + errType + ' error)')
    plt.grid()
    
    filename = 'comparison-'+'-'.join(plot_titles)+"-time.png" 
    plt.savefig(filename)
         
         
    plt.draw()
    
def plotTime(cStudyArr):
    plt.figure()
   
    
    
    for cStudy in cStudyArr:
        plt.plot(cStudy.compareFromPts_log_arr, np.log10(cStudy.getTimes()),  label=cStudy.plotTitle, marker='x')
        
       
    
    plt.legend(loc = 4) 
    plt.title("runtime vs resolution")
    plt.xlabel('log10(# grid points)')
    plt.ylabel('log10(runtime)')
    plt.grid()
    
#     filename = 'comparison-'+'-'.join(plot_titles)+"-time.png" 
#     plt.savefig(filename)
         
         
    plt.draw()
    
     

class convergenceStudy:
    
    
    
    #Constants - don't change! (unless I change some of the mushyLayer code)
    ANALYTIC = ModelOutput.ANALYTIC
    CALCULATED = ModelOutput.CALCULATED
    L1ERR = ModelOutput.L1ERR
    L2ERR = ModelOutput.L2ERR
    MAX = ModelOutput.MAX
    SUM = ModelOutput.SUM
    CONCENTRATION = ModelOutput.composition
    THETA = ModelOutput.theta
    POROSITY = ModelOutput.porosity
    LIQUID_CONC = ModelOutput.compositionLiquid
    
    VAR_NAMES = ModelOutput.VAR_NAMES
    
    PREAMBLE_LINES = ModelOutput.PREAMBLE_LINES
    NUM_VARS = ModelOutput.NUM_VARS
   
    

    
    def __init__(self, compareFromPts_arr, compareToPts_arr, compareFrom, compareTo, parallel, data_dir):
        self.machine_specific = MachineSpecific()
        
        
        self.compareFromPts_arr = compareFromPts_arr
        self.compareToPts_arr = compareToPts_arr
        self.compareFrom = compareFrom
        self.compareTo = compareTo
        self.parallel = parallel
        self.data_dir = data_dir
        
         
        parts = str.split(self.data_dir,'/')
        if parts[-1]:
            self.plotTitle = parts[-1]
        else:
            self.plotTitle = parts[-2]
        
        #this allows use to just specify one file to compare against if we're being lazy
        if (len(self.compareToPts_arr) != len(self.compareFromPts_arr)):
            if ( len(self.compareToPts_arr) == 1):
                self.compareToPts_arr = [self.compareToPts_arr[0]] * len(self.compareFromPts_arr)
            else:
                print("Error - the number of files to compare against is different than the number to compare from.")
                return
        
        
            
        self.doComparison()
        
        
    def doComparison(self):
        self.compareFromPts_log_arr = []
        self.L1 = []
        self.L2 = []
        self.Max = []
        self.Sum = []
        self.times = []
        self.chomboCompareOut = []
        
        pts_i=0

        while pts_i < len(self.compareFromPts_arr):
            compareFromPts = self.compareFromPts_arr[pts_i]
            compareToPts = self.compareToPts_arr[pts_i]
            
            chCompare = ModelOutput(compareFromPts, compareToPts,
                                    self.compareFrom, self.compareTo,
                                    self.machine_specific.get_home_dir() + self.data_dir,
                                    self.parallel)
                
             #   chCompare.readOutput(filename, exactFile, timingFilename, poutFilename, lines)  
                
                    
            #Store this analysis and move on to the next comparison  
            
            self.chomboCompareOut.append(chCompare)
            
    
            self.compareFromPts_log_arr.append(np.log10(compareFromPts))
            pts_i = pts_i + 1
            
    

        
        
    # Get all the variables needed for plotting
    def plotVars(self, var, errType):        
        allErrors = []
        
        for chComp in self.chomboCompareOut:
            print(chComp.errors)
            if var in chComp.errors:
                var_errors = chComp.errors[var]
                print(var_errors)
                allErrors.append(var_errors[errType])
            else:
                allErrors.append(float('nan'))
            

            
       
        secondOrderComp = []
        firstOrderComp = []        
        log_errArray = []
        
        i = 0
        for err in allErrors:
            #print(err)
             
            log_errArray.append(np.log10(err))
            firstOrderComp.append(log_errArray[0] - (self.compareFromPts_log_arr[i] - self.compareFromPts_log_arr[0]))
            secondOrderComp.append(log_errArray[0] - 2*(self.compareFromPts_log_arr[i] - self.compareFromPts_log_arr[0]))
            i=i+1
            
            
        return [log_errArray, firstOrderComp, secondOrderComp]
        
        
        
    #Go through the array of this error type for all vars, and extract the var we're interested in
    def plotErrorResolution(self, var, errType):
        plt.figure()
        [log_errArray, firstOrderComp, secondOrderComp] = self.plotVars(var, errType)
    
        plt.plot(self.compareFromPts_log_arr, log_errArray, label="convergence", marker='x')
        plt.plot(self.compareFromPts_log_arr, firstOrderComp, label="1st order", marker='x')
        plt.plot(self.compareFromPts_log_arr, secondOrderComp, label="2nd order", marker='x')
        plt.legend(loc = 3)
        plt.title(errType + " Error convergence of " + self.compareFrom + " solution to " + self.compareTo + " solution ("+self.VAR_NAMES[var]+")")
        plt.xlabel('log10(grid points)')
        plt.ylabel('log10(' + errType + ' error)')
        
       
        save_to_location = self.machine_specific.get_home_dir() + self.data_dir
        save_to_location += errType + "-" + self.VAR_NAMES[var] + "-" + "convergence.png" 
        plt.savefig(save_to_location)
        
        plt.draw()
        
    def getTimes(self):
        
        times = [ch.time for ch in self.chomboCompareOut]

        return times
        
    #Plot error vs runtime
    def plotErrorTime(self, var, errType):
        plt.figure()
        [log_errArray, firstOrderComp, secondOrderComp] = self.plotVars(var, errType)
        
        #print([ch.time for ch in self.chomboCompareOut])
        
        logTime = [np.log10(ch.time) for ch in self.chomboCompareOut]
    
        plt.plot(logTime, log_errArray, marker='x')
        
        #plt.legend(loc = 3)
        plt.title(errType + " Error (" + self.compareFrom + " solution vs " + self.compareTo + " solution) vs runtime ("+self.VAR_NAMES[var]+")")
        plt.xlabel('log10(runtime)')
        plt.ylabel('log10(' + errType + ' error)')
        
       
        save_to_location = self.machine_specific.get_home_dir() + self.data_dir
        save_to_location += errType + "-" + self.VAR_NAMES[var] + "-" + "time.png" 
        plt.savefig(save_to_location)
        plt.grid()
        
        plt.draw()
        
    def plotTimeResolution(self):
        plt.figure()
        
        logTime = [np.log10(ch.time) for ch in self.chomboCompareOut]
        plt.plot(self.compareFromPts_log_arr, logTime, marker='x')
        
        plt.title(self.data_dir + ", runtime vs grid_resolution")
        plt.xlabel('log10(# grid points)')
        plt.ylabel('log10(runtime)')
        
        plt.grid()
        plt.draw()
   
    
            