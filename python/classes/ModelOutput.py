import re
import numpy as np
from os.path import isfile
from subprocess import Popen
from glob import glob
from MachineSpecific import MachineSpecific
import matplotlib.pyplot as plt


def plotNumIters(outputArray):
    plt.figure()

    for out in outputArray:
        timestep_nums = [tstep[out.TIMESTEP] for tstep in out.timesteps]
        timestep_num_iters = [tstep[out.NUM_ITERATIONS] for tstep in out.timesteps]

        total_iters = np.sum(timestep_num_iters)
        lab = out.computedFilenameTitle() + " (" + str(total_iters) + ")"


        plt.plot(np.log10(timestep_nums), np.log10(timestep_num_iters), label=lab)

    plt.xlabel('log10(Timestep)')
    plt.ylabel('log10(Num iterations in timestep)')

    plt.title('num iterations')

    plt.legend(loc=3)
    plt.grid()

    plt.draw()

def chomboCompareProgram(parallel):
    #Machine specific vars
    machSpec = MachineSpecific()

    parallelProgram = "/compare2d.Linux.64.mpiCC.gfortran.DEBUG.MPI.ex"
    serialProgram = "/compare2d.Linux.64.g++.gfortran.DEBUG.ex"
    ChomboCompareDir = machSpec.get_home_dir() + "/Chombo/lib/util/ChomboCompare"

    if (parallel):
        program = ChomboCompareDir + parallelProgram
    else:
        program = ChomboCompareDir + serialProgram

    return program

class modelOutput:
    enthalpy = 0
    enthalpySolid = 1
    enthalpyEutectic = 2
    enthalpyLiquid = 3
    composition = 4
    theta = 5
    thetaLiquidus = 6 
    thetaSolidus = 7
    porosity = 8
    porosityEutectic = 9 
    compositionLiquid = 10
    compositionSolid = 11
    thetaForcing = 12
    liquidCompositionGrad = 13 
    solidFraction = 14
    steadyStateImbalance = 15 
    Tanalytic = 16
    ThetaLAnalytic = 17 
    solidFractionTrue = 18 
    thetaTrue = 19
    ThetaAnalytic = 20 
    enthalpyAnalytic = 21 
    enthalpyAdvection = 22
    thetaLaplacian = 23
    streamFunction = 24
    resid = 25
    concSource = 26
    ThetaDiffusion = 27 
    ThetaDiffusionN = 28
    thetaBcoef = 29
    ThetaBCoef = 30
    permeability = 31
    pressure = 32
    divU = 33
    ThetaFrameAdvection = 34 
    ThetaSSource = 35
    ThetaPorositySource = 36
    
    PREAMBLE_LINES = 5 
    NUM_VARS = 36
    
    L1ERR = "L1"
    L2ERR = "L2"
    MAX = "Max"
    SUM = "Sum"
    
    VAR_NAMES = {composition: 'bulk concentration',
                 theta: 'theta',
                 porosity: 'porosity',
                 compositionLiquid: 'liquid concentration'}
    
    TIMESTEP = "timestep"
    TIME = "time"
    DT = "dt"
    ITERATIONS = "iterations"
    NUM_ITERATIONS = "num iterations"
    
    ITERATION_NUM = "iteration_num"
    HEAT_RESIDUAL = "heat residual"
    CONC_RESIDUAL = "conc residual"
    
    ANALYTIC = "analytic"
    CALCULATED = "steady"
    

    
            
    def __init__(self, compareFromPts, compareToPts, compareFrom, compareTo, dataDir, parallel):
        self.data_dir = dataDir
        self.compareFrom = compareFrom
        self.compareTo = compareTo
        self.errors = {}
        self.machine_specific = MachineSpecific()
        self.parallel = parallel
        
        self.test_dir = self.machine_specific.get_home_dir() + '/convection-in-sea-ice/test'
        self.program = chomboCompareProgram(parallel)
        
        #First check if we've already done this comparison
        self.saveFile = self.data_dir + "err" + str(compareFromPts) + self.compareFrom + "-" + str(compareToPts) + self.compareTo + ".dat"
        
        #Create some initial defaults
        self.time = float('nan')
        
        if isfile(self.saveFile ):
            self.read()
        
        else:
        
            #If we haven't already done this comparison, do it now!
            #First get the files we need to read in
            
            self.computedFilename = self.data_dir + "bm1-" + str(compareFromPts) + "pts-" + self.compareFrom + ".2d.hdf5"
            
            if not isfile(self.computedFilename):
                print("File to compare , "+self.computedFilename+", doesn't exist. Skipping this comparison.")
                return
    
    
            if (self.compareTo == self.ANALYTIC):
                self.exactFilename = self.data_dir + "bm1-"+str(compareToPts)+"pts-analytic.2d.hdf5"
            elif(self.compareTo == self.CALCULATED):
                self.exactFilename = self.data_dir + "bm1-"+str(compareToPts)+"pts-steady.2d.hdf5"
    
            print("File to compare against: " + self.exactFilename)
    
            if (self.computedFilename == self.exactFilename):
                print("Skipping " + self.computedFilename)
                return
            
            if not isfile(self.exactFilename):
                print("File to compare against, "+self.exactFilename+", doesn't exist. Skipping this comparison.")
                return
    
            
        
            print("Comparing " + self.computedFilename)
    
            errFilename = self.data_dir + "err" + str(compareFromPts) + self.compareFrom + "-" + str(compareToPts) + self.compareTo + ".hdf5";
        

            process = Popen(self.program + " " + self.test_dir + "/bm1inputs.compare" + 
        " compare.computedRoot=" + self.computedFilename + 
        " compare.exactRoot=" + self.exactFilename + 
        " compare.errorRoot=" + errFilename, shell=True)
            (outputOneLevel, err) = process.communicate()
            
            exit_code = process.wait()
            print("exit code = " + str(exit_code))
            if exit_code == -6:
                print("Couldn't find one or more of the input files. Aborting.")
                return
            elif exit_code != 0:
                print("ChomboCompare didn't run properly. Aborting.")
                return
    
            #in parallel, the outputOneLevel is actually put in a file - pout.0
            if (self.parallel):
                lines = open('pout.0','r')
            else:
                lines = outputOneLevel.split("\n")
                
            #Let's also get how long this simulation took to run
            #Find the the most recent folder for this resolution, and get the
            #details from the time.table.0 file
            
            #paths = glob('*/')
            pathSearch = self.data_dir + str(compareFromPts) + "-*/"
            paths = glob(pathSearch)
            
            self.timingFilename = ""
            self.poutFilename = ""
            
            if len(paths) > 0:
             
                pathDict = {}
                pattern = self.data_dir + str(compareFromPts) + '-(\d+)/'
                
                for path in paths:
                    m = re.match(pattern, path)
                    if m:
                        pathDict[int(m.group(1))] = path
                 
                keys = pathDict.keys()
                keys.sort()
                 
                mostRecentDir = pathDict[keys[-1]]
                 
                self.timingFilename = mostRecentDir + 'time.table.0'
                self.poutFilename = mostRecentDir + 'pout.0'
                
            self.parseOutputFiles(lines)
    
    def parseOutputFiles(self, outputLines):
        #Read in outputOneLevel from lines to get errors
        i = 0
        var_i = 0
        for line in outputLines:
            #print (str(i) + ": " + line)
            #print (line)

            # Don't parse lines that don't contain errors
            if (i < self.PREAMBLE_LINES or i > (self.PREAMBLE_LINES + self.NUM_VARS)):
                i = i+1
                continue

            parts = line.split(":")
            err = parts[1].split(",")
            
            
            err_array = [float(x) for x in err[0:4]]
            err_dict = {self.L1ERR: err_array[0],
                        self.L2ERR: err_array[1],
                        self.MAX: err_array[2],
                        self.SUM: err_array[3],}
            
            self.errors[var_i] = err_dict
            
            i=i+1
            var_i = var_i+1
            
        # get the time taken to run
        if self.timingFilename and isfile(self.timingFilename):
            timingFile = open(self.timingFilename, 'r')
            timePattern = '\[0\]main\s+(\d+\.\d+)'
                    
            for line in timingFile:
                m = re.match(timePattern, line)
                if m:
                    self.time = float(m.group(1))
                    break
        else:
            self.time = -1
            
            
        # Get information we want from the pout file
        self.timesteps = []
        regex_float = "([-+]?\d*\.?\d*(?:e-\d+)?)"
        regex_timestep = "Timestep (\d+) Advancing solution from time "+regex_float+" with dt = "+regex_float
        regex_iteration = "\s+amrMushyLayer::timestep\(\) Iteration (\d+). Heat Residual = "+regex_float+" \(converged = ([01]), rate = "+regex_float+"\), Concentration residual = "+regex_float+" \(converged = ([01]), rate = "+regex_float+"\)"
        regex_converged = "amrMushyLayer::timestep ,\s+end time =\s+"+regex_float+"\( "+regex_float+" secs\) , dt = \s+"+regex_float+"\( "+regex_float+" secs\) "
        
        if isfile(self.poutFilename):
            thisTimestep = {}
            pout = open(self.poutFilename, 'r')
            for line in pout:
                m_timestep = re.match(regex_timestep, line)
                m_iteration = re.match(regex_iteration, line)
                m_converged = re.match(regex_converged, line)
                if m_timestep:
                                    
                    #Start another timestep
                    thisTimestep = {}
                    timestep = m_timestep.group(1)
                    time = m_timestep.group(2)
                    dt = m_timestep.group(3)
                    thisTimestep[self.TIMESTEP] = int(timestep)
                    thisTimestep[self.TIME] = float(time)
                    thisTimestep[self.DT] = float(dt)
                    
                    thisTimestep[self.ITERATIONS] = []
                    
                elif m_iteration:
                    iterationNum = m_iteration.group(1)
                    heatResidual = m_iteration.group(2)
                    concResidual = m_iteration.group(5)
                    
                    thisIteration = {self.ITERATION_NUM: int(iterationNum),
                                     self.HEAT_RESIDUAL: float(heatResidual),
                                     self.CONC_RESIDUAL: float(concResidual)}
                    
                    thisTimestep[self.ITERATIONS].append(thisIteration)
                    
                elif m_converged:
                    thisTimestep[self.NUM_ITERATIONS] = len(thisTimestep[self.ITERATIONS])
                    self.timesteps.append(thisTimestep)
                    
                
        #print(self.timesteps)
            
        #print(self.errors)
        self.write()
        
    def plotConvergenceWithSign(self, timesteps):
           
        plt.figure()
        
        for timestep in timesteps:
            if timestep >= len(self.timesteps):
                continue
            
            thisTimestep = self.timesteps[timestep]
            thisTimestepIterations = thisTimestep[self.ITERATIONS]
            
            iter_nums =   [iter[self.ITERATION_NUM] for iter in thisTimestepIterations]
            heat_residuals_sign = [iter[self.HEAT_RESIDUAL] for iter in thisTimestepIterations]
            #heat_residuals_log_sign = [np.sign(resid)*np.log10(abs(resid)) for resid in heat_residuals_sign]
         
         
            
            plt.plot(iter_nums, heat_residuals_sign, label=str(timestep))
            
            plt.xlabel('# iterations')
            plt.ylabel('enthalpy residual')
            
        plt.legend(loc = 3) 
        plt.title( str(self.computedFilename))
        plt.grid()
        plt.draw()
        
   
    def plotConvergence(self, timesteps, residual_var):
        plt.figure()
        
        for timestep in timesteps:
            if timestep >= len(self.timesteps):
                continue
            
            thisTimestep = self.timesteps[timestep]
            thisTimestepIterations = thisTimestep[self.ITERATIONS]
            
            iter_nums =   [iter[self.ITERATION_NUM] for iter in thisTimestepIterations]
            residuals = np.abs([iter[residual_var] for iter in thisTimestepIterations])
            
            
            min_vals = [1e-6 for resid in residuals]
            residuals = np.maximum(residuals, min_vals)
            print(residuals)
            
            log_heat_resid = np.log10(residuals)
            
            plt.plot(np.log10(iter_nums), log_heat_resid, label=str(timestep))
            
            plt.xlabel('log10(# iterations)')
            plt.ylabel('log10(' + residual_var +' residual)')
            
        plt.legend(loc = 3) 
        plt.title( str(self.computedFilenameTitle()))
        plt.grid()
        plt.draw()
        
    def plotNumIterations(self):
        timestep_nums = [tstep[self.TIMESTEP] for tstep in self.timesteps]
        timestep_num_iters = [tstep[self.NUM_ITERATIONS] for tstep in self.timesteps]
        
        total_iters = np.sum(timestep_num_iters)
        
        plt.figure()
        plt.plot(np.log10(timestep_nums), np.log10(timestep_num_iters))
        plt.xlabel('log10(Timestep)')
        plt.ylabel('log10(Num iterations in timestep)')
        plt.title(self.computedFilenameTitle() + ', Total iterations: ' + str(total_iters))
        plt.grid()
        plt.draw()
        
    def computedFilenameTitle(self):
        regex = '.*(\/[^\/]*\/[^\/]*\.hdf5)'
        m = re.match(regex, self.computedFilename)
        if m:
            return m.group(1)
        else:
            return self.computedFilename
        
    def read(self):
        #Read in
        #Set up regex strings
        
        
        lines = open(self.saveFile, 'r')
        
        self.exactFilename = lines.readline()
        self.computedFilename = lines.readline()
        self.timingFilename = lines.readline()
        self.time = float(lines.readline())
        self.timesteps = eval(lines.readline())
        
        regex = re.compile('(\d+)\s+\(([a-zA-Z_\s]+)\),\s+(.+)')
        #regex = re.compile('(\d+)\s+(.+)')
        for line in lines.readlines():
            m = regex.match(line)
            if m:
                var_i = int(m.group(1))
                var_name = m.group(2)
                errors = m.group(3)
                
                self.errors[var_i] = eval(errors)
            
        
#         print(self.errors)
            
        
    def write(self):
        #write to file
        f = open(self.saveFile, 'w')
        
        f.write(self.exactFilename + "\n")
        f.write(self.computedFilename + "\n")
        f.write(self.timingFilename + "\n")
        f.write(str(self.time) + "\n")
        f.write(str(self.timesteps))
        f.write("\n")
        
        for key in self.errors.keys():
            
            name = "Unknown var"
            if key in self.VAR_NAMES:
                name = self.VAR_NAMES[key]
                
            f.write(str(key) + " (" + name +"),  " + str(self.errors[key]) + "\n")
        
        
        f.close()