# Run all convergence tests for the methods paper
import os, sys
from colorama import Fore, Style
import getopt

from runAMRConvergenceTest import runTest
from mushyLayerRunUtils import getBaseOutputDir, getMatlabBaseCommand
from SlurmTask import SlurmTask

##########################################################################
# Make figures for methods paper
##########################################################################

def makeFigures():

    base_output_dir = getBaseOutputDir()
    matlab_command = getMatlabBaseCommand()

    print(Fore.GREEN + 'Make figures' + Style.RESET_ALL)

    dataFolderNu = os.path.join(base_output_dir, '/ConvectionDB-cfl0.15/chi0.4-Da1.0e-02-Ra1.0e+05/')
    dataFolderVariablePorosity = os.path.join(base_output_dir, '/FixedPorousHole-1proc/')

    porousMushyHoleFolder = os.path.join(base_output_dir, '/PorousMushyHole-t0.0005/')

    figureDirectory = base_output_dir

    makeFig5 = 'Fig5FixedPorosityConvectionPlots(\'' + dataFolderNu + '\', \'' + dataFolderVariablePorosity + '\', \'' + figureDirectory + '\')'
    makeFig6 = 'Fig6PorousHole(\'' + porousMushyHoleFolder + '\', \'' + figureDirectory + '\')'

    full_matlab_command = matlab_command + ' "' + makeFig5 + '; ' + makeFig6 + '; exit ;"'

    jobName = 'makeFigures'

    s = SlurmTask(base_output_dir, jobName, '', 4)

    #s.setDependency(job_ids)
    s.setCustomCommand(full_matlab_command)

    s.writeSlurmFile(jobName + '.sh')
    s.runTask(jobName + '.sh')
    print(Fore.GREEN + 'Submitted make figures job \n' + Fore.RESET)


if __name__ == "__main__":
    #testPorousMushyHole(sys.argv[1:])
    makeFigures()

