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

    dataFolderNu = os.path.join(base_output_dir, 'ConvectionDB-cfl0.15','chi0.4-Da1.0e-02-Ra1.0e+05')
    dataFolderVariablePorosity = os.path.join(base_output_dir, 'FixedPorousHole-1proc')

    porousMushyHoleFolder = os.path.join(base_output_dir, 'PorousMushyHole-t0.0002')

    figureDirectory = base_output_dir

    noFlowData = os.path.join(base_output_dir, 'NoFlow')
    fixedChillData = os.path.join(base_output_dir, 'FixedChill-t1.0e-01-Ra1e+06-Da5.0e-04-C2.00-Rel1.0e-04-0')

    figCommands = []
    figCommands.append('Fig4NoFlow(\'' + noFlowData + '\', \'' + figureDirectory + ' \')')
    figCommands.append('Fig5FixedPorosityConvectionPlots(\'' + dataFolderNu + '\', \'' + dataFolderVariablePorosity + '\', \'' + figureDirectory + '\')')
    figCommands.append('Fig6PorousHole(\'' + porousMushyHoleFolder + '\', \'' + figureDirectory + '\')')
    figCommands.append('Fig7FixedChill(\'' + fixedChillData + '\', [3000, 4800, 17000], \'' + figureDirectory + '\')')

    full_matlab_command = matlab_command + ' "'

    for c in figCommands:
        full_matlab_command = full_matlab_command + 'try \n' + c + '; \n catch e \n fprintf(\'Plotting failed\'); \n end \n'

    full_matlab_command = full_matlab_command + ' exit;"'
 #  + makeFig5 + '; ' + makeFig6 + '; exit ;"'

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

