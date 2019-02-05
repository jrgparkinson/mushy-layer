# Run all convergence tests for the methods paper
import os
from colorama import Fore, Style
from mushyLayerRunUtils import get_base_output_dir, get_matlab_base_command
from SlurmTask import SlurmTask

##########################################################################
# Make figures for methods paper
##########################################################################

def make_figures():

    base_output_dir = get_base_output_dir()
    matlab_command = get_matlab_base_command()

    print(Fore.GREEN + 'Make figures' + Style.RESET_ALL)

    data_folder_nu = os.path.join(base_output_dir, 'ConvectionDB-cfl0.15','chi0.4-Da1.0e-02-Ra1.0e+05')
    data_folder_variable_porosity = os.path.join(base_output_dir, 'FixedPorousHole-1proc')

    porous_mushy_hole_folder = os.path.join(base_output_dir, 'PorousMushyHole-t0.00015')

    figure_directory = base_output_dir

    no_flow_data = os.path.join(base_output_dir, 'NoFlow')
    fixed_chill_data = os.path.join(base_output_dir, 'FixedChill-t1.0e-01-Ra1e+06-Da5.0e-04-C2.00-Rel1.0e-04-0')

    fig_commands = ['Fig4NoFlow(\'' + no_flow_data + '\', \'' + figure_directory + 'Fig4BenchmarkNoFlow.pdf\')',
                   'Fig5FixedPorosityConvectionPlots(\'' + data_folder_nu + '\', \'' + data_folder_variable_porosity + '\', \'' + figure_directory + '\')',
                   'Fig6PorousHole(\'' + porous_mushy_hole_folder + '\', \'' + figure_directory + '\')',
                   'Fig7FixedChill(\'' + fixed_chill_data + '\', [3000, 4000, 7000], \'' + figure_directory + '\')',
                   'Fig7MakeVideo(\'' + fixed_chill_data + '\')']

    full_matlab_command = ''

    for c in fig_commands:
        full_matlab_command = full_matlab_command + 'try \n' + c + ';  \ncatch e \n fprintf(\'Plotting failed\'); \nend \n'

   
    # Write out as a script then run that
    scriptFile = os.path.join(base_output_dir, 'makeFigureScript.m')
    f = open(scriptFile, "w")
    full_matlab_command = "set(groot, 'defaultAxesFontName', 'Times'); \nset(groot, 'defaultTextFontName', 'Times'); \n" + full_matlab_command
    f.write(full_matlab_command)
    f.close()

    slurm_command = '  matlab -nodisplay -nosplash -nodesktop -r " makeFigureScript; exit; "'
    #slurmCommand = matlab_command + ' "' + full_matlab_command + ' exit;"'
 #  + makeFig5 + '; ' + makeFig6 + '; exit ;"'

    job_name = 'makeFigures'

    s = SlurmTask(base_output_dir, job_name, '', 4)

    #s.setDependency(job_ids)
    s.set_custom_command(slurm_command)

    s.write_slurm_file(job_name + '.sh')
    s.run_task(job_name + '.sh')
    print(Fore.GREEN + 'Submitted make figures job \n' + Fore.RESET)


if __name__ == "__main__":
    #testPorousMushyHole(sys.argv[1:])
    make_figures()

