from subprocess import Popen

BM1 = 1
BM2 = 2
FULL_PROBLEM = 3

exec_dir = {BM1: '/benchmark/bm1_src/', BM2: '/benchmark/bm2_src/',
            FULL_PROBLEM: '/src/'}

output_base_dir = {BM1: './katz-bm1/', BM2: './katz-bm2/', FULL_PROBLEM: './katz/'}

def mushylayer_run(output_dir, problem_type, params):
    program = '../../mushylayer' + exec_dir[problem_type] + 'mushyLayer'
    command = (program +
               ' -output_file ' + output_dir + 'out ')

    for key in params.keys():
        command = command + ' -' + key + ' ' + str(params[key]) + ' '
            #  ' -t_max ' + str(t_max) +
            #  ' -t_out ' + str(t_out) +
            #  ' -ni ' + str(ni) +
            # ' -nj ' + str(nj)
            # )

    print(command + "\n")

    process = Popen(command, shell=True)

    exit_code = process.wait()
    return exit_code


params = {}

# First run for figure 6
problem_type = FULL_PROBLEM
output_dir = output_base_dir[problem_type] + 'fig6/'
params['t_max'] = 120000  # max model time (seconds)
params['t_out'] = 30  # model output interval (seconds)

params['ni'] = 120
params['nj'] = 80

#mushylayer_run(output_dir, problem_type, params)

# Now try and recreate figure 9
params['Pi0'] = 1e-10
params['Da'] = 0.1

output_dir = output_base_dir[problem_type] + 'fig9-5/'

mushylayer_run(output_dir, problem_type, params)
