import os
import sys
import subprocess
import getopt


def usage():
    print('python install.py -g <if from github actions>')


def print_var(var_name, is_git):
    if is_git:
        print("::set-env name=" + var_name + "::{}".format(os.environ[var_name]))
    else:
        print('%s=%s' % (var_name, os.environ[var_name]))


def get_script_path():
    return os.path.dirname(os.path.realpath(sys.argv[0]))


def program_exists(program):

    process = subprocess.run(['which', program],
                                stdout=subprocess.PIPE)
    if process.stdout.readline():
        print('Found %s' % program)
        return True
    else:
        return False


def setup(git_build=False):

    if git_build:
        print('Building from GitHub')

    # Add folders to PYTHONPATH
    required_paths = [os.path.join(get_script_path(), 'test'),
                      os.path.join(get_script_path(), 'plotting')]

    print('Checking PYTHONPATH')

    # Initialise PYTHONPATH if it doesn't exist
    if 'PYTHONPATH' not in os.environ:
        os.environ['PYTHONPATH'] = ''

    for required_path in required_paths:
        if required_path not in os.environ['PYTHONPATH']:
            os.environ['PYTHONPATH'] += ':' + required_path

    print_var('PYTHONPATH', git_build)

    # Setup MUSHY_LAYER_DIR
    print('Checking MUSHY_LAYER_DIR')
    if 'MUSHY_LAYER_DIR' not in os.environ:
        os.environ['MUSHY_LAYER_DIR'] = get_script_path()
    print_var('MUSHY_LAYER_DIR', git_build)

    if git_build:
        return

    # Check for required software
    required_commands = ['perl', 'csh']

    for cmd in required_commands:
        try:
            process = subprocess.run(['which', cmd],
                             stdout=subprocess.PIPE)
            output = process.stdout.readline()
            if output:
                print('Found %s' % cmd)
            else:
                print('Cannot find command: %s' % cmd)
        except FileNotFoundError:
            print('Cannot find program: %s' % cmd)
            sys.exit(-1)

    # Also need a c++ compiler, and some sort of make
    cpp_progs = ['g++', 'mpiCC', 'icc', 'xIC']
    found_cpp_progs = [p for p in cpp_progs if program_exists(p)]
    if not found_cpp_progs:
        print('Warning - could not find a c++ compiler (looked for %s)' % ','.join(cpp_progs))

    make_progs = ['make', 'cmake', 'gmake']
    if not [program_exists(p) for p in make_progs]:
        print('Warning - could not find a make program (looked for %s)' % ','.join(make_progs))

    fortran_progs = ['gfortran', 'g95', 'g97', 'f90', 'f95', 'f77', 'pgf77', 'pgf90', 'ifc']
    found_fortran = [f for f in fortran_progs if program_exists(f)]
    if not found_fortran:
        print('Warning - could not find a Fortran compiler (looked for %s)' % ','.join(fortran_progs))

    lib_paths = []
    for environmental_vars in ['LD_LIBRARY_PATH', 'PATH']:
        if environmental_vars in os.environ:
            lib_paths += os.environ[environmental_vars].split(':')

    hdf5_path = False
    for lib in lib_paths:
        if 'hdf5' in lib and os.path.exists(lib):
            hdf5_path = lib
            break

    if hdf5_path:
        print('Found hdf5: %s ' % hdf5_path)
    else:
        print('Cannot find hdf5. If you have already installed it, please ensure it is on your LD_LIBRARY_PATH')
        install_hdf5 = input('Download and install HDF5 now? (y/n) ')
        if install_hdf5 == 'y':
            parent_folder = os.path.dirname(os.environ['MUSHY_LAYER_DIR'])
            hdf5_path = os.path.join(parent_folder, 'hdf5')
            install_path_response = input('Enter install path, or hit enter to accept [%s]' % hdf5_path)
            if install_path_response:
                hdf5_path = install_path_response

            url= 'https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.21/src/hdf5-1.8.21.tar.gz'
            subprocess.check_call(['wget', url])
            subprocess.check_call(['tar', '-zxvf', 'hdf5-1.8.21.tar.gz'])
            subprocess.check_call(['cd', 'hdf5-1.8.21'])
            subprocess.check_call(['./ configure', '--prefix=', hdf5_path])
            subprocess.check_call(['make'])
            subprocess.check_call(['make', 'install'])

            os.environ['LD_LIBRARY_PATH'] += ':' + hdf5_path
            print_var('LD_LIBRARY_PATH', git_build)

    shell = os.environ['SHELL']
    print('Your shell is: %s, add environmental variables to '
          '~/.%src for future ease of use' % (shell, shell.replace('/bin/', '')))

    # Check llapack
    # locate libblas.so
    installed = False
    proc = subprocess.run(['apt-cache', 'policy', 'liblapack3'])
    for line in iter(proc.stdout.readline, ''):
        if 'Installed' in line.decode():
            installed = True

    if installed:
        print('Found liblapack3')
    else:
        print('Could not find lapack/lblas installation')

    # Construct example options
    options = {'DIM': 2,
               'DEBUG': False,
               'OPT': True,
               'CXX': found_cpp_progs[0],
               'FC': found_fortran[0],
               'MPI': False,
               'USE_HDF': True,
               'HDFINCFLAGS': '-I' + hdf5_path +'/include',
               'HDFLIBFLAGS': '-L/' + hdf5_path + '/lib - lhdf5 - lz',
               'flibflags': ['-lblas', '-llapack', '-lgfortran'],
               'LAPACKLIBS': '-L/usr/lib -L/usr/lib/lapack -L/usr/lib/libblas -llapack -lblas',
               }


if __name__ == "__main__":

    git_build = 'GITHUB_WORKSPACE' in os.environ

    try:
        opts, args = getopt.getopt(sys.argv[1:], "g")
    except getopt.GetoptError as err:
        # print help information and exit:
        print(str(err))  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    for o, a in opts:
        if o == "-g":
            git_build = True

    setup(git_build)
