import os
import sys
import subprocess
import getopt

def usage():
    print('python install.py -g <if from github actions>')


def print_var(var_name, is_git):
    if is_git:
        print("::set-env name=" + var_name + ":{}".format(os.environ[var_name]))
    else:
        print('%s=%s' % (var_name, os.environ[var_name]))


def get_script_path():
    return os.path.dirname(os.path.realpath(sys.argv[0]))

def program_exists(program):

    process = subprocess.Popen(['which', program],
                                stdout=subprocess.PIPE)
    if process.stdout.readline():
        print('Found %s' % program)
        return True
    else:
        return False

if __name__ == "__main__":

    git_build = False
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

    # Add folders to PYTHONPATH
    required_paths = [os.path.join(get_script_path(), 'test'),
                      os.path.join(get_script_path(), 'plotting')]

    for required_path in required_paths:
        if 'PYTHONPATH' not in os.environ or required_path not in os.environ['PYTHONPATH']:
            os.environ['PYTHONPATH'] += ':' + required_path
    print_var('PYTHONPATH', git_build)

    # Setup MUSHY_LAYER_DIR
    if 'MUSHY_LAYER_DIR' not in os.environ:
        os.environ['MUSHY_LAYER_DIR'] = get_script_path()
    print_var('MUSHY_LAYER_DIR', git_build)

    if git_build:
        sys.exit(-1)

    # Check for required software
    required_commands = ['perl', 'csh']

    for cmd in required_commands:
        try:
            process = subprocess.Popen(['which', cmd],
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
    if not [program_exists(p) for p in cpp_progs]:
        print('Warning - could not find a c++ compiler (looked for %s)' % ','.join(cpp_progs))

    make_progs = ['make', 'cmake', 'gmake']
    if not [program_exists(p) for p in make_progs]:
        print('Warning - could not find a make program (looked for %s)' % ','.join(make_progs))

    fortran_progs = ['gfortran', 'g95', 'g97', 'f90', 'f95', 'f77', 'pgf77', 'pgf90', 'ifc']
    if not [program_exists(p) for p in fortran_progs]:
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

    hdf5_path = False

    if hdf5_path:
        print('Found hdf5: %s ' % hdf5_path)
    else:
        install_hdf5 = input('Cannot find hdf5. Download and install now? (y/n) ')
        if install_hdf5 == 'y':
            parent_folder = os.path.dirname(os.environ['MUSHY_LAYER_DIR'])
            install_path = os.path.join(parent_folder, 'hdf5')
            install_path_response = input('Enter install path, or hit enter to accept [%s]' % install_path)
            if install_path_response:
                install_path = install_path_response

            url= 'https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.21/src/hdf5-1.8.21.tar.gz'
            response = subprocess.Popen(['wget', url])
            response = subprocess.Popen(['tar', '-zxvf', 'hdf5-1.8.21.tar.gz'])
            response = subprocess.Popen(['cd', 'hdf5-1.8.21'])
            response = subprocess.Popen(['./ configure', '--prefix=', install_path])
            response = subprocess.Popen(['make'])
            response = subprocess.Popen(['make', 'install'])


