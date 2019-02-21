import getopt
import sys
import os
from mushyLayerRunUtils import read_inputs, get_final_chk_file, write_inputs,\
    get_executable

def create_refined_restart(argv):
    # Default vals:
    old_dir = ''
    new_dir = ''
    refinement = ''
    
    try:
        opts, args = getopt.getopt(argv, "n:p:r:")
    except getopt.GetoptError as err:
        print(str(err))
        print('create_refined_restart.py -n<new folder> -p<previous dir>'
              ' -r<refinement>')
        sys.exit(2)

    for opt, arg in opts:
        if opt in "-n":
            new_dir = str(arg)
        elif opt in "-p":
            old_dir = str(arg)
        elif opt in "-r":
            refinement = int(arg)
            
            
    prev_chk_file = os.path.join(old_dir, get_final_chk_file(old_dir))
    old_inputs_loc = os.path.join(old_dir, 'inputs')
    old_inputs = read_inputs(old_inputs_loc)
    new_box_size = int(old_inputs['main.max_grid_size']) * refinement
    
    out_file = os.path.join(new_dir, 'restart.2d.hdf5')


    new_inputs = {'inFile': prev_chk_file,
                  'run_inputs': old_inputs_loc,
                  'outFile': out_file,
                  'box_size': new_box_size,
                  'refinement': refinement}
                        
    new_inputs_loc = os.path.join(new_dir, 'inputsRefine')
    write_inputs(new_inputs_loc, new_inputs)
    
    # Run refine code
    exec_file = os.path.join(os.environ['MUSHY_LAYER_DIR'],
                                  'setupNewRun', get_executable('setupnewrun'))
   
    cmd = 'cd %s; %s %s' % (new_dir, exec_file, new_inputs_loc)
    print(cmd)
    os.system(cmd)


if __name__ == "__main__":
    create_refined_restart(sys.argv[1:])
