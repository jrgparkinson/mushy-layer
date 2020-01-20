import mushyLayerRunUtils as util
import os

mem_leak_path = os.path.join(util.get_mushy_layer_dir(), 'test', 'memory_leak')

# Create inputs file and directory
inputs_base = os.path.join(util.get_mushy_layer_dir(), 'test', 'regression', 'AMRDarcyBrinkman', 'inputs')
inputs = util.read_inputs(inputs_base)

inputs['main.max_step'] = 5
util.write_inputs('inputs', inputs)

exec_file = util.get_executable_name(exec_name='mushyLayer2d', return_full_path=True)

# , os.path.join(mem_leak_path, 'inputs')
cmd = 'valgrind --child-silent-after-fork=yes --xml=yes --xml-file=valgrind.xml %s inputs ' % (exec_file)
print(cmd)
os.system(cmd)
