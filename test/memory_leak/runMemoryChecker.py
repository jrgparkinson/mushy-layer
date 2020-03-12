import mushyLayerRunUtils as util
import os
from valgrind_parser import ValgrindLogParser
import shutil

mem_leak_path = os.path.join(util.get_mushy_layer_dir(), 'test', 'memory_leak')
VALGRIND_LOG = os.path.join(mem_leak_path, 'valgrind_log.txt')

shutil.copy('valgrind_regexes.json', '../../env/lib/python3.6/site-packages/valgrind_parser/data/')

# Create inputs file and directory
inputs_base = os.path.join(util.get_mushy_layer_dir(), 'test', 'regression', 'AMRDarcyBrinkman', 'inputs')
inputs = util.read_inputs(inputs_base)

inputs['main.max_step'] = 5
util.write_inputs('inputs', inputs)

exec_file = util.get_executable_name(exec_name='mushyLayer2d', return_full_path=True)

# , os.path.join(mem_leak_path, 'inputs')
#  --log-file="%s"
cmd = 'valgrind --child-silent-after-fork=yes --xml=yes --xml-file=valgrind.xml --leak-check=full ' \
      '--show-reachable=yes %s inputs > %s ' % (exec_file, VALGRIND_LOG)
# print(cmd)
os.system(cmd)

vlp = ValgrindLogParser(VALGRIND_LOG, html_report_location='valgrind.html')
vlp.generate_html_report()
