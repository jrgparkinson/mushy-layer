import os
import re
from shutil import copyfile
import subprocess
import sys
from MachineSpecific import MachineSpecific

machine_specific = MachineSpecific()


home_dir = machine_specific.get_home_dir()  #'/home/parkinsonjl/'
home_dir = '/home/parkinsonj/'

chombo_dir = os.path.join(home_dir, 'chombofork')
make_file = os.path.join(chombo_dir, 'lib/mk/Make.defs.local')
backup_make_file = os.path.join(chombo_dir, 'lib/mk/Make.defs.local.backup')
exec_dir = os.path.join(home_dir, 'convection-in-sea-ice/MushyLayer/execSubcycle')


# Backup existing file
#copyfile(make_file, backup_make_file)
command = "cp " + make_file + " " + backup_make_file
p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
for line in p.stdout.readlines():
    print line,
retval = p.wait()

print("Backup of old file complete")

try:

    print("Try to load new files")

    # Load original file and new file
    print("Try to load " + backup_make_file)
    F_backup = open(backup_make_file, "r")
    
    print("Try to load " + make_file)
    
    F_new = open(make_file, "w")
    
    print("Loaded new file to write to")

    old_file= F_backup.read()


    # Create new file and write it out
    new_file = re.sub("DIM.*2\n", 'DIM=3\n', old_file)
    F_new.write(new_file)

    # Close file handles
    F_new.close()
    F_backup.close()

    # Build code
    command = "cd " + exec_dir + "; make all"
    print(command)

    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    for line in p.stdout.readlines():
        print line,
    retval = p.wait()

except:
    print "Unexpected error:", sys.exc_info()[0]
    

# Copy back original file
#copyfile(backup_make_file, make_file)

command = "cp " + backup_make_file + " " + make_file
p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
for line in p.stdout.readlines():
    print line,
retval = p.wait()
