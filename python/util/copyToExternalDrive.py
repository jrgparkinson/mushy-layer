from timeout import timeout
import os, sys
import subprocess
from subprocess import Popen
import paramiko
import time
import re
from datetime import date, datetime

# Call this function from the command line like:
#  python copyToExternalDrive.py copy gyre4 mushyLayerLowC-periodic 60(minutes)
# determines how old folders must be to be copied. I.e. folders modified more
# recently than now - timeDelay will not be copied. This is to stop messing
# with runs that are currently in progress. 
#
# Example usage: 
#  python copyToExternalDrive.py copy gyre4 mushyLayerLowC-periodic 60
# python copyToExternalDrive.py copy ray mushyLayerLowC-periodic 120

# Note this useful command for backing up an entire directory: 
# rsync -avzr --exclude '.Trash-1001' --ignore-existing -e ssh parkinsonj@atmlxmaster.atm.ox.ac.uk:/network/group/aopp/oceans/AW002_PARKINSON_MUSH/optimalStates-highRes-new/ "/media/parkinsonjl/DATA/optimalStates-highRes-new/"
#  rsync -avzr --exclude '.Trash-1001' --update -e ssh parkinsonj@atmlxmaster.atm.ox.ac.uk:/network/group/aopp/oceans/AW002_PARKINSON_MUSH/ "/media/parkinsonjl/DATA/"
# rsync -avzr --exclude '.Trash-1001' --update -e ssh parkinsonj@atmlxmaster.atm.ox.ac.uk:/network/group/aopp/oceans/AW002_PARKINSON_MUSH/fixedChill-insulating/CR8.0RaC-1e4Le200KozenyPermeabilitypts256-0 "/media/parkinsonjl/DATA/fixedChill-insulating/CR8.0RaC-1e4Le200KozenyPermeabilitypts256-0"
# /home/parkinsonjl/mnt/sharedStorage/
# rsync -avzr --exclude '.Trash-1001' --ignore-existing -e ssh parkinsonj@atmlxmaster.atm.ox.ac.uk:/network/group/aopp/oceans/AW002_PARKINSON_MUSH/ "/media/parkinsonjl/DATA/"
remote_add = {'cloud':'some2962@cloud.atm.ox.ac.uk',
'gyre4': 'parkinsonj@gyre4.atm.ox.ac.uk',
'ray': 'parkinsonj@atmlxmaster.atm.ox.ac.uk',
'sharedStorage': 'parkinsonj@atmlxmaster.atm.ox.ac.uk'}
remote_base_dir = {'cloud':'/local/home/cloud/ice/Users/parkinson/',
'gyre4':'/local/home/gyre4/ice/Users/parkinsonj/convection-in-sea-ice/MushyLayer/run/',
'ray':'/home/parkinsonj/convection-in-sea-ice/MushyLayer/run/',
'sharedStorage': '/network/group/aopp/oceans/AW002_PARKINSON_MUSH/'}
hard_drive_dir = '/media/parkinsonjl/DATA/'
time_format = '%Y-%m-%d %H:%M'

# Copy all folders in dir which haven't been modified
# since timeDelay minutes ago to the cloud
def copyToDrive(server, folder, timeDelay):
    # get all folders within dir
    dir = remote_base_dir[server] + folder

    ls_command = 'ssh ' + remote_add[server] + ' cd ' + dir + '; ls -ld --time-style="+' + time_format + '" */'
    print ls_command
    ls = subprocess.Popen(['ssh', remote_add[server], 'cd ', dir, '; ls -ld --time-style="+' + time_format + '" */'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = ls.communicate()
    print out
    lines = out.split("\n")

    pattern = re.compile(".*(\d{4}-\d{2}-\d{2}\s+\d{2}:\d{2})\s+(.*)")
    #pattern = re.compile(".*(\d+)\s(.*)")

    foldersToCopy = []

    for line in lines:
        #print line
        # want to recover data modified and name
        
        m = re.match(pattern, line)

        if not m:
            continue

        date = m.group(1)
        date = date.replace('  ', ' ')
        folder_name = m.group(2)

        print(date + " " + folder_name)

        timestamp = time.mktime(datetime.strptime(date, time_format).timetuple())
        timestamp = timestamp/60 # get this in minutes


        now = time.time()/60
        if timestamp < (now - timeDelay):
            print(date + " " + folder_name + " " + str(timestamp) + " + offset <  " + str(now))
            foldersToCopy.append(folder_name)

        
    for f in foldersToCopy:
        # Check folder doesn't exist on hard drive
        folderName = f.replace('/','') # removing trailing slash
        copyTo = f
        i = 1
        while os.path.isdir(hard_drive_dir + folder + "/" + copyTo):
            copyTo = folderName + 'v' + str(i)
            i = i + 1

        full_remote_dir = remote_base_dir[server] + folder + "/" + f
        full_local_dir = hard_drive_dir + folder + "/" + copyTo
        # --remove-source-files 
        rsync_command = 'rsync -avz --ignore-existing -e ssh ' + \
                remote_add[server] + ':' + full_remote_dir + ' "' + full_local_dir + '"'
        print(rsync_command)                      

        exit_code = os.system(rsync_command)

        
    # Clearing up - remove empty remove folders 
    ls = subprocess.Popen(['ssh', remote_add[server], 'cd ', dir, "; find . -type d -empty"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = ls.communicate()
    #print out
    lines = out.split("\n")

    for line in lines:

        if line:
            delete_command = 'cd ' + dir + ' rm -rf ' + line
            print(delete_command)
            #ls = subprocess.Popen(['ssh', remote_add[server],  delete_command], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            #out, err = ls.communicate()
        
        

def folderExistsOnCloud(cloudFullDir):
    full_command = "ssh " + remote_add + " cd " + cloudFullDir
    print(full_command + '\n')

    command  = ["ssh", "%s" % remote_add, full_command]
    folderExists = True
        
    try:
        output = subprocess.check_output(command,
                                         stderr=subprocess.STDOUT,
                                         universal_newlines=True)
            
    except subprocess.CalledProcessError:
        folderExists = False
        
    return folderExists
    
if __name__ == '__main__':
    import argparse
    function_map = { 
        'copyToDrive': copyToDrive,
        'copy': copyToDrive,
    }
    
    parser = argparse.ArgumentParser()
    parser.add_argument( 'command', nargs=1 )
    parser.add_argument( 'server', nargs='+' )
    parser.add_argument( 'folder', nargs='+' )
    parser.add_argument( 'timeDelay', nargs='+' )
    args= parser.parse_args()
    function = function_map[args.command[0]]
    function( args.server[0], args.folder[0], float(args.timeDelay[0]) )
    
