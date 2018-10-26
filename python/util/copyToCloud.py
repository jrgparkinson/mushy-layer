from timeout import timeout
import os, sys
import subprocess
from subprocess import Popen
import paramiko
import time

# Call this function from the command line like:
#  python copyToCloud.py copyToCloud directory timeDelay
# where directory is the one to copy to the cloud, and time delay (minutes)
# determines how old folders must be to be copied. I.e. folders modified more
# recently than now - timeDelay will not be copied. This is to stop messing
# with runs that are currently in progress. 
#
# Example usage: 
#  python copyToCloud.py copyToCloud ./mushyLayer-insulating 60

remote_add = 'some2962@cloud.atm.ox.ac.uk'
cloud_base_dir = '/local/home/cloud/ice/Users/parkinson/'

# Copy all folders in dir which haven't been modified
# since timeDelay minutes ago to the cloud
def copyToCloud(dir, timeDelay):
    # get all folders within dir
    dirs = [x[0] for x in os.walk(dir)]
    
    for dir in dirs:
        
        # Check this dir doesn't haven't any sub directories
        subDirs  = os.walk(dir).next()[1]
        
        if subDirs:
            continue
        
        # Check this dir wasn't modified too recently
        
        # Gives time since epoch directory was last modified
        folderTime = os.path.getmtime(dir)
        currentTime = time.time()
        
        # in seconds
        diff = currentTime - folderTime
        print(diff)
        if diff/60 < timeDelay:
            #print (str(diff/60) + " is less than " + str(timeDelay))
            continue
        
        # All good to go now
    
      
        firstSlash = dir.find('/')
        dir = dir[firstSlash+1:]
        
        cloudFullDir = cloud_base_dir + dir
        print(cloudFullDir)
 
        # check that folder doesn't exist on cloud by trying to cd to it
        folderExists = folderExistsOnCloud(cloudFullDir)
        
        if folderExists:
            i = 2
            while folderExists:
                newDir = cloudFullDir + 'v' + str(i)
                folderExists = folderExistsOnCloud(newDir)
                i = i + 1
                
            cloudFullDir = newDir
            print('New dir:' + cloudFullDir)
            
        # Before copying, ensure base directory exists
        thisDirParts = dir.split('/')
        parentDir = '/'.join(thisDirParts[:-1])
        
        parentDirs = [parentDir]
        
        while not folderExistsOnCloud(cloud_base_dir + parentDir):
            dirParts = parentDir.split('/')
            parentDir = '/'.join(dirParts[:-1])
            parentDirs.append(parentDir)
            
        for parDir in parentDirs:
            mkdirOnCloud(cloud_base_dir + parDir)
            
        
        print(parentDir)
        
        
        cloudParts = cloudFullDir.split('/')
        
        cloudParentDir = '/'.join(cloudParts[:-1])
                           
        rsync_command = 'rsync -avz --remove-source-files -e ssh ' + \
                 dir + ' '+ remote_add + ':' + cloudParentDir
                      
        exit_code = os.system(rsync_command)
        
        # rync will remove source files but not source directories.
        # Do this ourselves;
        
        # First double check directory is indeed empty
        files = os.listdir(dir)
        if len(files) == 0:
            print "Removing empty folder:", dir
            os.rmdir(dir)
    
    
        

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

@timeout(5)
def mkdirOnCloud(dir):
    os.system('ssh '+ remote_add + ' mkdir -p ' + dir)
    
    
if __name__ == '__main__':
    import argparse
    function_map = { 
        'copyToCloud': copyToCloud,
        'copy': copyToCloud,
    }
    
    parser = argparse.ArgumentParser()
    parser.add_argument( 'command', nargs=1 )
    parser.add_argument( 'dir', nargs='+' )
    parser.add_argument( 'timeDelay', nargs='+' )
    args= parser.parse_args()
    function = function_map[args.command[0]]
    function( args.dir[0], float(args.timeDelay[0]) )
    
