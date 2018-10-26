import os, sys

LOCAL_ROOT = "/media/parkinsonjl/FREECOM HDD/"
REMOTE_ROOT = "/network/group/aopp/oceans/AW002_PARKINSON_MUSH/"


run_type = "optimalStates-highRes"
#run_type = "optimalStates-lowRes"
#run_type = "channelSpacing-periodic"
#run_type = "fixedChill-insulating"
#run_type = "darcyBrinkman-insulating"

local_full_path = os.path.join(LOCAL_ROOT, run_type)
remote_full_path = os.path.join(REMOTE_ROOT, run_type)

folders = os.walk(local_full_path).next()[1]

print(folders)

for folder in folders:
    cmd =  'cd "' + local_full_path + '"; rsync -avz --ignore-existing -e ssh "' + folder + '" raymaster:' + remote_full_path
    exit_code = os.system(cmd)
    
    # Dry run:
    #print(cmd)
    #exit_code = 0
    
    
    print('Synced ' + folder + ', with exit status ' + str(exit_code))
               
