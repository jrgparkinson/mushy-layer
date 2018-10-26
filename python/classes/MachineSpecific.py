import socket


class MachineSpecific:
    # Constants to specify what machine we're running the code on
    GYRE1 = 1
    GYRE2 = 2
    GYRE3 = 3
    GYRE4 = 4
    RAY1  = 5
    RAY2  = 6
    RAYMASTER = 7
    LAPTOP = 10
    OSX = 11

    CLOUD = 20
    UNKNOWN = -1
                
    homeDirs = {GYRE1: '/local/home/gyre1/ice/Users/parkinsonj',
                GYRE2: '/local/home/gyre2/ice/Users/parkinsonj',
                GYRE3: '/local/home/gyre3/ice/Users/parkinsonj',
                GYRE4: '/local/home/gyre4/ice/Users/parkinsonj',
                RAY1: '/local/home/parkinsonj',
                RAY2: '/local/home/parkinsonj',
                RAYMASTER: '/home/parkinsonj',
                LAPTOP: '/home/parkinsonjl',
                OSX: '/Users/parkinsonj',
                UNKNOWN: '~'}
                
              
    OPTIONS = 'gfortran.OPT.MPI'
    #OPTIONS = 'gfortran.DEBUG.MPI'
    
    programNames = {GYRE1: 'mushyLayer2d.Linux.64.mpiCC.'+OPTIONS+'.ex',
                    GYRE2: 'mushyLayer2d.Linux.64.mpiCC.'+OPTIONS+'.ex',
                    GYRE3: 'mushyLayer3d.Linux.64.mpiCC.'+OPTIONS+'.ex',
                    GYRE4: 'mushyLayer2d.Linux.64.mpiCC.'+OPTIONS+'.ex',
                    RAY1: 'mushyLayer2d.Linux.64.mpiCC.'+OPTIONS+'.ex',
                    RAY2: 'mushyLayer2d.Linux.64.mpiCC.'+OPTIONS+'.ex',
                    RAYMASTER: 'mushyLayer2d.Linux.64.mpiCC.'+OPTIONS+'.ex',
                    LAPTOP: 'mushyLayer2d.Linux.64.mpiCC.'+OPTIONS+'.ex',
                    #LAPTOP: 'mushyLayer2d.Linux.64.mpiCC.gfortran.DEBUG.MPI.ex',
                    OSX: 'mushyLayer2d.Darwin.64.g++.gfortran-mp-5.OPTHIGH.ex',
                    UNKNOWN: 'mushyLayer2d.Linux.64.mpiCC.'+OPTIONS+'.ex'}
    
    localDirs = {GYRE1: '/local/home/',
                GYRE2: '/local/home/',
                GYRE3: '/local/home/',
                GYRE4: '/local/home/',
                RAY1: '/local/home/',
                RAY2: '/local/home/',
                RAYMASTER: '/home/',
                LAPTOP: '/local',
                OSX: '/local',
                UNKNOWN: '/local'
                 }
    
    machineNames = {GYRE1: 'gyre1',
                GYRE2: 'gyre2',
                GYRE3: 'gyre3',
                GYRE4: 'gyre4',
                RAY1: 'atmlxint1',
                RAY2: 'atmlxint2',
                RAYMASTER: 'atmlxmaster',
                LAPTOP: 'Lenovo Thinkpad',
                OSX: 'Apple desktop',
                UNKNOWN: 'Unknown'
                    }

    def __init__(self):
        hostname = socket.gethostname()


        if (hostname == 'atmlxlap005'):
            self.machine = self.LAPTOP
        elif (hostname == 'gyre1'):
            self.machine = self.GYRE1
        elif (hostname == 'gyre2'):
            self.machine = self.GYRE2
        elif (hostname == 'gyre3'):
            self.machine = self.GYRE3
        elif (hostname == 'gyre4'):
            self.machine = self.GYRE4
        elif (hostname == 'atmlxint1'):
            self.machine = self.RAY1
        elif (hostname == 'atmlxint2'):
            self.machine = self.RAY2
        elif ('atmnode' in hostname) or 'atmlxmaster' == hostname:
            self.machine = self.RAYMASTER
        elif (hostname == 'atmosx11'):
            self.machine = self.OSX
        elif (hostname == 'cloud'):
            self.machine = self.CLOUD
        else:
            self.machine = self.UNKNOWN

    def get_home_dir(self):
        return self.homeDirs[self.machine]

    def get_program_name(self):
        return self.programNames[self.machine]
    
    def get_machine_name(self):
        return self.machineNames[self.machine]
    
    def local_dir(self):
        return self.localDirs[self.machine]

    def has_dropbox(self):
        if self.machine == self.LAPTOP:
            return True
        else:
            return False

    def is_gyre(self):
        if (self.machine == self.GYRE1 or self.machine == self.GYRE2 or
                self.machine == self.GYRE3 or self.machine == self.GYRE4):
            return True
        else:
            return False
