# Class to load a checkpoint file and perform operations on it
# e.g. find where brine channels are
import h5py
import numpy as np
#import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import pyplot
import copy
import os
from mushyLayerRunUtils import readInputs

NUM_COMPS = 'num_comps'
DATA = 'data'
DX = 'dx'
DT = 'dt'
REF_RATIO = 'ref_ratio'
BOXES = 'boxes'

class ChkFile:
    

    def __init__(self, filename):
        self.filename = filename
        
        self.directory = os.path.dirname(filename)
        self.inputs = readInputs(os.path.join(self.directory, 'inputs'))
        
        h5File = h5py.File(self.filename, 'r') 

        #print('Loaded HDF5 file with groups: ' + str(h5File.keys()))
        #print(' and attributes: ' + str(h5File.attrs.keys()))

      

        # Get all the different bits of data from the file
        #print(h5File)

        ChomboGroup = h5File['Chombo_global']
        globalAttrs = ChomboGroup.attrs

        
        attrs = h5File.attrs

        self.time = attrs['time']
        self.iteration = int(attrs['iteration'])
        self.max_level = int(attrs['max_level'])
        self.num_levels = int(attrs['num_levels'])
        self.regrid_interval = int(attrs['regrid_interval_0'])
        self.num_comps = int(attrs['num_components'])
        self.spaceDim = int(globalAttrs['SpaceDim'])

        # Now read all the component names
        self.data = {}
        allNames = attrs.values()
        for i in range(0, self.num_comps):
            name = attrs['component_' + str(i)]
            # Work out if it's a vector or scalar
            if name[0] == 'x' or name[0] == 'y' or name[0] == 'z' and name[1:] in allNames:
                # Vector. Don't add if we've already got it
                if name[1:] not in self.data.keys():
                    self.data[name[1:]] = {NUM_COMPS: 1}
                else:
                    self.data[name[1:]][NUM_COMPS] = self.data[name[1:]][NUM_COMPS] + 1
            else:
                self.data[name] = {NUM_COMPS: 1}

        #print(self.data)

        grp0 = h5File['level_0']
        #print('Group level_0 has keys: ' + str(grp0.keys()))
        #print('  and attributes: ' + str(grp0.attrs.keys()))

        self.levels = [None] * self.num_levels
        for level in range(0, self.num_levels):
            levelGroup = h5File['level_' + str(level)]
            groupAtts = levelGroup.attrs
            boxes = levelGroup['boxes']
            self.levels[level] = {DX: groupAtts['dx'], DT: groupAtts['dt'],
                                  REF_RATIO: groupAtts['ref_ratio'], BOXES: list(boxes)}

            # Some attributes on level 0 apply to whole hierarchy
            if level == 0:
                self.time = groupAtts['time']
                self.isPeriodic = [None] * self.spaceDim
                for dim_i in range(0, self.spaceDim):
                    self.isPeriodic[dim_i] = groupAtts['is_periodic_' + str(dim_i)]

                self.tag_buffer_size = groupAtts['tag_buffer_size'] 
                self.prob_domain = groupAtts['prob_domain']
                print('Problem domain: ' + str(self.prob_domain))
        

        # Now get the actual data for each field, on each level
        print(self.data.keys())
        for compName in self.data.keys():
            self.data[compName][DATA] = [None] * self.num_levels
            for level in range(0, self.num_levels):
                levelGroup = h5File['level_' + str(level)]
                #print(compName + ':offsets=0')
                offsets = levelGroup[compName + ':offsets=0']
                data = levelGroup[compName + ':datatype=0']
                data_atts = levelGroup[compName + '_attributes']



                # Reshape data
                dataUnshaped = data[()]

                shapedData = []
                # First work out what the shape should be!
                #print(self.levels[level][BOXES])
                offset = 0
                for box in self.levels[level][BOXES]:
                    #print(box)
                    # Box = [lo_i lo_j hi_i hi_j]
                    lo_i = box[0]
                    lo_j = box[1]
                    hi_i = box[2]
                    hi_j = box[3]
 
                    numRows = hi_j + 1 - lo_j
                    numCols = hi_i + 1 - lo_i
                    numCells = numRows*numCols
                    #print(str(numCells))
                    dataBox = dataUnshaped[offset:offset + numCells]
                    
                    reshapedData = dataBox.reshape((numRows, numCols))
                    shapedData.append(reshapedData)


                    offset = offset + numCells
                

                # Should really get data into a nice format like a np array
                self.data[compName][DATA][level] = shapedData
           

        #print(self.levels[0])
        
        self.computeDiagnosticVars()
         
               
        h5File.close()
        
        
    def computeDiagnosticVars(self):
        # First get parameters
        CR = float(self.inputs['parameters.compositionRatio'])
        St = float(self.inputs['parameters.stefan'])
        cp = float(self.inputs['parameters.specificHeatRatio'])
        pc = float(self.inputs['parameters.waterDistributionCoeff'])
        thetaEutectic = 0.0
        ThetaEutectic = 0.0
        
        H = self.getData('Enthalpy')
        S = self.getData('Bulk concentration')
        
        
        He = np.empty(H.shape)
        Hs = np.empty(H.shape)
        Hl = np.empty(H.shape)
        
        porosity = np.empty(H.shape)
        temperature = np.empty(H.shape)
        liquidSalinity = np.empty(H.shape)
        solidSalinity = np.empty(H.shape)
        
        # Now compute bounding energies
        cols = H.shape[0]
        rows = H.shape[1]
        for j in xrange(cols):
            
            for i in xrange(rows):
                eutecticPorosity = (CR + S[j, i])/(ThetaEutectic + CR)
                He[j, i] = eutecticPorosity*(St + thetaEutectic*(1-cp)) + cp*thetaEutectic
                Hs[j, i] = cp*(thetaEutectic + max(0.0, (-S[j, i] - CR)/pc))
                Hl[j, i] = St - S[j, i] + thetaEutectic + ThetaEutectic
                
        # Compute diagnostic variables      
        for j in xrange(cols):
            
            for i in xrange(rows):
                if H[j, i] <= Hs[j, i]:
                    porosity[j, i] = 0.0;
                    temperature[j, i] = H[j, i]/cp;
                    liquidSalinity[j, i] = 0.0
                    solidSalinity[j, i] = S[j, i]
                elif H[j, i] <= He[j, i]:
                    porosity[j, i] = (H[j, i] - thetaEutectic*cp)/(St + thetaEutectic*(1-cp))
                    temperature[j, i] = thetaEutectic
                    liquidSalinity[j, i] = ThetaEutectic
                    solidSalinity[j, i] = S[j, i]/(1-porosity[j, i])
                elif H[j, i] < Hl[j, i]:
                    a = CR*(cp-1) + St*(pc-1)
                    b = CR*(1-2*cp) + H[j, i]*(1-pc) - S[j, i] * (cp-1) - pc*St
                    c = (S[j, i] + CR)*cp + pc*H[j, i]
                    porosity[j, i] = (-b - np.sqrt(b*b-4*a*c))/(2*a)
                    
                    liquidSalinity[j, i] = (S[j, i] + CR*(1-porosity[j, i]))/(porosity[j, i] + pc*(1-porosity[j, i]))
                    temperature[j, i] =  -liquidSalinity[j, i] 
                    solidSalinity[j, i] = (pc*S[j, i] -CR*porosity[j, i])/(porosity[j, i] + pc*(1-porosity[j, i]))
                else:
                    porosity[j, i] = 1.0
                    temperature[j, i] = H[j, i]- St
                    liquidSalinity[j, i] = S[j, i]
                    solidSalinity[j, i] = 0.0
                    
        
        # Finally, save data into boxes
        self.saveData('Porosity', porosity)
        self.saveData('Temperature', temperature)
        self.saveData('LiquidSalinity', liquidSalinity)
        self.saveData('SolidSalinity', solidSalinity)
        
        
    def saveData(self, varName, lev0dat):
        width = self.prob_domain[2] +1 - self.prob_domain[0]
        height = self.prob_domain[3]  +1 - self.prob_domain[1]

        # If data doesn't exist, make it
        if not varName in self.data:
            H = self.data['Enthalpy']
            deepCopy = copy.deepcopy(H)
            self.data[varName] = deepCopy

        for box_i in range(0, len(self.levels[0]['boxes'])):
            box = self.levels[0]['boxes'][box_i]
    
            #print(box)
            lev0datBox = self.data[varName]['data'][0][box_i]
    
            ioffset = box[0]-self.prob_domain[0]
            joffset = box[1]-self.prob_domain[1]
    
            for i in range(0, box[2]+1-box[0]):
                for j in range(0, box[3]+1-box[1]):
                    #print(str(i) + ', ' + str(j))
                    #val = porosityBox[j][i]
                    lev0datBox[j][i] = lev0dat[j+joffset][i+ioffset]
                    #porosity[j+joffset][i+ioffset] = porosityBox[j][i]
        
    def getData(self, varName):
        # Reconstruct level 0 data as single np array
        width = self.prob_domain[2] +1 - self.prob_domain[0]
        height = self.prob_domain[3]  +1 - self.prob_domain[1]

        lev0Dat = np.empty([height, width])

        for box_i in range(0, len(self.levels[0]['boxes'])):
            box = self.levels[0]['boxes'][box_i]
    
            #print(box)
            lev0DatBox = self.data[varName]['data'][0][box_i]
    
            ioffset = box[0]-self.prob_domain[0]
            joffset = box[1]-self.prob_domain[1]
    
            for i in range(0, box[2]+1-box[0]):
                for j in range(0, box[3]+1-box[1]):
                    #print(str(i) + ', ' + str(j))
                    #val = lev0DatBox[j][i]
                    lev0Dat[j+joffset][i+ioffset] = lev0DatBox[j][i]
                    
                    
        return lev0Dat
        
        
    def channelProperties(self, doPlots = False):
        porosity = self.getData('Porosity')
        
        #width = self.prob_domain[2] +1 - self.prob_domain[0]
        width = porosity.shape[1]
        print('Domain width ' + str(width))

        # Iterate over porosity field 
        cols = porosity.shape[0]
        rows = porosity.shape[1]
        channels = [None] * cols
        chimneyPositions = []
        for j in xrange(cols):
            # chimneys in row
            chimneysInRow = 0
            averagePorosity = 0
            
            currentlyLiquid = False
            chimneyPosRow = []
    
            for i in xrange(rows):
                #print porosity[i,j]
                chi = porosity[j, i]
                averagePorosity = averagePorosity + chi
                if chi < 1.0 and currentlyLiquid:
                    # Just left liquid region
                    currentlyLiquid = False
                    
                elif chi > 0.999 and not currentlyLiquid:
                    # Entered a liquid region
                    chimneysInRow = chimneysInRow + 1
                    chimneyPosRow.append(i)
                    currentlyLiquid = True
                    
            averagePorosity = averagePorosity/rows
            if averagePorosity > 0.99:
                channels[j] = 0 # This region was entirely liquid, can't have a channel
            else:     
            	channels[j] = chimneysInRow
            	
            	
    	
    	        # First add chimney positions for chimneys we've already found
    	        #print(str(chimneyPositions))
    	        #print(str(len(chimneyPositions)))
    	
    	        if chimneysInRow == 0:
    	            continue
    	    
    	        for chimney_i in range(0, len(chimneyPositions)):
                    if chimney_i < len(chimneyPosRow):
    	                chimneyPositions[chimney_i].append(chimneyPosRow[chimney_i])
    	    
    	        # Now, add extra rows to chimneyPositions vector if needed
    	        for chimney_i in range(len(chimneyPositions), chimneysInRow):
                    if chimney_i < len(chimneyPosRow):
    	                chimneyPositions.append([chimneyPosRow[chimney_i]])
    
        print(str(chimneyPositions))
        # Filter out short channels
        chimneyLengths = []
        for chimney_i in range(0, len(chimneyPositions)):
            chimneyLengths.append(float(len(chimneyPositions[chimney_i])))
          
            
        avChimLength = np.mean(chimneyLengths)
        longEnoughChimneys = []
        for chimney_i in range(0, len(chimneyPositions)):
            if len(chimneyPositions[chimney_i]) > avChimLength/2:
                longEnoughChimneys.append(chimneyPositions[chimney_i])
            

        print(str(longEnoughChimneys))

        averageChanPositions = []
        relChanPositions = []
        for i in range(0, len(longEnoughChimneys)):
            averageChanPositions.append(np.mean(longEnoughChimneys[i]))
            relChanPositions.append(averageChanPositions[-1]/width)
    
        numChannels =len(longEnoughChimneys)
        #print('Number of channels: ' + str(numChannels))
        #print('Channel positions: ' + str(averageChanPositions))
        #print('Relative channel positions: ' + str(relChanPositions))     
        
        

        #doPlots = False
        if doPlots:
            cmap2 = mpl.colors.LinearSegmentedColormap.from_list('my_colormap',
                                           ['blue','black','red'],
                                           256)

            # tell imshow about color map so that only set colors are used
            img = pyplot.imshow(porosity,interpolation='nearest',
                    cmap = cmap2,origin='lower')

            # make a color bar
            pyplot.colorbar(img,cmap=cmap2 )

            pyplot.show() 
            
        return [numChannels, relChanPositions]
