# Class to load a checkpoint file and perform operations on it
# e.g. find where brine channels are
import h5py
import numpy as np
#import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import pyplot



class PltFile:
    
    NUM_COMPS = 'num_comps'
    DATA = 'data'
    DX = 'dx'
    DT = 'dt'
    REF_RATIO = 'ref_ratio'
    BOXES = 'boxes'
    

    def __init__(self, filename):
        self.filename = filename
        
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
        #self.regrid_interval = int(attrs['regrid_interval_0'])
        self.num_comps = int(attrs['num_components'])
        self.spaceDim = int(globalAttrs['SpaceDim'])

        # Now read all the component names
        self.data = {}
        allNames = attrs.values()
        
        # self.compNames is an ordered list of component names
        self.compNames = []
        for i in range(0, self.num_comps):
            name = attrs['component_' + str(i)]
            #print('component_' + str(i) + ': ' + name)
            # Work out if it's a vector or scalar
            actualName = name
            if name[0] == 'x' or name[0] == 'y' or name[0] == 'z' and actualName in allNames:
                # Vector. Don't add if we've already got it
                actualName = name[1:]
                if actualName not in self.data.keys():
                    self.data[actualName] = {self.NUM_COMPS: 1}
                else:
                    self.data[actualName][self.NUM_COMPS] = self.data[actualName][self.NUM_COMPS] + 1 
            else:
                self.data[name] = {self.NUM_COMPS: 1}
                
                
            self.data[actualName][self.DATA] = [None] * self.num_levels
            self.compNames.append(actualName)
                
        #print(self.data.keys())
        #print(self.compNames)
            

        #print(self.data)

        grp0 = h5File['level_0']
        #print('Group level_0 has keys: ' + str(grp0.keys()))
        #print('  and attributes: ' + str(grp0.attrs.keys()))


        self.levels = [None] * self.num_levels
        for level in range(0, self.num_levels):
            levelGroup = h5File['level_' + str(level)]
            groupAtts = levelGroup.attrs
            boxes = levelGroup['boxes']
            self.levels[level] = {self.DX: groupAtts['dx'], self.DT: groupAtts['dt'],
                                  self.REF_RATIO: groupAtts['ref_ratio'], self.BOXES: list(boxes)}

            # Some attributes on level 0 apply to whole hierarchy
            if level == 0:
                self.time = groupAtts['time']
                #self.isPeriodic = [None] * self.spaceDim
                #for dim_i in range(0, self.spaceDim):
                #    self.isPeriodic[dim_i] = groupAtts['is_periodic_' + str(dim_i)]

                #self.tag_buffer_size = groupAtts['tag_buffer_size'] 
                self.prob_domain = groupAtts['prob_domain']
                #print('Problem domain: ' + str(self.prob_domain))
                self.domainSize = [self.prob_domain[i]*self.levels[level][self.DX] for i in range(0, len(self.prob_domain))]
        
             # Now get the  data for each field, on each level
        
            data = levelGroup['data:datatype=0']
            data_offsets = levelGroup['data:offsets=0']
            data_atts = levelGroup['data_attributes']
           
            advVel = levelGroup['advVel:datatype=0']
            advVel_offsets = levelGroup['advVel:offsets=0']
            advVel_atts = levelGroup['advVel_attributes']
           
            dataUnshaped = data[()]
            # Data is sorted by box then by component
            
            
            #print(data_offsets[()])
            #print(dataUnshaped)
            #print('data length: ' + str(len(dataUnshaped)))
            
            numComps = 0
            for compName in self.data.keys():
                numComps = numComps + self.data[compName][self.NUM_COMPS]
            
            offset = 0
            for box in self.levels[level][self.BOXES]:
                #print(box)
                # Box = [lo_i lo_j hi_i hi_j]
                lo_i = box[0]
                lo_j = box[1]
                hi_i = box[2]
                hi_j = box[3]
    
                numRows = hi_j + 1 - lo_j
                numCols = hi_i + 1 - lo_i
                numCells = numRows*numCols*numComps
                #print(str(numCells))
                dataBox = dataUnshaped[offset:offset + numCells]
                
                offset = offset + numCells
                #print(dataBox)
                
                # Now split into individual components
                
    
                # data contains all fields on this level, sort into individual fields
                
                comp_offset_start = 0;
                #print(self.data.keys())
                for compName in self.compNames:

                    numBoxCells = numRows * numCols
                    #print('Num cells in a box: ' + str(numBoxCells))
                    comp_offset_finish = comp_offset_start + numBoxCells
                    
                    dataBoxComp = dataBox[comp_offset_start:comp_offset_finish]
                    
                    #if compName == 'Porosity':
                    #    print(dataBoxComp)
                    
                    #print(compData)
                    
                    comp_offset_start = comp_offset_finish
                    
                    
                    # Reshape data
                    #shapedData = []
                    # First work out what the shape should be!
                    #print(self.levels[level][BOXES])
                   
                    reshapedData = dataBoxComp.reshape((numRows, numCols))
                    
                    #if compName == 'Enthalpy':
                    cmap2 = mpl.colors.LinearSegmentedColormap.from_list('my_colormap',
                                           ['blue','black','red'],
                                           256)

                    rotatedData = np.flipud(reshapedData)
                    # tell imshow about color map so that only set colors are used
                    #img = pyplot.imshow(rotatedData,interpolation='nearest',
                    #            cmap = cmap2,origin='lower')
            
                    #pyplot.colorbar(img,cmap=cmap2 )
            
                    #pyplot.show() 
        
                    #shapedData.append(reshapedData)
  
                    # Should really get data into a nice format like a np array
                    if self.data[compName][self.DATA][level]:
                        self.data[compName][self.DATA][level].append(reshapedData)
                    else:
                        self.data[compName][self.DATA][level] = [reshapedData]
                    #print('Got data for ' + compName)
      


        
        #print(self.levels[0])
        
           
               
        h5File.close()
        
        
    def channelProperties(self, doPlots = False):
        # Reconstruct level 0 porosity as single np array
        width = self.prob_domain[2] +1 - self.prob_domain[0]
        height = self.prob_domain[3]  +1 - self.prob_domain[1]

        #porosity = np.empty([height, width])

        porosity = self.singleBox('Porosity')
                    
                    
        #print(porosity)

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
                        thisChimneyPositionInRow = chimneyPosRow[chimney_i]
    	                chimneyPositions[chimney_i].append(thisChimneyPositionInRow)
    	    
    	        # Now, add extra rows to chimneyPositions vector if needed
    	        for chimney_i in range(len(chimneyPositions), chimneysInRow):
    	            chimneyPositions.append([chimneyPosRow[chimney_i]])
    
        

        # Get an idea of the channel depths
        chanDepths = []
        for i in range(0, len(chimneyPositions)):
            chanDepths.append(len(chimneyPositions[i]))
            
        
        maxChanDepth = 0
        if chanDepths:
            maxChanDepth = max(chanDepths)
            

        averageChanPositions = []
        relChanPositions = []
        numChannels = 0 
        for i in range(0, len(chimneyPositions)):
            if chanDepths[i] > 0.5*maxChanDepth:
                averageChanPositions.append(np.mean(chimneyPositions[i]))
                relChanPositions.append(averageChanPositions[-1]/width)
                numChannels  = numChannels + 1
            #else:
                #print('Discarding channel as too short: '  + str(chimneyPositions[i]))
        #print('Number of channels: ' + str(numChannels))
        #print('Channel positions: ' + str(averageChanPositions))
        #print('Relative channel positions: ' + str(relChanPositions))     
        
        

        #doPlots = False
        if doPlots:
            print('Making plot')
            self.plotField('Porosity')
            
        channelSpacing = [(relChanPositions[i] - relChanPositions[i-1]) for i in range(1, len(relChanPositions))]
        #channelSpacing = channelSpacing*self.levels[0][DX]
            
        return [numChannels, relChanPositions, channelSpacing]
    
    def plotField(self, field):
        
        fieldArray = self.singleBox(field)
        
        cmap2 = mpl.colors.LinearSegmentedColormap.from_list('my_colormap',
                                           ['blue','black','red'],
                                           256)

        # tell imshow about color map so that only set colors are used
        img = pyplot.imshow(fieldArray,interpolation='nearest',
                    cmap = cmap2,origin='lower')

        # make a color bar
        pyplot.colorbar(img,cmap=cmap2 )

        pyplot.show() 
        
    def singleBox(self, field):
        width = self.prob_domain[2] +1 - self.prob_domain[0]
        height = self.prob_domain[3]  +1 - self.prob_domain[1]

        fieldArr = np.empty([height, width])
        
        for box_i in range(0, len(self.levels[0]['boxes'])):
            box = self.levels[0]['boxes'][box_i]
    
            #print(box)
            thisBox = self.data[field]['data'][0][box_i]
            
            
            #print(thisBox)
    
            ioffset = box[0]-self.prob_domain[0]
            joffset = box[1]-self.prob_domain[1]
    
            for i in range(0, box[2]+1-box[0]):
                for j in range(0, box[3]+1-box[1]):
                    #print(str(i) + ', ' + str(j))
                    val = thisBox[j][i]
                    fieldArr[j+joffset][i+ioffset] = thisBox[j][i]
                    
        return fieldArr


