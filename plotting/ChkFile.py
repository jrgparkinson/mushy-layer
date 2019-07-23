# Class to load a checkpoint file and perform operations on it
# e.g. find where brine channels are
import h5py
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot
import copy
import os
from mushyLayerRunUtils import read_inputs
import sys




class ChkFile:

    NUM_COMPS = 'num_comps'
    DATA = 'data'
    DX = 'dx'
    DT = 'dt'
    REF_RATIO = 'ref_ratio'
    BOXES = 'boxes'

    def __init__(self, filename, inputs_loc=None):
        self.filename = filename

        self.directory = os.path.dirname(filename)
        if not inputs_loc:
            inputs_loc = os.path.join(self.directory, 'inputs')

        if not os.path.exists(inputs_loc):
            print('Warning - ChkFile cannnot find inputs file, so unable to determine parameters')
            # sys.exit(-1)
            self.inputs = None
            return

        self.inputs = read_inputs(inputs_loc)


        h5_file = h5py.File(self.filename, 'r')

        # print('Loaded HDF5 file with groups: ' + str(h5File.keys()))
        # print(' and attributes: ' + str(h5File.attrs.keys()))

        # Get all the different bits of data from the file
        # print(h5File)

        chombo_group = h5_file['Chombo_global']
        global_attrs = chombo_group.attrs

        attrs = h5_file.attrs

        self.time = attrs['time']
        self.iteration = int(attrs['iteration'])
        self.max_level = int(attrs['max_level'])
        self.num_levels = int(attrs['num_levels'])
        self.regrid_interval = int(attrs['regrid_interval_0'])
        self.num_comps = int(attrs['num_components'])
        self.space_dim = int(global_attrs['SpaceDim'])

        # Now read all the component names
        self.data = {}
        all_names = []
        for i in range(0, self.num_comps):
            name = attrs['component_' + str(i)]
            name = name.decode('UTF-8')
            # Work out if it's a vector or scalar
            # if name[0] == 'x' or name[0] == 'y' or name[0] == 'z' and name[1:] in all_names:
            #     # Vector. Don't add if we've already got it
            #     if name[1:] not in self.data.keys():
            #         self.data[name[1:]] = {NUM_COMPS: 1}
            #     else:
            #         self.data[name[1:]][NUM_COMPS] = self.data[name[1:]][NUM_COMPS] + 1
            # else:
            #     self.data[name] = {NUM_COMPS: 1}
            self.data[name] = {self.NUM_COMPS: 1}
            all_names.append(name)

        # print(self.data)

        # grp0 = h5_file['level_0']
        # print('Group level_0 has keys: ' + str(grp0.keys()))
        # print('  and attributes: ' + str(grp0.attrs.keys()))

        self.levels = [None] * self.num_levels
        for level in range(0, self.num_levels):
            level_group = h5_file['level_' + str(level)]
            group_atts = level_group.attrs
            boxes = level_group['boxes']
            self.levels[level] = {self.DX: group_atts['dx'], self.DT: group_atts['dt'],
                                  self.REF_RATIO: group_atts['ref_ratio'], self.BOXES: list(boxes)}

            # Some attributes on level 0 apply to whole hierarchy
            if level == 0:
                self.time = group_atts['time']
                self.is_periodic = [None] * self.space_dim
                for dim_i in range(0, self.space_dim):
                    self.is_periodic[dim_i] = group_atts['is_periodic_' + str(dim_i)]

                self.tag_buffer_size = group_atts['tag_buffer_size']
                self.prob_domain = group_atts['prob_domain']
                # print('Problem domain: ' + str(self.prob_domain))

        # Now get the actual data for each field, on each level
        # print('Fields loaded: ' + str(self.data.keys()))
        for comp_name in self.data.keys():
            self.data[comp_name][self.DATA] = [None] * self.num_levels

            for level in range(0, self.num_levels):
                level_group = h5_file['level_' + str(level)]
                # print(compName + ':offsets=0')
                # offsets = level_group[comp_name + ':offsets=0']

                component = 0

                # Need to get data differently if this is a vector
                is_vector =  (comp_name[0] == 'x' or comp_name[0] == 'y' or comp_name[0] == 'z' and comp_name in all_names)
                if is_vector:
                    if comp_name[0] == 'x':
                        component = 0

                    elif comp_name[0] == 'y':
                        component = 1
                    elif comp_name[0] == 'z':
                        component = 2

                    # Hardwired to 2D for now
                    num_comps = self.space_dim

                    data = level_group[comp_name[1:] + ':datatype=0']
                else:
                    data = level_group[comp_name + ':datatype=0']
                    component = 0

                    num_comps = 1

                # data_atts = level_group[comp_name + '_attributes']

                # Reshape data
                data_unshaped = data[()]

                shaped_data = []
                # First work out what the shape should be!
                # print(self.levels[level][BOXES])

                offset = 0

                for box in self.levels[level][self.BOXES]:
                    # print(box)
                    # Box = [lo_i lo_j hi_i hi_j]
                    lo_i = box[0]
                    lo_j = box[1]
                    hi_i = box[2]
                    hi_j = box[3]

                    num_rows = hi_j + 1 - lo_j
                    num_cols = hi_i + 1 - lo_i
                    num_comp_cells = num_rows * num_cols
                    num_cells = num_comp_cells * num_comps
                    # print(str(numCells))
                    data_box = data_unshaped[offset:offset + num_cells]

                    comp_offset = num_comp_cells*component

                    data_comp_box = data_box[comp_offset:comp_offset+num_comp_cells]

                    reshaped_data = data_comp_box.reshape((num_rows, num_cols))
                    # reshaped_data = data_comp_box.reshape((num_cols, num_rows))
                    reshaped_data = np.array(reshaped_data).transpose()
                    shaped_data.append(reshaped_data)

                    offset = offset + num_cells

                # Should really get data into a nice format like a np array
                self.data[comp_name][self.DATA][level] = shaped_data

        # print(self.levels[0])

        self.compute_diagnostic_vars()

        h5_file.close()

    def compute_diagnostic_vars(self):
        # First get parameters

        if self.inputs is None:
            print('No inputs file available, so do not know parameters and cannot compute diagnostic variables')
            sys.exit(-1)

        conc_ratio = float(self.inputs['parameters.compositionRatio'])
        stefan = float(self.inputs['parameters.stefan'])
        cp = float(self.inputs['parameters.specificHeatRatio'])
        pc = float(self.inputs['parameters.waterDistributionCoeff'])
        theta_eutectic = 0.0

        enthalpy = self.get_data('Enthalpy')
        bulk_salinity = self.get_data('Bulk concentration')

        enthalpy_eutectic = np.empty(enthalpy.shape)
        enthalpy_solidus = np.empty(enthalpy.shape)
        enthalpy_liquidus = np.empty(enthalpy.shape)

        porosity = np.empty(enthalpy.shape)
        temperature = np.empty(enthalpy.shape)
        liquid_salinity = np.empty(enthalpy.shape)
        solid_salinity = np.empty(enthalpy.shape)

        # Now compute bounding energies
        cols = enthalpy.shape[0]
        rows = enthalpy.shape[1]
        for j in range(cols):

            for i in range(rows):
                eutectic_porosity = (conc_ratio + bulk_salinity[j, i]) / (theta_eutectic + conc_ratio)
                enthalpy_eutectic[j, i] = eutectic_porosity * (stefan + theta_eutectic * (1 - cp)) + cp * theta_eutectic
                enthalpy_solidus[j, i] = cp * (theta_eutectic + max(0.0, (-bulk_salinity[j, i] - conc_ratio) / pc))
                enthalpy_liquidus[j, i] = stefan - bulk_salinity[j, i] + theta_eutectic + theta_eutectic

        # Compute diagnostic variables
        for j in range(cols):

            for i in range(rows):
                if enthalpy[j, i] <= enthalpy_solidus[j, i]:
                    porosity[j, i] = 0.0
                    temperature[j, i] = enthalpy[j, i] / cp
                    liquid_salinity[j, i] = 0.0
                    solid_salinity[j, i] = bulk_salinity[j, i]
                elif enthalpy[j, i] <= enthalpy_eutectic[j, i]:
                    porosity[j, i] = (enthalpy[j, i] - theta_eutectic * cp) / (stefan + theta_eutectic * (1 - cp))
                    temperature[j, i] = theta_eutectic
                    liquid_salinity[j, i] = theta_eutectic
                    solid_salinity[j, i] = bulk_salinity[j, i] / (1 - porosity[j, i])
                elif enthalpy[j, i] < enthalpy_liquidus[j, i]:
                    a = conc_ratio * (cp - 1) + stefan * (pc - 1)
                    b = conc_ratio * (1 - 2 * cp) + enthalpy[j, i] * (1 - pc) \
                        - bulk_salinity[j, i] * (cp - 1) - pc * stefan
                    c = (bulk_salinity[j, i] + conc_ratio) * cp + pc * enthalpy[j, i]
                    porosity[j, i] = (-b - np.sqrt(b * b - 4 * a * c)) / (2 * a)

                    liquid_salinity[j, i] = (bulk_salinity[j, i] + conc_ratio * (1 - porosity[j, i])) / (
                            porosity[j, i] + pc * (1 - porosity[j, i]))
                    temperature[j, i] = -liquid_salinity[j, i]
                    solid_salinity[j, i] = (pc * bulk_salinity[j, i] - conc_ratio * porosity[j, i]) / (
                            porosity[j, i] + pc * (1 - porosity[j, i]))
                else:
                    porosity[j, i] = 1.0
                    temperature[j, i] = enthalpy[j, i] - stefan
                    liquid_salinity[j, i] = bulk_salinity[j, i]
                    solid_salinity[j, i] = 0.0

        # Finally, save data into boxes
        self.save_data('Porosity', porosity)
        self.save_data('Temperature', temperature)
        self.save_data('LiquidSalinity', liquid_salinity)
        self.save_data('SolidSalinity', solid_salinity)

    def save_data(self, var_name, lev0dat):
        # width = self.prob_domain[2] + 1 - self.prob_domain[0]
        # height = self.prob_domain[3] + 1 - self.prob_domain[1]

        # If data doesn't exist, make it
        if var_name not in self.data:
            enthalpy = self.data['Enthalpy']
            deep_copy = copy.deepcopy(enthalpy)
            self.data[var_name] = deep_copy

        for box_i in range(0, len(self.levels[0]['boxes'])):
            box = self.levels[0]['boxes'][box_i]


            level_0_data_box = self.data[var_name]['data'][0][box_i]

            # ND version:
            box_size = [box[self.space_dim + i] + 1 - box[i] for i in range(0, self.space_dim)]
            offsets = tuple([box[i] - self.prob_domain[i] for i in range(0, self.space_dim)])

            loopover = [range(s) for s in box_size]
            import itertools

            prod = itertools.product(*loopover)

            for idx in prod:
                idx_plus_offset = tuple(map(lambda x, y: x + y, idx, offsets))
                level_0_data_box[idx] = lev0dat[idx_plus_offset]

            # ioffset = box[0] - self.prob_domain[0]
            # joffset = box[1] - self.prob_domain[1]
            #
            # for i in range(0, box[2] + 1 - box[0]):
            #     for j in range(0, box[3] + 1 - box[1]):
            #         level_0_data_box[j][i] = lev0dat[j + joffset][i + ioffset]

    def get_data(self, var_name, rotate_dims = False):
        # Reconstruct level 0 data as single np array

        # N dimensional version:
        domain_size = [self.prob_domain[self.space_dim + i] + 1 - self.prob_domain[i] for i in range(0, self.space_dim)]
        lev0_dat = np.empty(domain_size)
        for box_i in range(0, len(self.levels[0]['boxes'])):
            box = self.levels[0]['boxes'][box_i]
            lev0_dat_box = self.data[var_name]['data'][0][box_i]

            box_size = [box[self.space_dim + i] + 1 - box[i] for i in range(0, self.space_dim)]

            offsets = tuple([box[i] - self.prob_domain[i] for i in range(self.space_dim)])
            # offsets = tuple([box[i] - self.prob_domain[i] for i in range(self.space_dim-1, -1, -1)])

            loopover = [range(s) for s in box_size]
            import itertools

            prod = itertools.product(*loopover)

            for idx in prod:
                idx_plus_offset = tuple(map(lambda x, y: x + y, idx, offsets))
                lev0_dat[idx_plus_offset] = lev0_dat_box[idx]

        if rotate_dims:
            lev0_dat = np.array(lev0_dat).transpose()


        # width = self.prob_domain[2] + 1 - self.prob_domain[0]
        # height = self.prob_domain[3] + 1 - self.prob_domain[1]
        #
        # lev0_dat = np.empty([height, width])
        #
        # for box_i in range(0, len(self.levels[0]['boxes'])):
        #     box = self.levels[0]['boxes'][box_i]
        #
        #     # print(box)
        #     lev0_dat_box = self.data[var_name]['data'][0][box_i]
        #
        #     ioffset = box[0] - self.prob_domain[0]
        #     joffset = box[1] - self.prob_domain[1]
        #
        #     for i in range(0, box[2] + 1 - box[0]):
        #         for j in range(0, box[3] + 1 - box[1]):
        #             lev0_dat[j + joffset][i + ioffset] = lev0_dat_box[j][i]




        return lev0_dat

    def channel_properties(self, do_plots=False):
        porosity = self.get_data('Porosity')

        # width = self.prob_domain[2] +1 - self.prob_domain[0]
        width = porosity.shape[1]
        # print('Domain width ' + str(width))

        # Iterate over porosity field
        cols = porosity.shape[0]
        rows = porosity.shape[1]
        channels = [None] * cols
        chimney_positions = []
        for j in range(cols):
            # chimneys in row
            chimneys_in_row = 0
            average_porosity = 0

            currently_liquid = False
            chimney_pos_row = []

            for i in range(rows):
                # print porosity[i,j]
                chi = porosity[j, i]
                average_porosity = average_porosity + chi
                if chi < 1.0 and currently_liquid:
                    # Just left liquid region
                    currently_liquid = False

                elif chi > 0.999 and not currently_liquid:
                    # Entered a liquid region
                    chimneys_in_row = chimneys_in_row + 1
                    chimney_pos_row.append(i)
                    currently_liquid = True

            average_porosity = average_porosity / rows
            if average_porosity > 0.99:
                channels[j] = 0  # This region was entirely liquid, can't have a channel
            else:
                channels[j] = chimneys_in_row

                # First add chimney positions for chimneys we've already found
                # print(str(chimneyPositions))
                # print(str(len(chimneyPositions)))

                if chimneys_in_row == 0:
                    continue

                for chimney_i in range(0, len(chimney_positions)):
                    if chimney_i < len(chimney_pos_row):
                        chimney_positions[chimney_i].append(chimney_pos_row[chimney_i])

                # Now, add extra rows to chimneyPositions vector if needed
                for chimney_i in range(len(chimney_positions), chimneys_in_row):
                    if chimney_i < len(chimney_pos_row):
                        chimney_positions.append([chimney_pos_row[chimney_i]])

        # print(str(chimney_positions))
        # Filter out short channels
        chimney_lengths = []
        for chimney_i in range(0, len(chimney_positions)):
            chimney_lengths.append(float(len(chimney_positions[chimney_i])))

        av_chim_length = np.mean(chimney_lengths)
        long_enough_chimneys = []
        for chimney_i in range(0, len(chimney_positions)):
            if len(chimney_positions[chimney_i]) > av_chim_length / 2:
                long_enough_chimneys.append(chimney_positions[chimney_i])

        # print(str(long_enough_chimneys))

        average_chan_positions = []
        rel_chan_positions = []
        for i in range(0, len(long_enough_chimneys)):
            average_chan_positions.append(np.mean(long_enough_chimneys[i]))
            rel_chan_positions.append(average_chan_positions[-1] / width)

        num_channels = len(long_enough_chimneys)
        # print('Number of channels: ' + str(numChannels))
        # print('Channel positions: ' + str(averageChanPositions))
        # print('Relative channel positions: ' + str(relChanPositions))

        # doPlots = False
        if do_plots:
            cmap2 = mpl.colors.LinearSegmentedColormap.from_list('my_colormap',
                                                                 ['blue', 'black', 'red'],
                                                                 256)

            # tell imshow about color map so that only set colors are used
            img = pyplot.imshow(porosity, interpolation='nearest',
                                cmap=cmap2, origin='lower')

            # make a color bar
            pyplot.colorbar(img, cmap=cmap2)

            pyplot.show()

        return [num_channels, rel_chan_positions]


    def get_mesh_grid(self, level=0, rotate_dims=False):

        dx = self.levels[level][self.DX]
        dy = dx

        components = list(self.data.keys())

        field_array = self.get_data(components[0])
        grid_size = field_array.shape

        # Make sure the grid stretches from the start of the first cell to the end of the last cell
        # this means we stretch dx slightly out of proportion, but it ensures plot limits are correct
        x_max = (grid_size[0] + 1) * dx
        y_max = (grid_size[1] + 1) * dy

        # grid_dx = x_max / grid_size[0]
        # grid_dy = y_max / grid_size[1]
        #
        # y, x = np.mgrid[slice(0, x_max, grid_dx),
        #                 slice(0, y_max, grid_dy)]

        y, x = np.mgrid[slice(dx/2, x_max-dx/2, dx),
                        slice(dy/2, y_max-dy/2, dy)]

        if rotate_dims:
            y_new = x.transpose()
            x_new = y.transpose()

            x = x_new
            y = y_new




        return x, y


    def get_permeability(self, permeability_function='kozeny', rotate_dims = False):

        porosity = self.get_data('Porosity', rotate_dims)

        if permeability_function == 'kozeny':

            # Cap max porosity just below one to avoid dividing by 0
            porosity = np.clip(porosity, 0, 1-10**(-10))

            liquid_permeability = porosity**3 / (1-porosity)**2

            hele_shaw_permeability = 1/float(self.inputs['parameters.nonDimReluctance'])

            total_permeability = (hele_shaw_permeability**(-1) + liquid_permeability**(-1))**(-1)

        elif permeability_function == 'cubic':
            total_permeability = porosity**3

        else:
            total_permeability = 1.0


        return total_permeability