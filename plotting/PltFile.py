import h5py
import numpy as np
from matplotlib import pyplot
import re
import os
# import yt
import xarray as xr
from shapely.ops import cascaded_union
from shapely.geometry import Polygon
import geopandas as gpd
from mushyLayerRunUtils import read_inputs

def compute_z(porosity_slice, y_slice, porosity):

    # valid_cells = np.argwhere(porosity_slice < porosity)
    #
    # j_ml_above = min(valid_cells)
    #
    z = np.interp(porosity, porosity_slice, y_slice)

    return z


class PltFile:
    """ Class to load a plot file and perform operations on it
         e.g. find where brine channels are"""
    NUM_COMPS = 'num_comps'
    DATA = 'data'
    DX = 'dx'
    DT = 'dt'
    REF_RATIO = 'ref_ratio'
    BOXES = 'boxes'

    YT = 'YT'
    NATIVE = 'native'
    XARRAY = 'xarray'

    indices = None
    reflect = None
    ds_amr = None

    FIELD_LABEL = {'Porosity': '$\chi$',
                   'Bulk concentration': '$\Theta$',
                   'Liquid concentration': '$\Theta_l$',
                   'Permeability': '$\Pi$',
                   'xAdvection velocity': '$u$',
                   'yAdvection velocity': '$w$',
                   'streamfunction': '$\psi$'}

    def __init__(self, filename, load_data=False):
        if not os.path.exists(filename):
            print('PltFile: file does not exist %s ' % filename)
            return

        self.filename = filename

        # Get the plot prefix and frame number assuming a format
        prefix_format = '(.*)-(\d+)\.2d\.hdf5'
        m = re.search(prefix_format, self.filename)

        if m and m.groups() and len(m.groups()) == 2:
            self.plot_prefix = m.groups(0)
            self.frame = m.groups(1)

        else:
            self.plot_prefix = None
            self.frame = -1

        self.data_loaded = False
        self.ds = None

        self.data_load_method = self.XARRAY

        if load_data:
            self.load_data()

        # Initialise to bogus values
        self.iteration = -1
        self.max_level = -1
        self.num_levels = -1
        self.num_comps = -1
        self.space_dim = -1
        self.comp_names = []
        self.levels = []
        self.time = -1
        self.prob_domain = None
        self.domain_size = []
        self.xarr_data = None
        self.level_outlines = []


        # Now read all the component names
        self.data = {}


        # get inputs
        output_folder = os.path.abspath(os.path.join(self.filename, '..'))
        inputs_file_loc = os.path.join(output_folder, 'inputs')
        if os.path.exists(inputs_file_loc):
            self.inputs = read_inputs(inputs_file_loc)
        else:
            self.inputs = None

    def load_data(self, zero_x=False):
        if self.data_loaded:
            return

        if self.data_load_method in (self.NATIVE, self.XARRAY):
            self.load_data_native(zero_x)
        # else:
        #     self.load_data_yt()

        self.data_loaded = True

    # def load_data_yt(self):
    #     self.ds = yt.load(self.filename)

    def unload_data(self):
        del self.ds_levels

        self.data_loaded = False



    # noinspection PyUnresolvedReferences
    def load_data_native(self, zero_x=False):

        print('Loading %s' % self.filename)

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
        # self.regrid_interval = int(attrs['regrid_interval_0'])
        self.num_comps = int(attrs['num_components'])
        self.space_dim = int(global_attrs['SpaceDim'])

        # Now read all the component names
        self.data = {}

        for i in range(0, self.num_comps):
            name = attrs['component_' + str(i)]
            name = name.decode('UTF-8')

            # Previously, we treated vectors and scalars differently,
            # now we just store vector components as separate scalar fields
            # retained previous code (commented out below) in case I ever want it
            actual_name = name

            self.data[name] = {self.NUM_COMPS: 1}

            self.data[actual_name][self.DATA] = [None] * self.num_levels
            self.comp_names.append(actual_name)



        ds_levels = []

        self.levels = [None] * self.num_levels
        for level in range(0, self.num_levels):
            level_group = h5_file['level_' + str(level)]
            group_atts = level_group.attrs
            boxes = level_group['boxes']
            lev_dx =  group_atts['dx']
            self.levels[level] = {self.DX: lev_dx, self.DT: group_atts['dt'],
                                  self.REF_RATIO: group_atts['ref_ratio'], self.BOXES: list(boxes)}

            # Some attributes on level 0 apply to whole hierarchy
            if level == 0:
                self.time = group_atts['time']

                self.prob_domain = group_atts['prob_domain']

                self.domain_size = [self.prob_domain[i] * self.levels[level][self.DX] for i in
                                    range(0, len(self.prob_domain))]

                # Moving to ND
                self.full_domain_size = self.domain_size
                for i in range(self.space_dim, self.space_dim + self.space_dim):
                    self.full_domain_size[i] = self.full_domain_size[i] + lev_dx
                    # self.fullDomainSize[3] = self.fullDomainSize[3] + lev_dx

            # Now get the  data for each field, on each level
            data = level_group['data:datatype=0']

            # Some other stuff we can get, but don't at the moment:
            # data_offsets = level_group['data:offsets=0']
            # data_atts = level_group['data_attributes']
            # advVel = level_group['advVel:datatype=0']
            # advVel_offsets = level_group['advVel:offsets=0']
            # advVel_atts = level_group['advVel_attributes']

            data_unshaped = data[()]
            # Data is sorted by box then by component

#            print(data_offsets[()])
            # print(dataUnshaped)
            # print('data length: ' + str(len(dataUnshaped)))

            num_comps = 0
            for comp_name in self.data.keys():
                num_comps = num_comps + self.data[comp_name][self.NUM_COMPS]

            offset = 0

            ds_boxes = []

            # Initialise with a box spanning the whole domain, then add data where it exists
            # Important to do it like this for refined levels, where the whole domain isn't covered with data
            # lev_dom_box_x = np.arange(self.fullDomainSize[0] + lev_dx / 2, self.fullDomainSize[2] - lev_dx / 2,
            #                           lev_dx)
            # lev_dom_box_y = np.arange(self.fullDomainSize[1] + lev_dx / 2, self.fullDomainSize[3] - lev_dx / 2,
            #                           lev_dx)
            # blank_data = np.empty((lev_dom_box_x.size, lev_dom_box_y.size))

            size = []
            for i in range(self.space_dim):
                lev_dom_box_dir = np.arange(self.full_domain_size[i] + lev_dx / 2, self.full_domain_size[i + self.space_dim] - lev_dx / 2,
                                            lev_dx)
                size.append(lev_dom_box_dir.size)


            blank_data = np.empty(tuple(size))

            blank_data[:] = np.nan

            # Use indexes rather than x, y for now - then convert to x,y later
            # this is to avoid issues with floating point arithmetic when merging datasets
            # (we can end up trying to merge datasets where x coordinates differ by ~ 10^{-10}, creating nonsense)


            index_coords_names = ['i', 'j', 'k', 'l', 'm'] # add more here if more dimensions
            coords = {}
            box_size = ()
            for d in range(self.space_dim):
                coords_dir = np.arange(self.prob_domain[d], self.prob_domain[self.space_dim + d] + 1)
                coords[index_coords_names[d]] = coords_dir
                box_size = box_size + (coords_dir.size, )  # append to tuple of sizes

            # i = np.arange(self.prob_domain[0], self.prob_domain[self.spaceDim]   + 1)
            # j = np.arange(self.prob_domain[1], self.prob_domain[self.spaceDim+1] + 1)

            # blank_data = np.empty((j.size, i.size))
            blank_data = np.empty(box_size)

            # ds_dom_box = xr.Dataset({},
            #                     coords={'x': lev_dom_box_x, 'y': lev_dom_box_y})

            ds_dom_box = xr.Dataset({}, coords=coords)
                                # coords={'i': i, 'j': j})


            for comp_name in self.comp_names:
                s = blank_data.shape
                if not s[1] == len(coords['i']):
                    blank_data = blank_data.T

                # ds_dom_box[comp_name] = xr.DataArray(blank_data, dims=['y', 'x'],
                #                           coords={'x': lev_dom_box_x, 'y': lev_dom_box_y, 'level': level})

                extended_coords = coords
                extended_coords['level'] = level
                dims = index_coords_names[:self.space_dim]
                dims= dims[::-1] # reverse list so we have k, j, i etc
                ds_dom_box[comp_name] = xr.DataArray(blank_data, dims=dims, # dims=['j', 'i'],
                                                     coords=extended_coords)

            ds_boxes.append(ds_dom_box)

            polygons = []

            for box in self.levels[level][self.BOXES]:
                # print(box)
                # Box = [lo_i lo_j hi_i hi_j]
                lo_i = box[0]
                lo_j = box[1]
                hi_i = box[self.space_dim]
                hi_j = box[self.space_dim + 1]



                lo_indices = [box[i] for i in range(self.space_dim)]
                hi_indices = [box[i] for i in range(self.space_dim, 2 * self.space_dim)]

                # 0.5 because cell centred
                x_box = lev_dx * (0.5 + np.arange(lo_i, hi_i + 1) )
                y_box = lev_dx * (0.5 + np.arange(lo_j, hi_j + 1) )
                if self.space_dim > 2:
                    lo_k = box[2]
                    hi_k = box[self.space_dim + 2]
                    z_box = lev_dx * (0.5 + np.arange(lo_k, hi_k + 1) )

                lo_vals =  [lev_dx *(0.5 + i) for i in lo_indices]
                hi_vals =  [lev_dx *(0.5 + i) for i in hi_indices]

                end_points = [[lo_vals[i]- lev_dx/2, hi_vals[i]+lev_dx/2] for i in range(self.space_dim)]

                polygon_vertices_2d = [(x_box[0] - lev_dx / 2, y_box[0] - lev_dx / 2),
                     (x_box[-1] + lev_dx / 2, y_box[0] - lev_dx / 2),
                     (x_box[-1] + lev_dx / 2, y_box[-1] + lev_dx / 2),
                     (x_box[0] - lev_dx / 2, y_box[-1] + lev_dx / 2)]

                polygon_vertices_2d_new = [(lo_vals[0] - lev_dx / 2, lo_vals[1] - lev_dx / 2),
                                       (hi_vals[0] + lev_dx / 2, lo_vals[1] - lev_dx / 2),
                                       (hi_vals[0] + lev_dx / 2, hi_vals[1] + lev_dx / 2),
                                       (lo_vals[0] - lev_dx / 2, hi_vals[1] + lev_dx / 2)]



                from itertools import product, permutations, chain
                # Construct vertices in n dimensions
                polygon_vertices_auto = list(product(*end_points))

                polygon_vertices_auto = sorted(polygon_vertices_auto, key=lambda x: np.arctan(x[1]/max(abs(x[0]), 0.0001)))



                # For plotting level outlines
                # polygons.append(Polygon(
                #     [(x_box[0]-lev_dx/2, y_box[0]-lev_dx/2),
                #      (x_box[-1]+lev_dx/2, y_box[0]-lev_dx/2),
                #      (x_box[-1]+lev_dx/2, y_box[-1]+lev_dx/2),
                #      (x_box[0]-lev_dx/2, y_box[-1]+lev_dx/2)]))

                # Construct vertices in n dimensions
                # polygon_vertices = []

                # polygons.append(Polygon(
                #     [(x_box[0] - lev_dx / 2, y_box[0] - lev_dx / 2),
                #      (x_box[-1] + lev_dx / 2, y_box[0] - lev_dx / 2),
                #      (x_box[-1] + lev_dx / 2, y_box[-1] + lev_dx / 2),
                #      (x_box[0] - lev_dx / 2, y_box[-1] + lev_dx / 2)]))

                poly = Polygon(polygon_vertices_auto)
                if poly.is_valid:
                    polygons.append(poly)


                # i_box = np.arange(lo_i, hi_i + 1)
                # j_box = np.arange(lo_j, hi_j + 1)

                n_cells_dir = [ hi_indices[d]+1 - lo_indices[d] for d in range(self.space_dim)]
                #  n_cells_dir = [hi_indices[d] + 1 - lo_indices[d] for d in range(self.space_dim-1, -1, -1)]

                # num_rows = hi_j + 1 - lo_j
                # num_cols = hi_i + 1 - lo_i
                # num_cells = num_rows * num_cols * num_comps
                num_box_cells =  np.prod(n_cells_dir)
                num_cells = num_box_cells * num_comps
                # print(str(num_cells))
                data_box = data_unshaped[offset:offset + num_cells]

                offset = offset + num_cells
                # print(data_box)

                # Now split into individual components
                # data contains all fields on this level, sort into individual fields

                comp_offset_start = 0
                # print(self.data.keys())

                coords = {}
                # box_size = ()
                for d in range(self.space_dim):
                    coords_dir = np.arange(lo_indices[d], hi_indices[d] + 1)
                    coords[index_coords_names[d]] = coords_dir

                    # box_size = box_size + (coords_dir.size,)  # append to tuple of sizes

                ds_box = xr.Dataset({}, coords=coords)
                                    # coords={'i': i_box, 'j': j_box})


                for comp_name in self.comp_names:


                    # print('Num cells in a box: ' + str(num_box_cells))
                    comp_offset_finish = comp_offset_start + num_box_cells

                    data_box_comp = data_box[comp_offset_start:comp_offset_finish]

                    comp_offset_start = comp_offset_finish


                    # reshaped_data =    data_box_comp.reshape((num_rows, num_cols))

                    reshaped_data = data_box_comp.reshape(tuple(n_cells_dir))

                    # I think we only need to transpose in 3D
                    if self.space_dim == 3 or self.space_dim == 2:
                        reshaped_data = reshaped_data.transpose()



                    # Should really get data into a nice format like a np array
                    if self.data[comp_name][self.DATA][level]:
                        self.data[comp_name][self.DATA][level].append(reshaped_data)
                    else:
                        self.data[comp_name][self.DATA][level] = [reshaped_data]

                    reshaped_data = np.array(reshaped_data)

                    # Check if scalar or vector
                    trimmed_comp_names = [n[1:] for n in self.comp_names]
                    field_type = 'scalar'
                    if comp_name[0] in ('x', 'y', 'z') and comp_name[1:] in trimmed_comp_names:
                        field_type = 'vector'

                    dim_list = index_coords_names[:self.space_dim]
                    # dim_list = dim_list[::-1]
                    extended_coords = coords
                    extended_coords['level'] = level

                    # print(comp_name)

                    xarr_component_box = xr.DataArray(reshaped_data, dims = dim_list, # ['j', 'i'],
                                                      coords=extended_coords, # {'i': i_box, 'j': j_box, 'level': level},
                                                      attrs={'field_type': field_type})

                    ds_box[comp_name] = xarr_component_box

                ds_boxes.append(ds_box)


            level_outline = gpd.GeoSeries(cascaded_union(polygons))
            self.level_outlines.append(level_outline)

            # ds_level = xr.merge(ds_boxes)

            # Update will replace in place
            first_box = 1
            ds_level = ds_boxes[first_box]
            for b in ds_boxes[first_box+1:]:

                ds_level =  ds_level.combine_first(b)



            # Create x,y,z, coordinates
            x_y_coords_names = ['x', 'y', 'z']

            for d in range(self.space_dim):
                ds_level.coords[x_y_coords_names[d]] = ds_level.coords[index_coords_names[d]] * lev_dx
            # ds_level.coords['x'] = ds_level.coords['i'] * lev_dx
            # ds_level.coords['y'] = ds_level.coords['j'] * lev_dx

            if zero_x:
                ds_level.coords['x'] = ds_level.coords['x'] - min(ds_level.coords['x'])

            # Swap i,j,k to x,y,z coordinates
            for d in range(self.space_dim):
                ds_level = ds_level.swap_dims({index_coords_names[d]: x_y_coords_names[d]})

            # ds_level = ds_level.swap_dims({'i': 'x'})
            #ds_level = ds_level.swap_dims({'j': 'y'})


            #TODO: should level be an attribute or co-ordinate? need to try with actual AMR data
            ds_level.attrs['level'] =  level


            ds_levels.append(ds_level)


        self.ds_levels = ds_levels

        h5_file.close()

    def plot_outlines(self, ax, colors=None):
        """ Plot all level outlines (except level 0)"""

        for level in np.arange(1, len(self.level_outlines)):
            self.plot_outline(ax, level, colors)

    def plot_outline(self, ax, level, colors=None):
        """ Plot level outline for a particular color"""

        # Default colors
        if not colors:
            colors = [[0, 0, 0, 1.0],
                      [1, 0, 0, 1.0],
                      [0, 1, 0, 1.0],
                      [0, 0, 1, 1.0]]

        ec = colors[level][:]

        outline = self.level_outlines[level]

        # Shrink outline slightly

        # outline
        #outline = outline.scale(0.99, 0.99)
        # dx = self.levels[level][self.DX]
        #domain = Polygon([(dx, dx), (1-dx, dx), (1-dx, 1-dx), (dx, 1-dx)])
        #intersect = gpd.sjoin(domain, outline, how="inner", op='intersection')
        #intersect.plot(ax=ax, edgecolor=ec, facecolor=[1,1,1,0], linewidth=2.0)


        outline = outline.scale(0.99,0.99)
        outline.plot(ax=ax, edgecolor=ec, facecolor=[1, 1, 1, 0], linewidth=3.0)

    def channel_properties(self, do_plots=False):
        # Reconstruct level 0 porosity as single np array
        width = self.prob_domain[2] + 1 - self.prob_domain[0]
        # height = self.prob_domain[3] + 1 - self.prob_domain[1]

        # porosity = np.empty([height, width])

        porosity = self.single_box('Porosity')

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
                # print(str(chimney_positions))
                # print(str(len(chimney_positions)))

                if chimneys_in_row == 0:
                    continue

                for chimney_i in range(0, len(chimney_positions)):
                    if chimney_i < len(chimney_pos_row):
                        this_chimney_position_in_row = chimney_pos_row[chimney_i]
                        chimney_positions[chimney_i].append(this_chimney_position_in_row)

                # Now, add extra rows to chimney_positions vector if needed
                for chimney_i in range(len(chimney_positions), chimneys_in_row):
                    chimney_positions.append([chimney_pos_row[chimney_i]])

        # Get an idea of the channel depths
        chan_depths = []
        for i in range(0, len(chimney_positions)):
            chan_depths.append(len(chimney_positions[i]))

        max_chan_depth = 0
        if chan_depths:
            max_chan_depth = max(chan_depths)

        average_chan_positions = []
        rel_chan_positions = []
        num_channels = 0
        for i in range(0, len(chimney_positions)):
            if chan_depths[i] > 0.5 * max_chan_depth:
                average_chan_positions.append(np.mean(chimney_positions[i]))
                rel_chan_positions.append(average_chan_positions[-1] / width)
                num_channels = num_channels + 1
            # else:
            # print('Discarding channel as too short: '  + str(chimney_positions[i]))
        # print('Number of channels: ' + str(num_channels))
        # print('Channel positions: ' + str(average_chan_positions))
        # print('Relative channel positions: ' + str(rel_chan_positions))

        # doPlots = False
        if do_plots:
            print('Making plot')
            self.plot_field('Porosity')

        channel_spacing = [(rel_chan_positions[i] - rel_chan_positions[i - 1]) for i in
                           range(1, len(rel_chan_positions))]
        # channel_spacing = channel_spacing*self.levels[0][DX]

        return [num_channels, rel_chan_positions, channel_spacing]


    def get_mesh_grid_n(self, arr, grow=0):
        x = np.array(arr.coords['x'])
        y = np.array(arr.coords['y'])

        x_min = x[0]
        x_max = x[-1]
        y_min = y[0]
        y_max = y[-1]

        nx = len(x) + grow
        ny = len(y) + grow

        x = np.linspace(x_min, x_max, nx)
        y = np.linspace(y_min, y_max, ny)


        #
        # y, x = np.mgrid[slice(float(y[0]), float(y[-1]) + dy, dy),
        #                 slice(float(x[0]), float(x[-1]) + dx, dx)                     ]
        #
        # y = self.scale_slice_transform(y)
        # x = self.scale_slice_transform(x, no_reflect=True)

        return x, y


    def get_mesh_grid_xarray(self, arr, grow=False):
        x = np.array(arr.coords['x'])
        y = np.array(arr.coords['y'])

        dx = float(x[1]-x[0])
        dy = float(y[1]-y[0])

        if grow:
            x = np.append(x, [float(x[-1])+dx])
            y = np.append(y, [float(y[-1])+dy])

            x = x-dx/2
            y = y-dx/2


        #
        # y, x = np.mgrid[slice(float(y[0]), float(y[-1]) + dy, dy),
        #                 slice(float(x[0]), float(x[-1]) + dx, dx)                     ]
        #
        # y = self.scale_slice_transform(y)
        # x = self.scale_slice_transform(x, no_reflect=True)

        return x, y


    def get_mesh_grid(self, level=0):

        dx = self.levels[level][self.DX]
        dy = dx

        components = list(self.data.keys())
        # components = [c.decode('UTF-8') for c in components]

        field_array = self.single_box(components[0])
        grid_size = field_array.shape

        # x_max = (gridSize[0]+0.5) * dx
        # y_max = (gridSize[1]+0.5) * dy

        # Y, X = np.mgrid[slice(dx/2, x_max, dx),
        #                slice(dy/2, y_max, dy)]

        # Make sure the grid stretches from the start of the first cell to the end of the last cell
        # this means we stretch dx slightly out of proportion, but it ensures plot limits are correct
        x_max = (grid_size[0] + 1) * dx
        y_max = (grid_size[1] + 1) * dy

        grid_dx = x_max / grid_size[0]
        grid_dy = y_max / grid_size[1]

        y, x = np.mgrid[slice(0, x_max, grid_dx),
                        slice(0, y_max, grid_dy)]

        y = self.scale_slice_transform(y)
        x = self.scale_slice_transform(x, no_reflect=True)

        return x, y

    def get_level_data(self, field, level=0, valid_only=False):

        if self.data_load_method == self.YT:
            pass
            #TODO: write this

        if field not in self.data.keys():
            print('Field: %s not found. The following fields do exist: ' % field)
            print(self.data.keys())
            return

        if self.data_load_method == self.XARRAY:

            # ds_lev = self.ds_amr.sel(level=level)
            #ds_lev = self.ds_levels[level].sel(level=level)
            ds_lev = self.ds_levels[level]

            ld = ds_lev[field]

            # Set covered cells to NaN
            # This is really slow, I'm sure there's a faster way
            if valid_only and level < self.num_levels - 1:


                coarseness = self.levels[level][self.REF_RATIO]
                fine_level = np.array(self.ds_levels[level+1][field])
                temp = fine_level.reshape((fine_level.shape[0]// coarseness, coarseness,
                                           fine_level.shape[1] // coarseness, coarseness))
                coarse_fine = np.sum(temp, axis=(1,3))

                isnan = np.isnan(coarse_fine)
                ld = ld.where(isnan == True)

        else:
            ld = self.single_box(field, level)

        ld = self.scale_slice_transform(ld)

        # If this is the x-component of a vector, and we're reflecting, we need to also make the field negative
        if self.reflect:
            if field[0] == 'x':
                ld = -ld
            elif field == 'streamfunction':
                ld = -ld

        return ld

    def scale_slice_transform(self, data, no_reflect=False):
        if self.indices:
            data = data[self.indices]

        if self.reflect and not no_reflect:
            data = np.flip(data, 1)

        return data

    def plot_field(self, field):

        self.load_data()

        if self.data_load_method == self.NATIVE:

            cmap = pyplot.get_cmap('PiYG')

            ld = self.get_level_data(field, 0)
            x, y = self.get_mesh_grid()
            img = pyplot.pcolormesh(x, y, ld, cmap=cmap)

            # make a color bar
            pyplot.colorbar(img, cmap=cmap, label='')

            pyplot.xlabel('$x$')
            pyplot.xlabel('$z$')

        else:
            pass
            # print(self.ds.field_list)

            # plot = yt.SlicePlot(self.ds, 'z', "Porosity")
            # plot.save()

    def single_box(self, field, level=0):
        width = self.prob_domain[2] + 1 - self.prob_domain[0]
        height = self.prob_domain[3] + 1 - self.prob_domain[1]

        field_arr = np.empty([height, width])

        for box_i in range(0, len(self.levels[0]['boxes'])):
            box = self.levels[level]['boxes'][box_i]
            # print(box)

            this_box = self.data[field]['data'][0][box_i]
            # print(this_box)

            ioffset = box[0] - self.prob_domain[0]
            joffset = box[1] - self.prob_domain[1]

            for i in range(0, box[2] + 1 - box[0]):
                for j in range(0, box[3] + 1 - box[1]):
                    # print(str(i) + ', ' + str(j))
                    field_arr[j + joffset][i + ioffset] = this_box[j][i]

        return field_arr


    def set_scale_slice_transform(self, indices, reflect):
        """ Describe how to extract data """
        self.indices = indices
        self.reflect = reflect

    def reset_scale_slice_transform(self):
        self.indices = None
        self.reflect = None


    def get_field_label(self, field_name):

        if field_name in self.FIELD_LABEL.keys():
            return self.FIELD_LABEL[field_name]
        else:
            return field_name













