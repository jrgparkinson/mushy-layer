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
import sys
from ChkFile import compute_channel_properties
from scipy.signal import find_peaks

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

    # List of names for each direction to be used later
    INDEX_COORDS_NAMES = ['i', 'j', 'k', 'l', 'm']  # add more here if more dimensions

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

            # Data is stored differently in plot files and chk files
            # much of what follows will be different as a result
            if 'data:datatype=0' in list(level_group.keys()):
                self.is_plot_file = True
            else:
                self.is_plot_file = False

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



            # Create a box which spans the whole domain
            # Initialise with a box spanning the whole domain, then add data where it exists
            # Important to do it like this for refined levels, where the whole domain isn't covered with data
            size = []
            for i in range(self.space_dim):
                lev_dom_box_dir = np.arange(self.full_domain_size[i] + lev_dx / 2,
                                            self.full_domain_size[i + self.space_dim] - lev_dx / 2,
                                            lev_dx)
                size.append(lev_dom_box_dir.size)

            blank_data = np.empty(tuple(size))
            blank_data[:] = np.nan

            # Create an empty dataset which spans entire domain, which we will use as a template for loading data into
            coords = {}
            box_size = ()
            for d in range(self.space_dim):
                coords_dir = np.arange(self.prob_domain[d], self.prob_domain[self.space_dim + d] + 1)
                coords[self.INDEX_COORDS_NAMES[d]] = coords_dir
                box_size = box_size + (coords_dir.size,)  # append to tuple of sizes

            blank_data = np.empty(box_size)
            ds_dom_box = xr.Dataset({}, coords=coords)

            # Use indexes rather than x, y for now - then convert to x,y later
            # this is to avoid issues with floating point arithmetic when merging datasets
            # (we can end up trying to merge datasets where x coordinates differ by ~ 10^{-10}, creating nonsense)

            # Create empty datasets spanning entire domain for each component
            for comp_name in self.comp_names:
                s = blank_data.shape

                # Last component should be x-direction
                if not s[-1] == len(coords['i']):
                    blank_data = blank_data.T

                extended_coords = coords
                extended_coords['level'] = level
                dims = self.INDEX_COORDS_NAMES[:self.space_dim]
                dims = dims[::-1]  # reverse list so we have k, j, i etc
                ds_dom_box[comp_name] = xr.DataArray(blank_data, dims=dims,  # dims=['j', 'i'],
                                                     coords=extended_coords)

            ds_boxes = [ds_dom_box]


            # Get level outlines
            polygons = []
            for box in self.levels[level][self.BOXES]:
                lo_indices = [box[i] for i in range(self.space_dim)]
                hi_indices = [box[i] for i in range(self.space_dim, 2 * self.space_dim)]

                # 0.5 because cell centred
                lo_vals = [lev_dx * (0.5 + i) for i in lo_indices]
                hi_vals = [lev_dx * (0.5 + i) for i in hi_indices]

                end_points = [[lo_vals[i] - lev_dx / 2, hi_vals[i] + lev_dx / 2] for i in range(self.space_dim)]

                # For plotting level outlines
                from itertools import product, permutations, chain
                # Construct vertices in n dimensions
                polygon_vertices_auto = list(product(*end_points))
                polygon_vertices_auto = sorted(polygon_vertices_auto,
                                               key=lambda x: np.arctan(x[1] / max(abs(x[0]), 0.0001)))

                poly = Polygon(polygon_vertices_auto)
                if poly.is_valid:
                    polygons.append(poly)

            level_outline = gpd.GeoSeries(cascaded_union(polygons))
            self.level_outlines.append(level_outline)

            # Data is sorted by box and by component, so need to know total number of components
            num_comps = 0
            for comp_name in self.data.keys():
                num_comps = num_comps + self.data[comp_name][self.NUM_COMPS]


            # Now get the  data for each field, on each level

            if self.is_plot_file:

                # For plt files, data is sorted by box then by component

                data = level_group['data:datatype=0']

                # Some other stuff we can get, but don't at the moment:
                # data_offsets = level_group['data:offsets=0']
                # data_atts = level_group['data_attributes']
                # advVel = level_group['advVel:datatype=0']
                # advVel_offsets = level_group['advVel:offsets=0']
                # advVel_atts = level_group['advVel_attributes']

                data_unshaped = data[()]



                offset = 0

                for box in self.levels[level][self.BOXES]:
                    lo_indices = [box[i] for i in range(self.space_dim)]
                    hi_indices = [box[i] for i in range(self.space_dim, 2 * self.space_dim)]


                    n_cells_dir = [hi_indices[d] + 1 - lo_indices[d] for d in range(self.space_dim)]

                    num_box_cells = np.prod(n_cells_dir)  # effectively nx * ny * nz * ...
                    num_cells = num_box_cells * num_comps  # also multiply by number of components

                    # Now split into individual components
                    # data contains all fields on this level, sort into individual fields
                    comp_offset_start = 0

                    coords = {}
                    # box_size = ()
                    for d in range(self.space_dim):
                        coords_dir = np.arange(lo_indices[d], hi_indices[d] + 1)
                        coords[self.INDEX_COORDS_NAMES[d]] = coords_dir

                    # Blank dataset for this box, which each component will be added to
                    ds_box = xr.Dataset({}, coords=coords)

                    for comp_name in self.comp_names:

                        # print('Num cells in a box: ' + str(num_box_cells))
                        comp_offset_finish = comp_offset_start + num_box_cells

                        indices = [offset+comp_offset_start, offset + comp_offset_finish]

                        comp_offset_start = comp_offset_finish

                        ds_box[comp_name] = self.get_box_comp_data(data_unshaped, level, indices, comp_name, n_cells_dir, coords)

                    # Move onto next box
                    offset = offset + num_cells

                    ds_boxes.append(ds_box)




            else:

                # For chk files, data is sotred by component then by box
                # Loop over components and get data for that component from each box


                box_offset = 0


                for box in self.levels[level][self.BOXES]:
                    lo_indices = [box[i] for i in range(self.space_dim)]
                    hi_indices = [box[i] for i in range(self.space_dim, 2 * self.space_dim)]


                    n_cells_dir = [hi_indices[d] + 1 - lo_indices[d] for d in range(self.space_dim)]

                    num_box_cells = np.prod(n_cells_dir)  # effectively nx * ny * nz * ...
                    # num_cells = num_box_cells * num_comps  # also multiply by number of components

                    # Now split into individual components
                    # data contains all fields on this level, sort into individual fields
                    comp_offset_start = 0

                    coords = {}
                    # box_size = ()
                    for d in range(self.space_dim):
                        coords_dir = np.arange(lo_indices[d], hi_indices[d] + 1)
                        coords[self.INDEX_COORDS_NAMES[d]] = coords_dir

                    # Blank dataset for this box, which each component will be added to
                    ds_box = xr.Dataset({}, coords=coords)



                    for comp_name in self.comp_names:
                        component = 0

                        # Need to get data differently if this is a vector
                        is_vector = (comp_name[0] == 'x' or comp_name[0] == 'y' or comp_name[
                            0] == 'z' and  sum([comp_name[1:] in x for x in self.comp_names]) == self.space_dim)
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

                        data_unshaped = data[()]

                        # print('Num cells in a box: ' + str(num_box_cells))
                        #comp_offset_finish = comp_offset_start + num_box_cells

                        component_offset = 0 # this is just 0 because we only get data for this component

                        # For vectors, we may have an offset?
                        if num_comps > 1:
                            component_offset = component*num_box_cells

                        start_index = component_offset + box_offset
                        end_index = start_index + num_box_cells

                        indices = [start_index, end_index]

                        #comp_offset_start = comp_offset_finish

                        ds_box[comp_name] = self.get_box_comp_data(data_unshaped, level, indices, comp_name, n_cells_dir, coords)

                        # component_offset = component_offset + (num_box_cells*num_comps)

                    # Move onto next box
                    box_offset = box_offset + num_box_cells

                    ds_boxes.append(ds_box)



            # ds_level = xr.merge(ds_boxes)

            # Update will replace in place
            first_box = 1
            ds_level = ds_boxes[first_box]
            for b in ds_boxes[first_box+1:]:

                ds_level =  ds_level.combine_first(b)


            # Create x,y,z, coordinates
            x_y_coords_names = ['x', 'y', 'z']

            for d in range(self.space_dim):
                ds_level.coords[x_y_coords_names[d]] = ds_level.coords[self.INDEX_COORDS_NAMES[d]] * lev_dx
            # ds_level.coords['x'] = ds_level.coords['i'] * lev_dx
            # ds_level.coords['y'] = ds_level.coords['j'] * lev_dx

            if zero_x:
                ds_level.coords['x'] = ds_level.coords['x'] - min(ds_level.coords['x'])

            # Swap i,j,k to x,y,z coordinates
            for d in range(self.space_dim):
                ds_level = ds_level.swap_dims({self.INDEX_COORDS_NAMES[d]: x_y_coords_names[d]})

            # ds_level = ds_level.swap_dims({'i': 'x'})
            #ds_level = ds_level.swap_dims({'j': 'y'})


            #TODO: should level be an attribute or co-ordinate? need to try with actual AMR data
            ds_level.attrs['level'] =  level


            ds_levels.append(ds_level)


        self.ds_levels = ds_levels

        h5_file.close()

    def get_box_comp_data(self, data_unshaped, level, indices, comp_name, n_cells_dir, coords):

        data_box_comp = data_unshaped[indices[0]:indices[1]]

        reshaped_data = data_box_comp.reshape(tuple(n_cells_dir))

        # I think we only need to transpose in 3D
        if self.space_dim == 3 or self.space_dim == 2:
            reshaped_data = reshaped_data.transpose()
            
        reshaped_data = np.array(reshaped_data)

        # Check if scalar or vector
        trimmed_comp_names = [n[1:] for n in self.comp_names]
        field_type = 'scalar'
        if comp_name[0] in ('x', 'y', 'z') and comp_name[1:] in trimmed_comp_names:
            field_type = 'vector'

        dim_list = self.INDEX_COORDS_NAMES[:self.space_dim]
        # dim_list = dim_list[::-1]
        extended_coords = coords
        extended_coords['level'] = level

        xarr_component_box = xr.DataArray(reshaped_data, dims=dim_list,  # ['j', 'i'],
                                          coords=extended_coords,  # {'i': i_box, 'j': j_box, 'level': level},
                                          attrs={'field_type': field_type})



        return xarr_component_box

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

        porosity = np.array(self.get_level_data('Porosity'))

        return compute_channel_properties(porosity, do_plots)


        # Reconstruct level 0 porosity as single np array
        # width = self.prob_domain[2] + 1 - self.prob_domain[0]
        # # height = self.prob_domain[3] + 1 - self.prob_domain[1]
        #
        # # porosity = np.empty([height, width])
        #
        # # porosity = self.single_box('Porosity')
        # porosity = self.get_level_data('Porosity')
        #
        # # Iterate over porosity field
        # cols = porosity.shape[0]
        # rows = porosity.shape[1]
        # channels = [None] * cols
        # chimney_positions = []
        # for j in range(cols):
        #     # chimneys in row
        #     chimneys_in_row = 0
        #     average_porosity = 0
        #
        #     currently_liquid = False
        #     chimney_pos_row = []
        #
        #     for i in range(rows):
        #         # print porosity[i,j]
        #         chi = porosity[j, i]
        #         average_porosity = average_porosity + chi
        #         if chi < 1.0 and currently_liquid:
        #             # Just left liquid region
        #             currently_liquid = False
        #
        #         elif chi > 0.999 and not currently_liquid:
        #             # Entered a liquid region
        #             chimneys_in_row = chimneys_in_row + 1
        #             chimney_pos_row.append(i)
        #             currently_liquid = True
        #
        #     average_porosity = average_porosity / rows
        #     if average_porosity > 0.99:
        #         channels[j] = 0  # This region was entirely liquid, can't have a channel
        #     else:
        #         channels[j] = chimneys_in_row
        #
        #         # First add chimney positions for chimneys we've already found
        #         # print(str(chimney_positions))
        #         # print(str(len(chimney_positions)))
        #
        #         if chimneys_in_row == 0:
        #             continue
        #
        #         for chimney_i in range(0, len(chimney_positions)):
        #             if chimney_i < len(chimney_pos_row):
        #                 this_chimney_position_in_row = chimney_pos_row[chimney_i]
        #                 chimney_positions[chimney_i].append(this_chimney_position_in_row)
        #
        #         # Now, add extra rows to chimney_positions vector if needed
        #         for chimney_i in range(len(chimney_positions), chimneys_in_row):
        #             chimney_positions.append([chimney_pos_row[chimney_i]])
        #
        # # Get an idea of the channel depths
        # chan_depths = []
        # for i in range(0, len(chimney_positions)):
        #     chan_depths.append(len(chimney_positions[i]))
        #
        # max_chan_depth = 0
        # if chan_depths:
        #     max_chan_depth = max(chan_depths)
        #
        # average_chan_positions = []
        # rel_chan_positions = []
        # num_channels = 0
        # for i in range(0, len(chimney_positions)):
        #     if chan_depths[i] > 0.5 * max_chan_depth:
        #         average_chan_positions.append(np.mean(chimney_positions[i]))
        #         rel_chan_positions.append(average_chan_positions[-1] / width)
        #         num_channels = num_channels + 1
        #     # else:
        #     # print('Discarding channel as too short: '  + str(chimney_positions[i]))
        # # print('Number of channels: ' + str(num_channels))
        # # print('Channel positions: ' + str(average_chan_positions))
        # # print('Relative channel positions: ' + str(rel_chan_positions))
        #
        # # doPlots = False
        # if do_plots:
        #     print('Making plot')
        #     self.plot_field('Porosity')
        #
        # channel_spacing = [(rel_chan_positions[i] - rel_chan_positions[i - 1]) for i in
        #                    range(1, len(rel_chan_positions))]
        # # channel_spacing = channel_spacing*self.levels[0][DX]
        #
        # return [num_channels, rel_chan_positions, channel_spacing]


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


    def get_mesh_grid(self, level=0, rotate_dims=False):

        dx = self.levels[level][self.DX]


        components = list(self.data.keys())
        # components = [c.decode('UTF-8') for c in components]

        field_array = self.single_box(components[0])
        grid_size = field_array.shape




        # Make sure the grid stretches from the start of the first cell to the end of the last cell
        # this means we stretch dx slightly out of proportion, but it ensures plot limits are correct
        x_max = (grid_size[0] + 1) * dx
        y_max = (grid_size[1] + 1) * dx

        grid_dx = x_max / grid_size[0]
        grid_dy = y_max / grid_size[1]



        y, x = np.mgrid[slice(0, x_max, grid_dx),
                        slice(0, y_max, grid_dy)]

        coord_max = [(grid_size[i] + 1) * dx for i in range(0, self.space_dim)]
        grid_spacing = [coord_max[i] / grid_size[i] for i in range(0, self.space_dim)]
        grids = np.mgrid[[slice(0, coord_max[i], grid_spacing[i]) for i in range(0, self.space_dim)]]

        if self.space_dim == 3:
            x = grids[0]
            y = grids[1]
            z = grids[2]

            if rotate_dims:
                x = x.transpose()
                y = y.transpose()
                z = z.transpose()

            z = self.scale_slice_transform(z)
            x = self.scale_slice_transform(x, no_reflect=True)
            y = self.scale_slice_transform(y, no_reflect=True)

            return x, y, z


        else:


            if rotate_dims:
                y_new = x.transpose()
                x_new = y.transpose()

                x = x_new
                y = y_new

            y = self.scale_slice_transform(y)
            x = self.scale_slice_transform(x, no_reflect=True)

            return x, y

    # Added for compatibility with ChkFile interface
    def get_data(self, var_name, rotate_dims=False):
        data = self.get_level_data(var_name)

        if data is None:
            return None

        data = np.array(data)

        if rotate_dims:
            data = data.transpose()

        return data

    def get_level_data(self, field, level=0, valid_only=False):

        if self.data_load_method == self.YT:
            pass
            #TODO: write this

        if not self.data_loaded:
            print('Data not loaded')
            return

        available_comps = list(self.ds_levels[0].keys())

        if field not in available_comps:
            print('Field: %s not found. The following fields do exist: ' % field)
            # print(self.data.keys())
            print(available_comps)

            if field == 'LiquidSalinity' and 'Liquid concentration' in available_comps:
                print('Using Liquid concentration instead')
                field = 'Liquid concentration'
            else:
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
        # width = self.prob_domain[2] + 1 - self.prob_domain[0]
        # height = self.prob_domain[3] + 1 - self.prob_domain[1]
        #
        # field_arr = np.empty([height, width])
        #
        # for box_i in range(0, len(self.levels[0]['boxes'])):
        #     box = self.levels[level]['boxes'][box_i]
        #     # print(box)
        #
        #     this_box = self.data[field]['data'][0][box_i]
        #     # print(this_box)
        #
        #     ioffset = box[0] - self.prob_domain[0]
        #     joffset = box[1] - self.prob_domain[1]
        #
        #     for i in range(0, box[2] + 1 - box[0]):
        #         for j in range(0, box[3] + 1 - box[1]):
        #             # print(str(i) + ', ' + str(j))
        #             field_arr[j + joffset][i + ioffset] = this_box[j][i]
        ds = self.get_level_data(field, level)
        field_arr = np.array(ds)

        return field_arr


    def set_scale_slice_transform(self, indices, reflect=False):
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


        for level in range(0, self.num_levels):
            enthalpy_ds = self.get_level_data('Enthalpy', level)
            bulk_salinity_ds = self.get_level_data('Bulk concentration', level)

            enthalpy_eutectic = enthalpy_ds.copy(deep=True).rename('Enthalpy eutectic')
            enthalpy_solidus = enthalpy_ds.copy(deep=True).rename('Enthalpy solidus')
            enthalpy_liquidus = enthalpy_ds.copy(deep=True).rename('Enthalpy liquidus')

            # porosity = enthalpy.copy(deep=True).rename('Porosity')
            # temperature = enthalpy.copy(deep=True).rename('Temperature')
            # liquid_salinity = enthalpy.copy(deep=True).rename('Liquid concentration')
            # solid_salinity = enthalpy.copy(deep=True).rename('Solid concentration')

            enthalpy = np.array(enthalpy_ds)
            bulk_salinity = np.array(bulk_salinity_ds)

            enthalpy_eutectic = np.array(enthalpy_eutectic)
            enthalpy_solidus = np.array(enthalpy_solidus)
            enthalpy_liquidus = np.array(enthalpy_liquidus)

            porosity = np.empty(enthalpy.shape)
            temperature = np.empty(enthalpy.shape)
            liquid_salinity = np.empty(enthalpy.shape)
            solid_salinity = np.empty(enthalpy.shape)

            # Now compute bounding energies
            print('Computing bounding energy')
            cols = enthalpy.shape[0]
            rows = enthalpy.shape[1]
            # for j in range(cols):
            #
            #     for i in range(rows):
            # replaced loop over [j,i] with loop over [idx] for 3D compatibility

            for idx, value in np.ndenumerate(enthalpy):

                eutectic_porosity = (conc_ratio + bulk_salinity[idx]) / (theta_eutectic + conc_ratio)
                enthalpy_eutectic[idx] = eutectic_porosity * (stefan + theta_eutectic * (1 - cp)) + cp * theta_eutectic
                enthalpy_solidus[idx] = cp * (theta_eutectic + max(0.0, (-bulk_salinity[idx] - conc_ratio) / pc))
                enthalpy_liquidus[idx] = stefan - bulk_salinity[idx] + theta_eutectic + theta_eutectic

            # eutectic_porosity = (conc_ratio + bulk_salinity) / (theta_eutectic + conc_ratio)
            # enthalpy_eutectic = eutectic_porosity * (stefan + theta_eutectic * (1 - cp)) + cp * theta_eutectic
            # enthalpy_solidus = cp * (theta_eutectic + np.max(0.0, (-bulk_salinity - conc_ratio) / pc))
            # enthalpy_liquidus = stefan - bulk_salinity + theta_eutectic + theta_eutectic

            print('Computing diagnostic variables')
            # Compute diagnostic variables
            # for j in range(cols):
            #
            #     for i in range(rows):

            # replaced loop over [j,i] with loop over [idx] for 3D compatibility
            for idx, value in np.ndenumerate(enthalpy):
                if enthalpy[idx] <= enthalpy_solidus[idx]:
                    porosity[idx] = 0.0
                    temperature[idx] = enthalpy[idx] / cp
                    liquid_salinity[idx] = 0.0
                    solid_salinity[idx] = bulk_salinity[idx]
                elif enthalpy[idx] <= enthalpy_eutectic[idx]:
                    porosity[idx] = (enthalpy[idx] - theta_eutectic * cp) / (stefan + theta_eutectic * (1 - cp))
                    temperature[idx] = theta_eutectic
                    liquid_salinity[idx] = theta_eutectic
                    solid_salinity[idx] = bulk_salinity[idx] / (1 - porosity[idx])
                elif enthalpy[idx] < enthalpy_liquidus[idx]:
                    a = conc_ratio * (cp - 1) + stefan * (pc - 1)
                    b = conc_ratio * (1 - 2 * cp) + enthalpy[idx] * (1 - pc) \
                        - bulk_salinity[idx] * (cp - 1) - pc * stefan
                    c = (bulk_salinity[idx] + conc_ratio) * cp + pc * enthalpy[idx]
                    porosity[idx] = (-b - np.sqrt(b * b - 4 * a * c)) / (2 * a)

                    liquid_salinity[idx] = (bulk_salinity[idx] + conc_ratio * (1 - porosity[idx])) / (
                            porosity[idx] + pc * (1 - porosity[idx]))
                    temperature[idx] = -liquid_salinity[idx]
                    solid_salinity[idx] = (pc * bulk_salinity[idx] - conc_ratio * porosity[idx]) / (
                            porosity[idx] + pc * (1 - porosity[idx]))
                else:
                    porosity[idx] = 1.0
                    temperature[idx] = enthalpy[idx] - stefan
                    liquid_salinity[idx] = bulk_salinity[idx]
                    solid_salinity[idx] = 0.0

            # Finally, save data into boxes

            ds_porosity = enthalpy_ds.copy(deep=True).rename('Porosity')
            ds_porosity.values = porosity

            ds_T = enthalpy_ds.copy(deep=True).rename('Temperature')
            ds_T.values = temperature

            ds_sl = enthalpy_ds.copy(deep=True).rename('Liquid concentration')
            ds_sl.values = liquid_salinity

            ds_ss = enthalpy_ds.copy(deep=True).rename('Solid concentration')
            ds_ss.values = solid_salinity


            self.ds_levels[level]['Porosity'] = ds_porosity
            self.ds_levels[level]['Temperature'] = ds_T
            self.ds_levels[level]['Liquid concentration'] = ds_sl
            self.ds_levels[level]['Solid concentration'] = ds_ss


            # available_comps = list(self.ds_levels[0].keys())
            # print('Available comps: %s' % str(available_comps))

    def get_permeability(self, permeability_function='kozeny', rotate_dims=False):

        porosity = self.get_data('Porosity', rotate_dims=rotate_dims)
        permeability = np.empty(porosity.shape)

        if permeability_function == 'kozeny':

            # Cap max porosity just below one to avoid dividing by 0
            porosity = np.clip(porosity, 0, 1 - 10 ** (-10))

            liquid_permeability = porosity ** 3 / (1 - porosity) ** 2

            hele_shaw_permeability = 1 / float(self.inputs['parameters.nonDimReluctance'])

            total_permeability = (hele_shaw_permeability ** (-1) + liquid_permeability ** (-1)) ** (-1)

        elif permeability_function == 'cubic':
            total_permeability = porosity ** 3

        else:
            total_permeability = 1.0

        return total_permeability


    def num_channels(self, z_ml):

        bulk_salinity = self.get_level_data('Bulk concentration')

        peak_height_scaling = 2.0
        separation = 2 # minimum pixel separation

        # dom = self.prob_domain
        # min_length = min(dom)
        # separation = float(min_length) /


        if self.space_dim == 2:
            slice = bulk_salinity.sel(y = z_ml, method='nearest')
            slice_arr = np.array(slice)


            peak_height = float(slice_arr.max()) / peak_height_scaling

            peaks, _ = find_peaks(slice_arr, height=peak_height, distance=separation)

            num_peaks = len(peaks)

            return num_peaks
        else:
            from skimage.feature import peak_local_max

            slice = bulk_salinity.sel(z = z_ml, method='nearest')

            slice_arr = np.array(slice)
            slice_arr = slice_arr - slice_arr.min()

            peak_height = float(slice_arr.max()) / peak_height_scaling
            # print('threshold_abs = %s' % peak_height)

            coordinates = peak_local_max(slice_arr, min_distance=separation, threshold_abs=peak_height)

            num_peaks = len(coordinates)
            # print('Num channels: %d' % num_peaks)

            # import matplotlib.pyplot as plt
            # fig = plt.figure()
            # ax = fig.gca()
            # bulk_s = ax.pcolormesh(slice_arr)
            # ax.plot(coordinates[:, 1], coordinates[:, 0], 'r.')
            #
            # plt.colorbar(bulk_s)
            # plt.tight_layout()
            # plt.show()

            return num_peaks

        return np.nan












