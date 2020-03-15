import socket
import matplotlib
import numpy as np
import os
import sys
from scipy.signal import find_peaks
from skimage.feature import peak_local_max
from chombopy.plotting import PltFile


class MushyPltFile(PltFile):
    """ Class to load a plot file and perform operations on it
         e.g. find where brine channels are"""

    FIELD_LABEL = {'Porosity': r'$\chi$',
                   'Bulk concentration': r'$\Theta$',
                   'Liquid concentration': r'$\Theta_l$',
                   'Permeability': r'$\Pi$',
                   'xAdvection velocity': '$u$',
                   'yAdvection velocity': '$w$',
                   'streamfunction': r'$\psi$'}

    # def __init__(self):
    #     super().__init__()

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
        if 'parameters.waterDistributionCoeff' in self.inputs.keys():
            pc = float(self.inputs['parameters.waterDistributionCoeff'])
        else:
            pc = 1e-5

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
            for idx, _ in np.ndenumerate(enthalpy):
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
            for idx, _ in np.ndenumerate(enthalpy):
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

            ds_temperature = enthalpy_ds.copy(deep=True).rename('Temperature')
            ds_temperature.values = temperature

            ds_sl = enthalpy_ds.copy(deep=True).rename('Liquid concentration')
            ds_sl.values = liquid_salinity

            ds_ss = enthalpy_ds.copy(deep=True).rename('Solid concentration')
            ds_ss.values = solid_salinity

            self.ds_levels[level]['Porosity'] = ds_porosity
            self.ds_levels[level]['Temperature'] = ds_temperature
            self.ds_levels[level]['Liquid concentration'] = ds_sl
            self.ds_levels[level]['Solid concentration'] = ds_ss

            # available_comps = list(self.ds_levels[0].keys())
            # print('Available comps: %s' % str(available_comps))

    def get_permeability(self, permeability_function='kozeny', rotate_dims=False):

        rotate_dims = self.get_rotate_dims(rotate_dims)

        porosity = self.get_data('Porosity', rotate_dims=rotate_dims)

        if permeability_function == 'kozeny':

            # Cap max porosity just below one to avoid dividing by 0
            porosity = np.clip(porosity, 0, 1 - 10 ** (-10))

            liquid_permeability = porosity ** 3 / (1 - porosity) ** 2

            if 'parameters.heleShawPermeability' in self.inputs.keys():
                hele_shaw_permeability = self.inputs['parameters.heleShawPermeability']
            elif 'parameters.nonDimReluctance' in self.inputs.keys():
                hele_shaw_permeability = 1 / float(self.inputs['parameters.nonDimReluctance'])
            else:
                print('Unable to determine Hele-Shaw permeability, cannot find either parameters.heleShawPermeability '
                      'or parameters.nonDimReluctance in the inputs file. ')
                sys.exit(-1)

            total_permeability = (hele_shaw_permeability ** (-1) + liquid_permeability ** (-1)) ** (-1)

        elif permeability_function == 'cubic':
            total_permeability = porosity ** 3

        else:
            total_permeability = 1.0

        return total_permeability

    def num_channels(self, z_ml):

        bulk_salinity = self.get_level_data('Bulk concentration')

        peak_height_scaling = 2.0
        separation = 5  # minimum pixel separation

        # dom = self.prob_domain
        # min_length = min(dom)
        # separation = float(min_length) /

        if self.space_dim == 2:

            separation = 5  # minimum pixel separation

            bulk_s_slice = bulk_salinity.sel(y=z_ml, method='nearest')
            vel_slice = np.array(self.get_level_data('yDarcy velocity').sel(y=z_ml, method='nearest'))
            slice_arr = np.array(bulk_s_slice) * vel_slice

            prominence = float(slice_arr.max()) / 4.0

            peaks, _ = find_peaks(slice_arr, prominence=prominence, distance=separation)
            num_peaks = len(peaks)

            # print('Num peaks = %d' % num_peaks)
            # x = np.array(slice.coords['x'])
            # import matplotlib.pyplot as plt
            # fig = plt.figure()
            # ax = fig.gca()
            # ax.plot(x, slice_arr)
            # x_peaks = x[peaks]
            # ax.plot(x_peaks, slice_arr[peaks] , 'r.')
            # plt.tight_layout()
            # plt.show()

            return num_peaks
        else:

            bulk_s_slice = bulk_salinity.sel(z=z_ml, method='nearest')

            vel_slice = np.array(self.get_level_data('zDarcy velocity').sel(z=z_ml, method='nearest'))

            slice_arr = np.array(bulk_s_slice) * vel_slice
            # slice_arr = slice_arr - slice_arr.min()

            peak_height = float(slice_arr.max()) / peak_height_scaling
            # print('threshold_abs = %s' % peak_height)

            coordinates = peak_local_max(slice_arr, min_distance=separation,
                                         threshold_abs=peak_height, exclude_border=False)

            num_peaks = len(coordinates)

            # print('Num channels: %d' % num_peaks)
            # import matplotlib.pyplot as plt
            # fig = plt.figure()
            # ax = fig.gca()
            # bulk_s = ax.pcolormesh(slice_arr)
            # ax.plot(coordinates[:, 1], coordinates[:, 0], 'r+')
            # plt.colorbar(bulk_s)
            # plt.tight_layout()
            # plt.show()

            return num_peaks

    def compute_mush_liquid_interface(self):

        porosity = self.get_level_data('Porosity')
        if self.space_dim == 2:
            z = porosity.coords['y']
        else:
            z = porosity.coords['z']

        porosity = np.array(porosity)
        z_interface = np.nan
        mushy_indices = np.argwhere(porosity < 1.0)

        if mushy_indices.shape[0] > 0:
            # Need to be more careful here
            # There can be liquid cells inside/above the mushy layer which we don't want to consider
            # Start from the bottom of the domain and trace cells upwards until we hit mushy cells

            # Rather than highest liquid cell, consider lowest mushy cell
            # iterate over everything but vertical index (which is in the last index)
            iterate_shape = porosity.shape[:-1]
            import itertools
            liquid_z_indices = []
            # for ix in np.arange(0, s[1]):

            for ix in itertools.product(*[range(s) for s in iterate_shape]):
                ind = list(ix)

                # Don't know how to do this in 2d/3d, but should be possible
                if self.space_dim == 2:
                    vals = np.argwhere(porosity[ind[0], :] < 1.0)
                elif self.space_dim == 3:
                    vals = np.argwhere(porosity[ind[0], ind[1], :] < 1.0)
                else:
                    print('Unknown number of dimensions: %d' % self.space_dim)
                    sys.exit(-1)

                if vals.any():
                    val = np.amin(vals)
                    liquid_z_indices.append(val)

            # Eutectic interface:
            solid_indices = []
            # for ix in np.arange(0, s[1]):
            #     vals = np.argwhere(porosity[:, ix] < 1e-6)
            for ix in itertools.product(*[range(s) for s in iterate_shape]):
                ind = list(ix)

                # Don't know how to do this in 2d/3d, but should be possible
                if self.space_dim == 2:
                    vals = np.argwhere(porosity[ind[0], :] < 1e-6)
                elif self.space_dim == 3:
                    vals = np.argwhere(porosity[ind[0], ind[1], :] < 1e-6)
                else:
                    print('Unknown dimension: %d' % self.space_dim)
                    sys.exit(-1)

                if vals.any():
                    val = np.amin(vals)
                    solid_indices.append(val)

            # median_eutectic_boundary = np.median(solid_indices)
            # eutectic_interface = median_eutectic_boundary * dx
            # eutectic_boundary_cell = int(median_eutectic_boundary)

            # liquid_indices = [np.amin(np.argwhere(porosity[:, ix] < 1.0)) for ix in np.arange(0, s[1])]

            if liquid_z_indices:
                median_liquid_boundary = int(np.median(liquid_z_indices))
            else:
                median_liquid_boundary = 0

            # height = abs(dx * float(median_eutectic_boundary - median_liquid_boundary))

            z_interface = z[median_liquid_boundary]
            # ice_ocean_boundary_cell = int(median_liquid_boundary)

            # If there's a mushy layer, compute confinement depth
            # confinement_depth = dx * (np.max(liquid_z_indices) - ice_ocean_interface)
            # data['channel_confinement_depth'] = confinment_depth

        return z_interface

    def skip_component_import_names(self):
        # Some of my files have wierd xEnthalpy yEnthalpy fields, which we should skip
        return ["xEnthalpy", "yEnthalpy"]

    def should_negate_field_upon_reflection(self, field):
        if field[0] == "x" or field == "streamfunction":
            return True

        return False


def latexify(fig_width=None, fig_height=None):
    """Set up matplotlib's RC params for LaTeX plotting.
    Call this before plotting a figure.

    Parameters
    ----------
    fig_width : float, optional, inches
    fig_height : float,  optional, inches

    Note that, in the thesis template, the standard column width is just over 5.8 inches and font size is 12pt
    """

    # code adapted from http://www.scipy.org/Cookbook/Matplotlib/LaTeX_Examples

    # Width and max height in inches for IEEE journals taken from
    # computer.org/cms/Computer.org/Journal%20templates/transactions_art_guide.pdf

    default_fig_width = 5.0

    # Pretty sure this is the standard caption font size in latex
    font_size = 9
    linewidth = 1

    if fig_width is None:
        # fig_width = 3.39 if columns == 1 else 6.9  # width in inches
        fig_width = default_fig_width

    if fig_height is None:
        golden_mean = (np.sqrt(5) - 1.0) / 2.0  # Aesthetic ratio
        fig_height = fig_width * golden_mean  # height in inches

    # Need the mathsrfs package for \mathscr if text.usetex = True
    params = {'text.latex.preamble': ['\\usepackage{gensymb}', '\\usepackage{mathrsfs}', '\\usepackage{amsmath}'],
              'axes.labelsize': font_size, 'axes.titlesize': font_size, 'legend.fontsize': font_size,
              'xtick.labelsize': font_size, 'ytick.labelsize': font_size, 'font.size': font_size,
              'xtick.direction': 'in', 'ytick.direction': 'in', 'lines.markersize': 3, 'lines.linewidth': linewidth,
              'text.usetex': True, 'figure.figsize': [fig_width, fig_height], 'font.family': 'serif', 'backend': 'ps'}

    if 'osx' in socket.gethostname():
        # params['text.usetex'] = False
        params['pgf.texsystem'] = 'pdflatex'
        os.environ['PATH'] = os.environ['PATH'] + ':/Library/TeX/texbin/:/usr/local/bin/'

    matplotlib.rcParams.update(params)


def latexify2(fig_width=5.0, fig_height=4.0):
    """
    Have ended up with two different latex formatting functions, unifying them here. Old version commented below.
    :param fig_width:
    :param fig_height:
    :return:
    """
    latexify(fig_width, fig_height)

    # font_size = 12
    #
    # params = {'backend': 'ps',
    #           'text.latex.preamble': ['\\usepackage{gensymb}', '\\usepackage{mathrsfs}'],
    #           'axes.labelsize': font_size,  # fontsize for x and y labels (was 10)
    #           'axes.titlesize': font_size,
    #           'legend.fontsize': font_size,  # was 10
    #           'xtick.labelsize': font_size,
    #           'ytick.labelsize': font_size,
    #           'lines.markersize': 3,
    #           'lines.linewidth': 1,
    #           'text.usetex': True,
    #           'figure.figsize': [fig_width, fig_height],
    #           'font.family': 'serif'
    #           }
    #
    # mpl.rcParams.update(params)
