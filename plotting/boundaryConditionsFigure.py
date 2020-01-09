# Script to load an inputs file and make a figure
# describing the boundary conditions being used
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from PltFile import latexify
from mushyLayerRunUtils import read_inputs, string_to_array
import sys
import getopt


def make_bc_fig(inputs_file, ndim=2):

    inputs = read_inputs(inputs_file)

    n_cells = string_to_array(inputs['main.num_cells'])
    nx = n_cells[0]
    ny = n_cells[1]

    max_dom_length = max(nx, ny)
    scaled_nx = nx / max_dom_length
    scaled_ny = ny / max_dom_length

    # Start making figure
    window_width = 9.0
    window_height = window_width *0.9* scaled_ny / scaled_nx
    # window_height=12.0
    latexify(fig_width=window_width, fig_height=window_height)
    fig = plt.figure()


    # Add domain rectangle
    ax_width = 0.5
    ax_height = ax_width*(scaled_ny/scaled_nx)*(window_width/window_height)
    ax = fig.add_axes([(1-ax_width)/2.0, 0.25, ax_width, ax_height])
    ax.set_axis_off()


    domain_box = patches.Rectangle((0, 0), 1, 1, fill=False, transform = ax.transAxes, clip_on=False )

    ax.add_patch(domain_box)

    # Add BCs to each side
    padding = 0.02
    for dim in range(0, ndim):
        for side in range(0, 2):

            # Defaults
            horiz_align = 'center'
            vert_align = 'top'
            rotate = 0

            plus_minus = 1
            if side == 0:
                plus_minus = -1

            if dim == 0:
                # X direction bcs (left/right)
                # rotate = 90
                # ypos = 0.0

                rotate = 0
                ypos = scaled_ny/2.0
                vert_align = 'bottom'

                if side == 0:
                    horiz_align = 'right'
                else:
                    horiz_align = 'left'
            else:
                ypos = side + padding*plus_minus

            if dim == 1:
                # Y direction bcs (top/bottom)
                xpos = scaled_nx/2.0

                rotate = 0
                horiz_align = 'center'

                if side == 0:
                    vert_align = 'top'
                else:
                    vert_align = 'bottom'

            else:
                xpos = side + padding*plus_minus


            bc_text = make_bc_text(inputs, dim, side)

            ax.text(xpos, ypos, bc_text, horizontalalignment=horiz_align, verticalalignment=vert_align,
                    rotation = rotate, transform=ax.transAxes)

    # TODO: label axis extents and add num cells label (e.g 64x64)

    # Also add dimensionless parameters to the middle of the domain
    dim_params = 'Dynamics: $Rm_S = %g, Rm_T = %g$ \n     $\Pi_H = %g, Da = %g, Pr=%g$, \n ' \
                 'Material properties: $Le = %g, c_p = %g, k = %g$, \n' \
                 'Thermodynamics: $\mathscr{C}=%g, \mathscr{S}=%g,$ \n' \
                 'Phase diagram: $\Gamma=%g, C_i = %g,$ \n     $C_e = %g, T_e=%g$ \n' % (
                     inputs['parameters.rayleighComp'],
                     inputs['parameters.rayleighTemp'],
                     1.0 / inputs['parameters.nonDimReluctance'],
                     inputs['parameters.darcy'],
                     inputs['parameters.prandtl'],
                     inputs['parameters.lewis'],
                     inputs['parameters.specificHeatRatio'],
                     inputs['parameters.heatConductivityRatio'],
                     inputs['parameters.compositionRatio'],
                     inputs['parameters.stefan'],
                     inputs['parameters.liquidusSlope'],
                     inputs['parameters.initialComposition'],
                     inputs['parameters.eutecticComposition'],
                     inputs['parameters.eutecticTemp']
                 )

    if 'heatSource.size' in inputs:
        dim_params = dim_params + '+ heat source $Q = \\frac{Q_0}{\sigma \sqrt{2 \pi}} ' \
                                  '\exp\left[ - 0.5 \left( \\frac{x-x_c}{\sigma} \\right)^2 \\right]  ' \
                                  '0.5 \left( 1 + \\tanh\left[10 (z-(H-h)) \\right]) \\right)$, \n' \
                                  'where $Q_0=%s,\sigma=%s,x_c=%s,h=%s$' % (
                     inputs['heatSource.size'], inputs['heatSource.width'], inputs['heatSource.xpos'],
                     inputs['heatSource.depth'])

    ax.text(0.5, 0.4, dim_params, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

    dom_height = float(inputs['main.domain_height'])
    dom_width = dom_height * n_cells[0]/n_cells[1]
    ax.text(0.5, 0.85, 'Domain: [%g, %g] with %d x %d cells' % (dom_width, dom_height, n_cells[0], n_cells[1]), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

    figure_full_path = inputs_file + '-auto-generated-visualisation.pdf'
    print('Saved to %s' % figure_full_path)
    plt.savefig(figure_full_path, format='pdf')

    plt.show()



def make_bc_text(inputs, directory, side):
    scalar_options = ['Dirichlet', 'Neumann', 'InflowOutflow', 'OnlyInflow', 'Robin', 'VariableFlux',
                      'FixedTemperature', 'TemperatureFlux', 'TemperatureFluxRadation']
    scalars = {'enthalpy': 'H', 'bulkConcentration': '\Theta', 'vel': '\mathbf{U}'}
    BC_TYPES = {'bulkConcentration': scalar_options,
                'enthalpy': scalar_options,
                'vel': ['$\mathbf{U} = 0$', 'Inflow',
                        'Outflow',
                        'OutflowNormal',  # only a normal velocity
                        'InflowOutflow',  # both inflow and outflow possible
                        'noShear',
                        'Symmetry ($\mathbf{U} \cdot \mathbf{n} = 0$) ',
                        'Plume inflow',
                        'Outflow with enforced pressure gradient',
                        'Pressure head']}

    sides = ['Lo', 'Hi']
    dirs = ['x', 'y']

    side_text = sides[side]

    #bc_text = dirs[dir] + sides[side]

    variables = ['bulkConcentration', 'enthalpy', 'vel']

    var_texts = []

    for v in variables:
        bc_type_name = 'bc.%s%s' % (v, side_text)
        bc_val_name = 'bc.%s%sVal' % (v, side_text)

        bc_type = string_to_array(inputs[bc_type_name])
        bc_type_this_dir = bc_type[directory]

        bc_type_description = BC_TYPES[v][bc_type_this_dir]

        # Now we know what the BC is, need to display it sensibly

        # Default:
        this_var_text = bc_type_description

        if v == 'vel':
            this_var_text = bc_type_description

        else:

            bc_val = string_to_array(inputs[bc_val_name], conversion=lambda x: float(x))[directory]

            perp_dir = directory + 1
            if perp_dir > 1:
                perp_dir = 0
            perp_dir_string = dirs[perp_dir]

            if bc_type_description == 'Dirichlet':
                this_var_text = 'Fixed: $%s = %.2g$' % (scalars[v], bc_val)
            elif bc_type_description == 'Neumann':
                this_var_text = 'No flux: $\mathbf{n} \cdot \\nabla %s = 0$' % scalars[v]
            elif bc_type_description == 'VariableFlux':

                this_var_text = 'Variable flux: $\mathbf{n} \cdot \\nabla %s = $ \n' \
                                '$%g (1+\\textrm{tanh}(50(%s-0.75)))$' % (scalars[v], bc_val*0.5, perp_dir_string)

            elif bc_type_description == 'FixedTemperature':
                no_flux_limit = string_to_array(inputs['bc.NoFluxLimit%s' % side_text], conversion=lambda x: float(x))[directory]
                this_var_text = '$ T = %g \; (%s > %g),$ \n $ \mathbf{n} \cdot \\nabla T = 0 \; (%s < %g)$' % (
                                                                                                                bc_val,
                                                                                                                perp_dir_string,
                                                                                                                no_flux_limit,
                                                                                                                perp_dir_string,
                                                                                                                no_flux_limit)

            elif bc_type_description == 'TemperatureFlux':
                no_flux_limit = string_to_array(inputs['bc.NoFluxLimit%s' % side_text], conversion=lambda x: float(x))[directory]
                this_var_text = '$ \mathbf{n} \cdot \\nabla  T = %g \; (%s > %g),$ \n $ \mathbf{n} \cdot \\nabla T = 0 \; (%s < %g)$' % (
                                                                                                                bc_val,
                                                                                                                perp_dir_string,
                                                                                                                no_flux_limit,
                                                                                                                perp_dir_string,
                                                                                                                no_flux_limit)


        var_texts.append(this_var_text)

    bc_text = ', \n'.join(var_texts)


    return bc_text


if __name__ == "__main__":

    inputs_file = '/home/parkinsonjl/mushy-layer/execSubcycle/enceladus/inputs'
    inputs_file = shared_storage.get_dir('enceladus/256x128-Rm200-HeatSourceSize0.2-TopBCFixedH-similarToNoFlow/inputs')

    arg = sys.argv[1:]
    try:
        opts, args = getopt.getopt(arg, "i:")
    except getopt.GetoptError:
        print(
            'analyseFixedChill.py -i <inputs file=?>')
        sys.exit(2)
    for opt, arg in opts:
        if opt in "-i":
            inputs_file = str(arg)

    print('Make BC figure for inputs file:\n%s' % inputs_file)

    make_bc_fig(inputs_file)
