import numpy as np

import solcore.analytic_solar_cells as ASC
from solcore.light_source import LightSource
from solcore.state import State
from solcore.optics import solve_beer_lambert, solve_tmm, solve_rcwa, rcwa_options, solve_external_optics
from solcore.absorption_calculator import RCWASolverError
from solcore.structure import Layer, Junction, TunnelJunction

try:
    import solcore.poisson_drift_diffusion as PDD

    a = PDD.pdd_options
except AttributeError:
    PDD.pdd_options = {}

default_options = State()
pdd_options = PDD.pdd_options
asc_options = ASC.db_options

def merge_dicts(*dict_args):
    """
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    """
    result = State()
    for dictionary in dict_args:
        result.update(dictionary)
    return result


# General
default_options.T_ambient = 298
default_options.T = 298

# Illumination spectrum
default_options.wavelength = np.linspace(300, 1800, 251) * 1e-9
default_options.light_source = LightSource(source_type='standard', version='AM1.5g', x=default_options.wavelength,
                                           output_units='photon_flux_per_m')

# IV control
default_options.voltages = np.linspace(0, 1.2, 100)
default_options.mpp = False
default_options.light_iv = False
default_options.internal_voltages = np.linspace(-6, 4, 1000)
default_options.position = None
default_options.radiative_coupling = False

# Optics control
default_options.optics_method = 'BL'
default_options.recalculate_absorption = False

default_options = merge_dicts(default_options, ASC.db_options, PDD.pdd_options, rcwa_options)


def solar_cell_solver(solar_cell, task, user_options=None):
    """ Solves the properties of a solar cell object, either calculating its optical properties (R, A and T), its quantum efficiency or its current voltage characteristics in the dark or under illumination. The general options for the solvers are passed as dicionaries.

    :param solar_cell: A solar_cell object
    :param task: Task to perform. It has to be "optics", "iv", "qe", "equilibrium" or "short_circuit". The last two only work for PDD junctions
    :param user_options: A dictionary containing the options for the solver, which will overwrite the default options.
    :return: None
    """
    if type(user_options) in [State, dict]:
        options = merge_dicts(default_options, user_options)
    else:
        options = merge_dicts(default_options)

    prepare_solar_cell(solar_cell, options)
    options.T = solar_cell.T

    if task == 'optics':
        solve_optics(solar_cell, options)
    elif task == 'iv':
        solve_iv(solar_cell, options)
    elif task == 'qe':
        solve_qe(solar_cell, options)
    elif task == 'equilibrium':
        solve_equilibrium(solar_cell, options)
    elif task == 'short_circuit':
        solve_short_circuit(solar_cell, options)
    else:
        raise ValueError(
            'ERROR in "solar_cell_solver":\n\tValid values for "task" are: "optics", "iv", "qe", "equilibrium" and "short_circuit".')


def solve_optics(solar_cell, options):
    """ Solves the optical properties of the structure, calculating the reflectance, absorptance and transmitance. The "optics_method" option controls which method is used to calculate the optical properties of the solar cell:

    - None: The calculation is skipped. Only useful for solar cells involving just "2-diode" kind of junctions.
    - BL: Uses the Beer-Lambert law to calculate the absorption in each layer. Front surface reflexion has to provided externally. It is the default method and the most flexible one.
    - TMM: Uses a transfer matrix calculation to obtain the RAT. Not valid for DB or 2D junction
    - RCWA: Uses the rigorous wave coupled analysisto obtain the RAT. This allows to include 2D photonic crystals in the structure, for example. Not valid for DB or 2D junctions
    - external: The reflection and absorption profiles are provided externally by the user, and therefore no calculation is performed by Solcore.

    :param solar_cell: A solar_cell object
    :param options: Options for the optics solver
    :return: None
    """
    print('Solving optics of the solar cell...')

    calculated = hasattr(solar_cell[0], 'absorbed')
    recalc = options.recalculate_absorption if 'recalculate_absorption' in options.keys() else False
    if not calculated or recalc:

        if options.optics_method is None:
            print('Warning: Not solving the optics of the solar cell.')
        elif options.optics_method == 'external':
            solve_external_optics(solar_cell, options)
        elif options.optics_method == 'BL':
            solve_beer_lambert(solar_cell, options)
        elif options.optics_method == 'TMM':
            solve_tmm(solar_cell, options)
        elif options.optics_method == 'RCWA':
            if solve_rcwa is not None:
                solve_rcwa(solar_cell, options)
            else:
                raise RCWASolverError("RCWA optical solver not available!!")
        else:
            raise ValueError(
                'ERROR in "solar_cell_solver":\n\tOptics solver method must be None, "external", "BL", "TMM" or "RCWA".')

    else:
        print('Already calculated reflection, transmission and absorption profile - not recalculating. '
              'Set recalculate_absorption to True in the options if you want absorption to be calculated again.')


def solve_iv(solar_cell, options):
    """ Calculates the IV at a given voltage range, providing the IVs of the individual junctions in addition to the total IV

    :param solar_cell: A solar_cell object
    :param options: Options for the solvers
    :return: None
    """
    solve_optics(solar_cell, options)

    print('Solving IV of the junctions...')
    for j in solar_cell.junction_indices:

        if solar_cell[j].kind == 'PDD':
            PDD.iv_pdd(solar_cell[j], options)
        elif solar_cell[j].kind == 'DA':
            ASC.iv_depletion(solar_cell[j], options)
        elif solar_cell[j].kind == '2D':
            ASC.iv_2diode(solar_cell[j], options)
        elif solar_cell[j].kind == 'DB':
            ASC.iv_detailed_balance(solar_cell[j], options)
        else:
            raise ValueError(
                'ERROR in "solar_cell_solver":\n\tJunction {} has an invalid "type". It must be "PDD", "DA", "2D" or "DB".'.format(
                    j))

    print('Solving IV of the tunnel junctions...')
    for j in solar_cell.tunnel_indices:

        if solar_cell[j].kind == 'resistive':
            # The tunnel junction is modeled as a simple resistor
            ASC.resistive_tunnel_junction(solar_cell[j], options)
        elif solar_cell[j].kind == 'parametric':
            # The tunnel junction is modeled using a simple parametric model
            ASC.parametric_tunnel_junction(solar_cell[j], options)
        elif solar_cell[j].kind == 'external':
            # The tunnel junction is modeled using a simple parametric model
            ASC.external_tunnel_junction(solar_cell[j], options)
        elif solar_cell[j].kind == 'analytic':
            print('Sorry, the analytical tunnel junction model is not implemented, yet.')
        else:
            raise ValueError(
                'ERROR in "solar_cell_solver":\n\tTunnel junction {} has an invalid "type". It must be "parametric", "analytic", "external" or "resistive".'.format(
                    j))

    print('Solving IV of the total solar cell...')
    ASC.iv_multijunction(solar_cell, options)


def solve_qe(solar_cell, options):
    """ Calculates the QE of all the junctions

    :param solar_cell: A solar_cell object
    :param options: Options for the solvers
    :return: None
    """

    solve_optics(solar_cell, options)

    print('Solving QE of the solar cell...')
    for j in solar_cell.junction_indices:
        if solar_cell[j].kind == 'PDD':
            PDD.qe_pdd(solar_cell[j], options)
        elif solar_cell[j].kind == 'DA':
            ASC.qe_depletion(solar_cell[j], options)
        elif solar_cell[j].kind == '2D':
            # We solve this case as if it were DB. Therefore, to work it needs the same inputs in the Junction object
            wl = options.wavelength
            ASC.qe_detailed_balance(solar_cell[j], wl)
        elif solar_cell[j].kind == 'DB':
            wl = options.wavelength
            ASC.qe_detailed_balance(solar_cell[j], wl)
        else:
            raise ValueError(
                'ERROR in "solar_cell_solver":\n\tJunction {} has an invalid "type". It must be "PDD", "DA", "2D" or "DB".'.format(
                    j))


def solve_equilibrium(solar_cell, options):
    """ Uses the PDD solver to calculate the properties of all the all the junctions under equilibrium

    :param solar_cell: A solar_cell object
    :param options: Options for the solvers
    :return: None
    """
    for j in solar_cell.junction_indices:

        if solar_cell[j].kind == 'PDD':
            PDD.equilibrium_pdd(solar_cell[j], options)
        else:
            print('WARNING: Only PDD junctions can be solved in "equilibrium".')


def solve_short_circuit(solar_cell, options):
    """ Uses the PDD solver to calculate the properties of all the all the junctions under short circuit

    :param solar_cell: A solar_cell object
    :param options: Options for the solvers
    :return: None
    """

    solve_optics(solar_cell, options)

    for j in solar_cell.junction_indices:

        if solar_cell[j].kind == 'PDD':
            PDD.short_circuit_pdd(solar_cell[j], options)
        else:
            print('WARNING: Only PDD junctions can be solved in "short_circuit".')


def prepare_solar_cell(solar_cell, options):
    """ This function scans all the layers and junctions of the cell, calculating the relative position of each of them with respect the front surface (offset).
    This information will later be use by the optical calculators, for example. It also processes the 'position' option, which determines the spacing used if the
    solver is going to calculate depth-dependent absorption.

    :param solar_cell: A solar_cell object
    :param options: an options (State) object with user/default options
    :return: None
    """
    offset = 0
    layer_widths = []
    for j, layer_object in enumerate(solar_cell):

        # Independent layers, for example in a AR coating
        if type(layer_object) is Layer:
            layer_widths.append(layer_object.width)

        # Each Tunnel junctions can also have some layers with a given thickness.
        elif type(layer_object) is TunnelJunction:
            junction_width = 0
            for i, layer in enumerate(layer_object):
                junction_width += layer.width
                layer_widths.append(layer.width)
            solar_cell[j].width = junction_width

        # For each junction, and layer within the junction, we get the layer width.
        elif type(layer_object) is Junction:

            try:
                kind = solar_cell[j].kind
            except AttributeError as err:
                print('ERROR preparing the solar cell: Junction {} has no kind!'.format(j))
                raise err

            # This junctions will not, typically, have a width
            if kind in ['2D', 'DB']:
                layer_widths.append(1e-6)
                # 2D and DB junctions do not often have a width (or need it) so we set an arbitrary width
                if not hasattr(layer_object, 'width'):
                    solar_cell[j].width = 1e-6  # 1 µm

            else:
                junction_width = 0
                for i, layer in enumerate(layer_object):
                    layer_widths.append(layer.width)
                    junction_width += layer.width
                solar_cell[j].width = junction_width

        solar_cell[j].offset = offset
        offset += solar_cell[j].width

    solar_cell.width = offset

    process_position(solar_cell, options, layer_widths)



def process_position(solar_cell, options, layer_widths):
    """
    To control the depth spacing, the user can pass:
        - a vector which specifies each position (in m) at which the depth should be calculated
        - a single number which specifies the spacing (in m) to generate the position vector, e.g. 1e-9 for 1 nm spacing
        - a list of numbers which specify the spacing (in m) to be used in each layer. This list can have EITHER the length
        of the number of individual layers + the number of junctions in the cell object, OR the length of the total number of individual layers including layers inside junctions.

    :param solar_cell: a SolarCell object
    :param options: aan options (State) object with user/default options
    :param layer_widths: list of widths of the individual layers in the stack, treating the layers within junctions as individual layers
    :return: None
    """

    if options.position is None:
        options.position = [max(1e-10, width/5000) for width in layer_widths]

        layer_offsets = np.insert(np.cumsum(layer_widths), 0, 0)
        options.position = np.hstack([np.arange(layer_offsets[j],
                                                layer_offsets[j] + layer_width,
                                                options.position[j]) for j, layer_width
                                      in enumerate(layer_widths)])

    elif isinstance(options.position, int) or isinstance(options.position, float):
        options.position = np.arange(0, solar_cell.width, options.position)

    elif isinstance(options.position, list) or isinstance(options.position, np.ndarray):
        if len(options.position) == 1:
            options.position = np.arange(0, solar_cell.width, options.position[0])

        if len(options.position) == len(solar_cell):
            options.position = np.hstack([np.arange(layer_object.offset, layer_object.offset + layer_object.width, options.position[j]) for j, layer_object in enumerate(solar_cell)])

        elif len(options.position) == len(layer_widths):
            layer_offsets = np.insert(np.cumsum(layer_widths), 0, 0)
            options.position = np.hstack([np.arange(layer_offsets[j], layer_offsets[j] + layer_width, options.position[j]) for j, layer_width in enumerate(layer_widths)])
