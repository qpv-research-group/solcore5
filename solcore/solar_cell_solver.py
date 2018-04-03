import numpy as np

import solcore.analytic_solar_cells as ASC
from solcore.light_source import LightSource
from solcore.state import State
from solcore.optics import solve_beer_lambert, solve_tmm, solve_rcwa, rcwa_options, solve_external_optics
from solcore.structure import Layer, Junction, TunnelJunction

try:
    import solcore.poisson_drift_diffusion as PDD

    a = PDD.pdd_options
except AttributeError:
    PDD.pdd_options = {}

default_options = State()


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
default_options.position = np.linspace(0, 4660, 4661)
default_options.radiative_coupling = False

# Optics control
default_options.optics_method = 'BL'

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

    prepare_solar_cell(solar_cell)
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
    if options.optics_method is None:
        print('Warning: Not solving the optics of the solar cell.')
    elif options.optics_method == 'external':
        solve_external_optics(solar_cell, options)
    elif options.optics_method == 'BL':
        solve_beer_lambert(solar_cell, options)
    elif options.optics_method == 'TMM':
        solve_tmm(solar_cell, options)
    elif options.optics_method == 'RCWA':
        solve_rcwa(solar_cell, options)
    else:
        raise ValueError(
            'ERROR in "solar_cell_solver":\n\tOptics solver method must be None, "external", "BL", "TMM" or "RCWA".')


def solve_iv(solar_cell, options):
    """ Calculates the IV at a given voltage range, providing the IVs of the individual junctions in addition to the total IV

    :param solar_cell: A solar_cell object
    :param options: Options for the solvers
    :return: None
    """
    solve_optics(solar_cell, options)

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

    ASC.iv_multijunction(solar_cell, options)


def solve_qe(solar_cell, options):
    """ Calculates the QE of all the junctions

    :param solar_cell: A solar_cell object
    :param options: Options for the solvers
    :return: None
    """

    solve_optics(solar_cell, options)

    for j in solar_cell.junction_indices:
        if solar_cell[j].kind is 'PDD':
            PDD.qe_pdd(solar_cell[j], options)
        elif solar_cell[j].kind is 'DA':
            ASC.qe_depletion(solar_cell[j], options)
        elif solar_cell[j].kind is '2D':
            # We solve this case as if it were DB. Therefore, to work it needs the same inputs in the Junction object
            wl = options.wavelength
            ASC.qe_detailed_balance(solar_cell[j], wl)
        elif solar_cell[j].kind is 'DB':
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

        if solar_cell[j].kind is 'PDD':
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

        if solar_cell[j].kind is 'PDD':
            PDD.short_circuit_pdd(solar_cell[j], options)
        else:
            print('WARNING: Only PDD junctions can be solved in "short_circuit".')


def prepare_solar_cell(solar_cell):
    """ This function scans all the layers and junctions of the cell, calculating the relative possition of each of them with respect the front surface (offset). This information will later be use by the optical calculators, for example.

    :param solar_cell: A solar_cell object
    :return: None
    """
    offset = 0
    for j, layer_object in enumerate(solar_cell):

        # Independent layers, for example in a AR coating
        if type(layer_object) is Layer:
            pass

        # Each Tunnel junctions can also have some layers with a given thickness.
        elif type(layer_object) is TunnelJunction:
            junction_width = 0
            for i, layer in enumerate(layer_object):
                junction_width += layer.width
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
                # 2D and DB junctions do not often have a width (or need it) so we set an arbitrary width
                if not hasattr(layer_object, 'width'):
                    solar_cell[j].width = 1e-6  # 1 µm

            else:
                junction_width = 0
                for i, layer in enumerate(layer_object):
                    junction_width += layer.width
                solar_cell[j].width = junction_width

        solar_cell[j].offset = offset
        offset += solar_cell[j].width

    solar_cell.width = offset
