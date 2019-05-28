from solcore.structure import Layer, Junction, TunnelJunction
from solcore.absorption_calculator import calculate_rat_rcwa, calculate_absorption_profile_rcwa

import numpy as np
import types
from solcore.state import State

rcwa_options = State()
rcwa_options.size = [500, 500]
rcwa_options.orders = 4
rcwa_options.theta = 0
rcwa_options.phi = 0
rcwa_options.pol = 'u'


def solve_rcwa(solar_cell, options):
    """ Calculates the reflection, transmission and absorption of a solar cell object using the rigorous coupled-wave analysis solver.

    :param solar_cell: A solar_cell object
    :param options: Options for the solver
    :return: None
    """
    wl = options.wavelength
    solar_cell.wavelength = options.wavelength

    # We include the shadowing losses
    initial = (1 - solar_cell.shading) if hasattr(solar_cell, 'shading') else 1

    # Now we calculate the absorbed and transmitted light. We first get all the relevant parameters from the objects
    all_layers = []
    for j, layer_object in enumerate(solar_cell):

        # Attenuation due to absorption in the AR coatings or any layer in the front that is not part of the junction
        if type(layer_object) is Layer:
            all_layers.append(layer_object)

        # For each junction, and layer within the junction, we get the absorption coefficient and the layer width.
        elif type(layer_object) in [TunnelJunction, Junction]:
            for i, layer in enumerate(layer_object):
                all_layers.append(layer)

    # With all the information, we create the optical stack
    stack = all_layers

    position = options.position * 1e9
    # angle_theta = options.angle_theta if 'angle_theta' in options.keys() else 0
    # angle_phi = options.angle_phi if 'angle_phi' in options.keys() else 0
    # pol = options.pol if 'pol' in options.keys() else 'u'
    # size = options.size if 'size' in options.keys() else [500, 500]
    # orders = options.orders if 'orders' in options.keys() else 4
    substrate = solar_cell.substrate

    print('Calculating RAT...')
    RAT = calculate_rat_rcwa(stack, options.size, options.orders, wl * 1e9, theta=options.theta,
                             phi=options.phi, pol=options.pol, substrate=substrate)

    print('Calculating absorption profile...')
    out = calculate_absorption_profile_rcwa(stack, options.size, options.orders, wl * 1e9, RAT,
                                            dist=position, theta=options.theta, phi=options.phi, pol=options.pol, substrate=substrate)

    # With all this information, we are ready to calculate the differential absorption function
    diff_absorption, all_absorbed = calculate_absorption_rcwa(out, initial)

    # Each building block (layer or junction) needs to have access to the absorbed light in its region.
    # We update each object with that information.
    for j in range(len(solar_cell)):
        solar_cell[j].diff_absorption = diff_absorption
        solar_cell[j].absorbed = types.MethodType(absorbed, solar_cell[j])

        layer_positions = options.position[(options.position >= solar_cell[j].offset) & (
                options.position < solar_cell[j].offset + solar_cell[j].width)]
        layer_positions = layer_positions - np.min(layer_positions)
        solar_cell[j].layer_absorption = np.trapz(solar_cell[j].absorbed(layer_positions), layer_positions, axis=0)

    solar_cell.reflected = RAT['R'] * initial
    solar_cell.absorbed = sum([solar_cell[x].layer_absorption for x in np.arange(len(solar_cell))])
    solar_cell.transmitted = initial - solar_cell.reflected - solar_cell.absorbed


def absorbed(self, z):
    out = self.diff_absorption(self.offset + z) * (z < self.width)
    return out.T


def calculate_absorption_rcwa(tmm_out, initial=1):
    all_z = tmm_out['position'] * 1e-9
    all_abs = initial * tmm_out['absorption'] / 1e-9

    def diff_absorption(z):
        idx = all_z.searchsorted(z)
        idx = np.where(idx <= len(all_z) - 2, idx, len(all_z) - 2)
        try:
            z1 = all_z[idx]
            z2 = all_z[idx + 1]

            f = (z - z1) / (z2 - z1)

            out = f * all_abs[:, idx] + (1 - f) * all_abs[:, idx + 1]

        except IndexError:
            out = all_abs[:, idx]

        return out

    all_absorbed = np.trapz(diff_absorption(all_z), all_z)

    return diff_absorption, all_absorbed
