from solcore.structure import Layer, Junction, TunnelJunction

import numpy as np
import types


def solve_external_optics(solar_cell, options):
    """ Calculates the reflection, transmission and absorption of a solar cell from external data, preparing the structure for further calculations. The external data that must be provided in the solar cell definition is:

    external_reflected: a function that provides the fraction of reflected light at the specified wavelengths (in m)
    external_absorbed: a function that provides an array with the differential absorption at a depth Z (in m) at the specified wavelengths.

    Note that both are functions of a single variable - wavelength in the first case and position in the second - and that the latter provides as output an array at the wavelengths indicated in the options.

    :param solar_cell:
    :param options:
    :return:
    """
    solar_cell.wavelength = options.wavelength

    # We include the shadowing losses
    initial = (1 - solar_cell.shading) if hasattr(solar_cell, 'shading') else 1

    # We try to get the external attributes
    try:
        solar_cell.reflected = solar_cell.external_reflected * initial
        diff_absorption = solar_cell.external_absorbed
    except AttributeError as err:
        raise err

    # We calculate the total amount of light absorbed in the solar cell, integrating over its whole thickness with a step of 1 nm
    all_z = options.position
    all_absorbed = np.trapz(diff_absorption(all_z), all_z)

    # Each building block (layer or junction) needs to have access to the absorbed light in its region.
    # We update each object with that information.
    for j in range(len(solar_cell)):
        solar_cell[j].diff_absorption = diff_absorption
        solar_cell[j].absorbed = types.MethodType(absorbed, solar_cell[j])

    solar_cell.transmitted = (1 - solar_cell.external_reflected - all_absorbed) * initial
    solar_cell.absorbed = all_absorbed * initial


def absorbed(self, z):
    out = self.diff_absorption(self.offset + z) * (z < self.width)
    return out.T
