from solcore.structure import Layer, Junction, TunnelJunction

import numpy as np
import types


def solve_external_optics(solar_cell, options):
    """ Calculates the reflection, transmission and absorption of a solar cell from external data, preparing the structure for further calculations. The external data that must be provided is:

    external_reflected: a function that provides the fraction of reflected light at the specified wavelengths (in m)
    external_absorbed: a function that provides an array with the differential absorption at a depth Z (in m) at the specified wavelengths.

    Note that both are functions of a single variable - wavelength in the first case and position in the second - and that the latter provides as output an array at the wavelengths indicated in the options.

    :param solar_cell:
    :param options:
    :return:
    """

    wl = options.wavelength

    # We include the shadowing losses
    initial = (1 - solar_cell.shading) if hasattr(solar_cell, 'shading') else 1

    # Now we calculate the absorbed and transmitted light. We first get all the relevant parameters from the objects
    widths = []
    alphas = []
    offset = 0
    for j, layer_object in enumerate(solar_cell):

        # Attenuation due to absorption in the AR coatings or any layer in the front that is not part of the junction
        if type(layer_object) is Layer:
            widths.append(layer_object.width)

        # Each Tunnel junctions can also have some layers absorbing light.
        elif type(layer_object) is TunnelJunction:
            junction_width = 0
            for i, layer in enumerate(layer_object):
                junction_width += layer.width
                widths.append(layer.width)

            solar_cell[j].width = junction_width

        # For each junction, and layer within the junction, we get the absorption coeficient and the layer width.
        elif type(layer_object) is Junction:

            try:
                kind = solar_cell[j].kind
            except AttributeError as err:
                print('ERROR in the external optics calculator: Junction {} has no kind!'.format(j))
                raise err

            if kind in ['2D', 'DB']:
                try:
                    w = layer_object.width
                    widths.append(w)
                except AttributeError as err:
                    print(
                        'ERROR in the external optics calculator: Junction {} of kind {} has no width defined!'.format(
                            j, kind))
                    raise err

            else:
                junction_width = 0

                try:
                    for i, layer in enumerate(layer_object):
                        junction_width += layer.width
                        widths.append(layer.width)

                    solar_cell[j].width = junction_width

                except TypeError as err:
                    print('ERROR in the external optics calculator: Junction {} of kind {} has no layers!'.format(j,
                                                                                                                  kind))
                    raise err

        solar_cell[j].offset = offset
        offset += layer_object.width

    # We try to get the external attributes
    try:
        solar_cell.reflected = solar_cell.external_reflected * initial
        diff_absorption = solar_cell.external_absorbed
    except AttributeError as err:
        raise err

    # We calculate the total amount of light absorbed in the solar cell, integrating over its whole thickness with a step of 1 nm
    all_z = np.arange(0, offset, 1e-9)
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
