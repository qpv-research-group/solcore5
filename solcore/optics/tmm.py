from solcore.structure import Layer, Junction, TunnelJunction
import solcore.analytic_solar_cells as ASC
from solcore.absorption_calculator import calculate_rat, OptiStack, calculate_absorption_profile
from solcore.optics.beer_lambert import calculate_absorption_beer_lambert
from solcore import si

import numpy as np
import types
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


def solve_tmm(solar_cell, options):
    """ Calculates the reflection, transmission and absorption of a solar cell object using the transfer matrix method.
    Internally, it creates an OptiStack and then it calculates the optical properties of the whole structure.
    A substrate can be specified in the SolarCell object, which is treated as a semi-infinite transmission medium.
    Shading can also be specified (as a fraction).
    Relevant options are 'wl' (the wavelengths, in m), the incidence angle 'theta' (in degrees), the polarization 'pol' ('s',
    'p' or 'u'), 'position' (locations in m at which depth-dependent absorption is calculated), 'no_back_reflexion' and 'BL_correction'.
    'no_back_reflexion' sets whether reflections from the back surface are suppressed (if set to True, the default),
    or taken into account (if set to False).
    If 'BL_correction' is set to True, the absorption at very high attenuation
    (absorption coefficient * width > 150 so attenuation is 1e65 times) is calculated using the Beer-Lambert law.
    This is to avoid absorption going to infinity at short wavelengths. If 'BL_correction' is set to True, thick layers
    (thickness > 10*maximum wavelength) are also treated incoherently using the Beer-Lambert law, to avoid the calculation of unphysical
    oscillations.


    :param solar_cell: A SolarCell object
    :param options: Options for the solver
    :return: None
    """
    wl = options.wavelength
    BL_correction = options.BL_correction if 'BL_correction' in options.keys() else True
    theta = options.theta if 'theta' in options.keys() else 0 # angle IN DEGREES
    pol = options.pol if 'pol' in options.keys() else 'u' # angle IN DEGREES

    # We include the shadowing losses
    initial = (1 - solar_cell.shading) if hasattr(solar_cell, 'shading') else 1

    # Now we calculate the absorbed and transmitted light. We first get all the relevant parameters from the objects
    all_layers = []

    # the next three are for Beer-Lambert correction only
    widths = []
    alphas = []
    n_layers_junction = []
    for j, layer_object in enumerate(solar_cell):

        # Attenuation due to absorption in the AR coatings or any layer in the front that is not part of the junction
        if type(layer_object) is Layer:
            all_layers.append(layer_object)
            widths.append(layer_object.width)
            alphas.append(layer_object.material.alpha(wl))
            n_layers_junction.append(1)

        # For each junction, and layer within the junction, we get the absorption coefficient and the layer width.
        elif type(layer_object) in [TunnelJunction, Junction]:
            n_layers_junction.append(len(layer_object))
            for i, layer in enumerate(layer_object):
                all_layers.append(layer)
                widths.append(layer.width)
                alphas.append(layer.material.alpha(wl))

    # With all the information, we create the optical stack
    no_back_reflexion = options.no_back_reflexion if 'no_back_reflexion' in options.keys() else True
    attn = np.multiply(np.array(widths), np.array(alphas).T).T
    #byBL = (attn/np.cos(theta*np.pi/180)) > 150
    byBL = attn > 150
    BL_from = len(all_layers)
    if BL_correction:
        print('Using Beer-Lambert absorption profile for optically thick layers (exp(-OD) < 1e-65 at normal incidence)')
        solar_cell.byBL = byBL
        if any(widths > 10*np.max(wl)): # assume it's safe to ignore interference effects
            BL_from = np.where(np.array(widths) > si('10um'))[0][0]
            print('Ignoring layer ' + str(BL_from) + ' and all layers below it in TMM calculation')
            byBL[BL_from:,:] = True

    else:
        byBL[:] = False

    stack = OptiStack(all_layers[0:BL_from], no_back_reflexion=no_back_reflexion, substrate=solar_cell.substrate)
    full_stack = OptiStack(all_layers, no_back_reflexion=no_back_reflexion, substrate=solar_cell.substrate)

    if 'coherency_list' in options.keys():
        coherency_list = options.coherency_list
        coherent = False
        assert len(coherency_list) == full_stack.num_layers, \
            'Error: The coherency list must have as many elements (now {}) as the ' \
            'number of layers (now {}).'.format(len(coherency_list), full_stack.num_layers)
        coherency_list_trunc = coherency_list[0:BL_from]
    else:
        coherency_list = None
        coherency_list_trunc = None
        coherent = True

    position = options.position * 1e9
    profile_position = position[position < sum(stack.widths)]

    print('Calculating RAT...')
    RAT = calculate_rat(full_stack, wl * 1e9, angle=theta,
                        coherent=coherent, coherency_list=coherency_list, no_back_reflexion=no_back_reflexion,
                        pol=pol)

    print('Calculating absorption profile...')
    out = calculate_absorption_profile(stack, wl * 1e9, dist=profile_position,
                                       angle=theta, no_back_reflexion=no_back_reflexion,
                                       pol=pol, coherent=coherent,
                                       coherency_list=coherency_list_trunc)

    # With all this information, we are ready to calculate the differential absorption function
    diff_absorption, all_absorbed = calculate_absorption_tmm(out, initial)

    # Each building block (layer or junction) needs to have access to the absorbed light in its region.
    # We update each object with that information.
    previous_abs = 0
    layer = 0
    for j in range(len(solar_cell)):

        layer_positions = options.position[(options.position >= solar_cell[j].offset) & (
                options.position < solar_cell[j].offset + solar_cell[j].width)]
        layer_positions = layer_positions - np.min(layer_positions)
        #fraction = initial - RAT['R']*initial - previous_abs
        fraction = initial*(1-RAT['R']) - previous_abs
        diff_absorption_BL, transmitted_BL, all_absorbed_BL = calculate_absorption_beer_lambert(widths[layer:layer+n_layers_junction[j]],
                                                                                                alphas[layer:layer+n_layers_junction[j]],
                                                                                                fraction)
        solar_cell[j].diff_absorption_BL = diff_absorption_BL
        solar_cell[j].diff_absorption = diff_absorption
        # solar_cell[j].where_BL = byBL[(layer[0]-1)]
        solar_cell[j].where_BL = byBL[layer:(layer + n_layers_junction[j])]

        # print(solar_cell[j].where_BL)

        solar_cell[j].absorbed = types.MethodType(absorbed, solar_cell[j])
        solar_cell[j].layer_absorption = np.trapz(solar_cell[j].absorbed(layer_positions), layer_positions, axis=0)
        previous_abs = previous_abs + solar_cell[j].layer_absorption
        layer = layer + n_layers_junction[j]

    solar_cell.reflected = RAT['R'] * initial
    solar_cell.absorbed = sum([solar_cell[x].layer_absorption for x in np.arange(len(solar_cell))])
    solar_cell.transmitted = initial - solar_cell.reflected - solar_cell.absorbed


def absorbed(self, z):
    out = self.diff_absorption(self.offset + z) * (z < self.width)

    if type(self) is Layer:
        replace = np.broadcast_to(self.where_BL[0], (len(z), len(self.where_BL[0]))).T

    else:
        replace = np.broadcast_to(np.any(self.where_BL, axis=0), (len(z), len(self.where_BL[0]))).T

    out_BL = (self.diff_absorption_BL(z).T) * (z < self.width)
    out[replace] = 0
    out_BL[np.logical_not(replace)] = 0
    out = out + out_BL

    return out.T


def calculate_absorption_tmm(tmm_out, initial=1):
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
