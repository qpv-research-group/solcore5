from solcore.structure import Layer, Junction, TunnelJunction
from solcore.absorption_calculator import calculate_rat, OptiStack, calculate_absorption_profile

import numpy as np
import types
from warnings import warn


def solve_tmm(solar_cell, options):
    """ Calculates the reflection, transmission and absorption of a solar cell object using the transfer matrix method.
    Internally, it creates an OptiStack and then it calculates the optical properties of the whole structure.
    A substrate can be specified in the SolarCell object, which is treated as a semi-infinite transmission medium.
    Shading can also be specified (as a fraction).

    Relevant options are 'wl' (the wavelengths, in m), the incidence angle 'theta' (in degrees), the polarization 'pol' ('s',
    'p' or 'u'), 'position' (locations in m at which depth-dependent absorption is calculated), 'no_back_reflection' and 'BL_correction'.
    'no_back_reflection' sets whether reflections from the back surface are suppressed (if set to True, the default),
    or taken into account (if set to False).

    If 'BL_correction' is set to True, thick layers (thickness > 10*maximum wavelength) are treated incoherently using
    the Beer-Lambert law, to avoid the calculation of unphysical interference oscillations in the R/A/T spectra.

    A coherency_list option can be provided: this should have elements equal to the total number of layers (if a Junction
    contains multiple Layers, each should have its own entry in the coherency list). Each element is either 'c' for coherent
    treatment of that layer or 'i' for incoherent treatment.

    :param solar_cell: A SolarCell object
    :param options: Options for the solver
    :return: None
    """
    wl = options.wavelength
    BL_correction = options.BL_correction if 'BL_correction' in options.keys() else True
    theta = options.theta if 'theta' in options.keys() else 0 # angle IN DEGREES
    pol = options.pol if 'pol' in options.keys() else 'u'
    zero_threshold = options.get("zero_threshold", 1e-5)

    # We include the shadowing losses
    initial = (1 - solar_cell.shading) if hasattr(solar_cell, 'shading') else 1

    # Now we calculate the absorbed and transmitted light. We first get all the relevant parameters from the objects
    all_layers = []
    widths = []
    n_layers_junction = []

    for j, layer_object in enumerate(solar_cell):

        # Attenuation due to absorption in the AR coatings or any layer in the front that is not part of the junction
        if type(layer_object) is Layer:
            all_layers.append(layer_object)
            widths.append(layer_object.width)
            n_layers_junction.append(1)

        # For each junction, and layer within the junction, we get the absorption coefficient and the layer width.
        elif type(layer_object) in [TunnelJunction, Junction]:
            n_layers_junction.append(len(layer_object))
            for i, layer in enumerate(layer_object):
                all_layers.append(layer)
                widths.append(layer.width)

    no_back_reflection = options.no_back_reflection if 'no_back_reflection' in options.keys() else True

    # With all the information, we create the optical stack
    if 'no_back_reflexion' in options.keys():
        warn('The no_back_reflexion warning is deprecated. Use no_back_reflection instead.', FutureWarning)
        no_back_reflection = options.no_back_reflexion


    full_stack = OptiStack(all_layers, no_back_reflection=no_back_reflection, substrate=solar_cell.substrate, incidence=solar_cell.incidence)

    if 'coherency_list' in options.keys():
        coherency_list = options.coherency_list
        coherent = False
        assert len(coherency_list) == full_stack.num_layers, \
            'Error: The coherency list must have as many elements (now {}) as the ' \
            'number of layers (now {}).'.format(len(coherency_list), full_stack.num_layers)
    else:
        coherency_list = None
        coherent = True

    if BL_correction and any(widths > 10*np.max(wl)): # assume it's safe to ignore interference effects
            make_incoherent = np.where(np.array(widths) > 10*np.max(wl))[0]
            print('Treating layer(s) ' + str(make_incoherent).strip('[]') + ' incoherently')
            if not 'coherency_list' in options.keys():
                coherency_list = np.array(len(all_layers)*['c'])
                coherent = False
            else:
                coherency_list = np.array(coherency_list)
            coherency_list[make_incoherent] = 'i'
            coherency_list = coherency_list.tolist()

    position = options.position * 1e9
    profile_position = position[position < sum(full_stack.widths)]

    print('Calculating RAT...')
    RAT = calculate_rat(full_stack, wl * 1e9, angle=theta,
                        coherent=coherent, coherency_list=coherency_list, no_back_reflection=no_back_reflection,
                        pol=pol)

    print('Calculating absorption profile...')
    out = calculate_absorption_profile(full_stack, wl * 1e9, dist=profile_position,
                                       angle=theta, no_back_reflection=no_back_reflection,
                                       pol=pol, coherent=coherent,
                                       coherency_list=coherency_list, zero_threshold=zero_threshold, RAT_out=RAT)

    # With all this information, we are ready to calculate the differential absorption function
    diff_absorption, all_absorbed = calculate_absorption_tmm(out, initial)

    # Each building block (layer or junction) needs to have access to the absorbed light in its region.
    # We update each object with that information.
    layer = 0
    A_per_layer = np.array(RAT['A_per_layer'][1:-1]) # first entry is R, last entry is T
    for j in range(len(solar_cell)):

        solar_cell[j].layer_absorption = initial*np.sum(A_per_layer[layer:(layer+n_layers_junction[j])], axis = 0)

        solar_cell[j].diff_absorption = diff_absorption
        solar_cell[j].absorbed = types.MethodType(absorbed, solar_cell[j])

        layer = layer + n_layers_junction[j]

    solar_cell.reflected = RAT['R'] * initial
    solar_cell.absorbed = sum([solar_cell[x].layer_absorption for x in np.arange(len(solar_cell))])
    solar_cell.transmitted = initial - solar_cell.reflected - solar_cell.absorbed


def absorbed(self, z):
    out = self.diff_absorption(self.offset + z) * (z < self.width)

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
