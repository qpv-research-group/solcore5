from solcore.structure import Layer, Junction, TunnelJunction
import solcore.analytic_solar_cells as ASC
from solcore.absorption_calculator import calculate_rat, OptiStack, calculate_absorption_profile

from solcore.absorption_calculator.transfer_matrix import calculate_inc_absorption_profile

import numpy as np
import types
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


def solve_tmm(solar_cell, options):
    """ Calculates the reflection, transmission and absorption of a solar cell object using the transfer matrix method. Internally, it creates an OptiStack and then it calculates the optical properties of the whole structure.

    :param solar_cell: A solar_cell object
    :param options: Options for the solver
    :return: None
    """
    wl = options.wavelength

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
    no_back_reflexion = options.no_back_reflexion if 'no_back_reflexion' in options.keys() else True
    stack = OptiStack(all_layers, no_back_reflexion=no_back_reflexion)

    position = options.position * 1e9
    
    try:
        angle = options.angle
    except:
        angle = 0

    try:
        # Incoming and outgoing layers have to be stated explicitly as incoherent.
        coherency_list = ['i'] + options.coherency_list + (['i','i'] if no_back_reflexion else ['i'])
    except:
        coherency_list = None
    
    if coherency_list is None:
        if angle!=0:
            raise ValueError('Non normal incident is not implemented yet for coherent calculations')
        print('Calculating RAT...')
        RAT = calculate_rat(stack, wl * 1e9, coherent=True, no_back_reflexion=no_back_reflexion)
        print('Calculating absorption profile...')
        out = calculate_absorption_profile(stack, wl * 1e9, dist=position, no_back_reflexion=no_back_reflexion)
    
        # With all this information, we are ready to calculate the differential absorption function
        diff_absorption, all_absorbed = calculate_absorption_tmm(out)
    
        # Each building block (layer or junction) needs to have access to the absorbed light in its region.
        # We update each object with that information.
        for j in range(len(solar_cell)):
            solar_cell[j].diff_absorption = diff_absorption
            solar_cell[j].absorbed = types.MethodType(absorbed, solar_cell[j])
            
    else:
        print('Calculating RAT...')
        RAT = calculate_rat(stack,
                            wl * 1e9,
                            coherent=False,
                            coherency_list=coherency_list,
                            no_back_reflexion=no_back_reflexion,
                            angle=angle)
        
        print('Calculating absorption profile...')
        absorption_per_layer, absorption_position_resolved = \
            calculate_inc_absorption_profile(stack, 
                                             wl * 1e9, 
                                             coherency_list, 
                                             angle=angle, 
                                             no_back_reflexion=no_back_reflexion)
        
        # Each building block (layer or junction) needs to have access to the absorbed light in its region.
        # We update each object with that information.
        all_absorbed = absorption_per_layer.sum(axis=0)
        for j in range(len(solar_cell)):
            solar_cell[j].absorbed = absorption_position_resolved[j]
            solar_cell[j].absorbed_in_layer = absorption_per_layer[j]
 
    solar_cell.reflected = RAT['R'] * initial
    solar_cell.transmitted = (1 - RAT['R'] - all_absorbed) * initial
    solar_cell.absorbed = all_absorbed * initial

def absorbed(self, z):
    out = self.diff_absorption(self.offset + z) * (z < self.width)
    return out.T


def calculate_absorption_tmm(tmm_out):
    all_z = tmm_out['position'] * 1e-9
    all_abs = tmm_out['absorption'] / 1e-9

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
