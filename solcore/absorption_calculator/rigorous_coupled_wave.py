import numpy as np
import tmm
from solcore.structure import Layer
from solcore.absorption_calculator import OptiStack

try:
    import S4
except ModuleNotFoundError:
    raise


def calculate_rat_rcwa(structure, size, orders, wavelength, theta=0, phi=0, pol='u', substrate=None):
    """ Calculates the reflected, absorbed and transmitted intensity of the structure for the wavelengths and angles
    defined using an RCWA method implemented using the S4 package.

    :param structure: A solcore Structure object with layers and materials or a OptiStack object.
    :param size: list with 2 entries, size of the unit cell (right now, can only be rectangular
    :param orders: number of orders to retain in the RCWA calculations.
    :param wavelength: Wavelengths (in nm) in which calculate the data.
    :param theta: polar incidence angle (in degrees) of the incident light. Default: 0 (normal incidence)
    :param phi: azimuthal incidence angle in degrees. Default: 0
    :param pol: Polarisation of the light: 's', 'p' or 'u'. Default: 'u' (unpolarised).
    :param substrate: semi-infinite transmission medium
    :return: A dictionary with the R, A and T at the specified wavelengths and angle.
    """

    num_wl = len(wavelength)

    # write a separate function that makes the OptiStack structure into an S4 object, defined materials etc.
    S, stack_OS, shape_mats_OS = initialise_S(structure, size, orders, substrate)

    output = {'R': np.zeros(num_wl), 'A': np.zeros(num_wl), 'T': np.zeros(num_wl),
              'A_layer': np.zeros((num_wl, len(stack_OS.get_widths())-2))}


    if pol in 'sp':
        if pol == 's':
            s = 1
            p = 0
        elif pol == 'p':
            s = 0
            p = 1

        S.SetExcitationPlanewave((theta, phi), s, p, 0)
        for i, wl in enumerate(wavelength):  # set the material values and indices in here
            update_epsilon(S, stack_OS, shape_mats_OS, wl)
            S.SetFrequency(1 / wl)
            out = rcwa_rat(S, len(stack_OS.get_widths()))
            output['R'][i] = out['R']
            output['A'][i] = 1 - out['R'] - out['T']
            output['T'][i] = out['T']
            output['A_layer'][i] = rcwa_absorption_per_layer(S, len(stack_OS.get_widths()))

    else:
        for i, wl in enumerate(wavelength):  # set the material values and indices in here
            S.SetFrequency(1 / wl)
            update_epsilon(S, stack_OS, shape_mats_OS, wl)
            S.SetExcitationPlanewave((theta, phi), 0, 1, 0)  # p-polarization
            out_p = rcwa_rat(S, len(stack_OS.get_widths()))
            S.SetExcitationPlanewave((theta, phi), 1, 0, 0)  # s-polarization
            out_s = rcwa_rat(S, len(stack_OS.get_widths()))

            output['R'][i] = 0.5 * (out_p['R'] + out_s['R'])  # average
            output['T'][i] = 0.5 * (out_p['T'] + out_s['T'])
            output['A'][i] = 1 - output['R'][i] - output['T'][i]
            # output['all_p'].append(out_p['power_entering_list'])
            # output['all_s'].append(out_s['power_entering_list'])
            output['A_layer'][i] = rcwa_absorption_per_layer(S, len(stack_OS.get_widths()))
    return output


def rcwa_rat(S, n_layers):
    below = 'layer_' + str(n_layers)  # identify which layer is the transmission medium
    R = 1 - sum(
        S.GetPowerFlux('layer_2'))  # GetPowerFlux gives forward & backward Poynting vector, so sum to get power flux
    # layer_2 is the top layer of the structure (layer_1 is incidence medium)
    T = sum(S.GetPowerFlux(below))
    return {'R': np.real(R), 'T': np.real(T)}


def initialise_S(stack, size, orders, substrate=None):
    S = S4.New(size, orders)
    S.SetOptions(  # these are the default
        LatticeTruncation='Circular',
        PolarizationDecomposition=False,
        PolarizationBasis='Default'#,
        #WeismannFormulation = False
    )
    geom_list = [layer.geometry for layer in stack]
    geom_list.insert(0, {})  # incidence medium
    geom_list.append({})  # transmission medium

    ## Materials for the shapes need to be defined before you can do .SetRegion
    shape_mats = necessary_materials(geom_list)
    Layers = []
    for x in shape_mats:
        Layers.append(Layer(0, x))

    shape_mats_OS = OptiStack(Layers)

    for i1 in range(len(shape_mats)):  # create the materials needed for all the shapes in S4
        S.SetMaterial('shape_mat_' + str(i1 + 1), 1)

    ## Make the layers
    stack_OS = OptiStack(stack, no_back_reflection=False, substrate=substrate)
    widths = stack_OS.get_widths()

    for i1 in range(len(widths)):  # create 'dummy' materials for base layers including incidence and transmission media
        S.SetMaterial('layer_' + str(i1 + 1), 1)  # This is not strictly necessary but it means S.SetExcitationPlanewave
        # can be done outside the wavelength loop in calculate_rat_rcwa
    for i1 in range(len(widths)):  # set base layers
        layer_name = 'layer_' + str(i1 + 1)
        if widths[i1] == float('Inf'):
            S.AddLayer(layer_name, 0, layer_name)  # Solcore4 has incidence and transmission media widths set to Inf;
            # in S4 they have zero width
        else:
            S.AddLayer(layer_name, widths[i1], layer_name)  # keep base unit m, not nm

        geometry = geom_list[i1]

        if bool(geometry):
            for shape in geometry:
                mat_name = 'shape_mat_' + str(shape_mats.index(shape['mat']) + 1)
                if shape['type'] == 'circle':
                    S.SetRegionCircle(layer_name, mat_name, shape['center'], shape['radius'])
                elif shape['type'] == 'ellipse':
                    S.SetRegionEllipse(layer_name, mat_name, shape['center'], shape['angle'], shape['halfwidths'])
                elif shape['type'] == 'rectangle':
                    S.SetRegionRectangle(layer_name, mat_name, shape['center'], shape['angle'], shape['halfwidths'])
                elif shape['type'] == 'polygon':
                    S.SetRegionPolygon(layer_name, mat_name, shape['center'], shape['angle'], shape['vertices'])

    return S, stack_OS, shape_mats_OS


def necessary_materials(geom_list):
    shape_mats = []
    for i1, geom in enumerate(geom_list):
        if bool(geom):
            shape_mats.append([x['mat'] for x in geom])
    return list(set([val for sublist in shape_mats for val in sublist]))


def update_epsilon(S, stack_OS, shape_mats_OS, wl):
    for i1 in range(len(stack_OS.get_widths())):
        S.SetMaterial('layer_' + str(i1 + 1), stack_OS.get_indices(wl)[i1] ** 2)
    for i1 in range(len(shape_mats_OS.widths)):  # initialise the materials needed for all the shapes in S4
        S.SetMaterial('shape_mat_' + str(i1 + 1), shape_mats_OS.get_indices(wl)[i1 + 1] ** 2)

    return S


def calculate_absorption_profile_rcwa(structure, size, orders, wavelength, rat_output,
                                      z_limit=None, steps_size=2, dist=None, theta=0, phi=0, pol='u', substrate=None):
    """ It calculates the absorbed energy density within the material. From the documentation:

    'In principle this has units of [power]/[volume], but we can express it as a multiple of incoming light power
    density on the material, which has units [power]/[area], so that absorbed energy density has units of 1/[length].'

    Integrating this absorption profile in the whole stack gives the same result that the absorption obtained with
    calculate_rat as long as the spacial mesh (controlled by steps_thinest_layer) is fine enough. If the structure is
    very thick and the mesh not thin enough, the calculation might diverege at short wavelengths.

    For now, it only works for normal incident, coherent light.

    :param structure: A solcore structure with layers and materials.
    :param size: list with 2 entries, size of the unit cell (right now, can only be rectangular
    :param orders: number of orders to retain in the RCWA calculations.
    :param wavelength: Wavelengths (in nm) in which calculate the data.
    :param rat_output: output from calculate_rat_rcwa
    :param z_limit: Maximum value in the z direction at which to calculate depth-dependent absorption (nm)
    :param steps_size: if the dist is not specified, the step size in nm to use in the depth-dependent calculation
    :param dist: the positions (in nm) at which to calculate depth-dependent absorption
    :param theta: polar incidence angle (in degrees) of the incident light. Default: 0 (normal incidence)
    :param phi: azimuthal incidence angle in degrees. Default: 0
    :param pol: Polarisation of the light: 's', 'p' or 'u'. Default: 'u' (unpolarised).
    :param substrate: semi-infinite transmission medium
    :return: A dictionary containing the positions (in nm) and a 2D array with the absorption in the structure as a
    function of the position and the wavelength.
    """

    num_wl = len(wavelength)

    if dist is None:
        if z_limit is None:
            stack = OptiStack(structure)
            z_limit = np.sum(np.array(stack.widths))
        dist = np.arange(0, z_limit, steps_size)

    output = {'position': dist, 'absorption': np.zeros((num_wl, len(dist)))}

    S, stack_OS, shape_mats_OS = initialise_S(structure, size, orders, substrate)

    if pol in 'sp':
        if pol == 's':
            s = 1
            p = 0
        elif pol == 'p':
            s = 0
            p = 1

        S.SetExcitationPlanewave((theta, phi), s, p, 0)
        for i, wl in enumerate(wavelength):
            update_epsilon(S, stack_OS, shape_mats_OS, wl)
            S.SetFrequency(1 / wl)
            A = rat_output['A'][i]
            for j, d in enumerate(dist):
                layer, d_in_layer = tmm.find_in_structure_with_inf(stack_OS.get_widths(),
                                                                   d)  # don't need to change this
                layer_name = 'layer_' + str(layer + 1)  # layer_1 is air above so need to add 1
                data = rcwa_position_resolved(S, layer_name, d_in_layer, A)
                output['absorption'][i, j] = data

    else:
        for i, wl in enumerate(wavelength):  # set the material values and indices in here
            print(i)
            update_epsilon(S, stack_OS, shape_mats_OS, wl)
            S.SetFrequency(1 / wl)
            A = rat_output['A'][i]

            for j, d in enumerate(dist):
                layer, d_in_layer = tmm.find_in_structure_with_inf(stack_OS.get_widths(),
                                                                   d)  # don't need to change this
                layer_name = 'layer_' + str(layer + 1)  # layer_1 is air above so need to add 1
                S.SetExcitationPlanewave((theta, phi), 0, 1, 0)  # p-polarization
                data_p = rcwa_position_resolved(S, layer_name, d_in_layer, A)
                S.SetExcitationPlanewave((theta, phi), 1, 0, 0)  # p-polarization
                data_s = rcwa_position_resolved(S, layer_name, d_in_layer, A)
                output['absorption'][i, j] = 0.5 * (data_p + data_s)

    return output


def rcwa_position_resolved(S, layer, depth, A):
    if A > 0:
        delta = 1e-9
        power_difference = np.real(
            sum(S.GetPowerFlux(layer, depth - delta)) - sum(S.GetPowerFlux(layer, depth + delta)))
        return power_difference / (2 * delta)  # absorbed energy density normalised to total absorption
    else:
        return 0

def rcwa_absorption_per_layer(S, n_layers):
    # layer 1 is incidence medium, layer n is the transmission medium
    A = np.empty(n_layers-2)
    for i1, layer in enumerate(np.arange(n_layers-2)+2):
        A[i1] = np.real(sum(S.GetPowerFlux('layer_' + str(layer))) - sum(S.GetPowerFlux('layer_' + str(layer+1))))

    A = [x if x > 0 else 0 for x in A]

    return A