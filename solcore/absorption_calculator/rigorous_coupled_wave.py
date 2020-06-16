import numpy as np
import tmm
from solcore.structure import Layer
from solcore.absorption_calculator import OptiStack
from joblib import Parallel, delayed
from copy import deepcopy
try:
    import S4
except ModuleNotFoundError:
    raise

default_options = dict(LatticeTruncation='Circular',
                        DiscretizedEpsilon=False,
                        DiscretizationResolution=8,
                        PolarizationDecomposition=False,
                        PolarizationBasis='Default',
                        LanczosSmoothing=False,
                        SubpixelSmoothing=False,
                        ConserveMemory=False,
                        WeismannFormulation=False)


def calculate_rat_rcwa(structure, size, orders, wavelength, incidence, substrate, theta=0, phi=0, pol='u',
                       parallel=False, n_jobs=-1, user_options={}):
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
    # write a separate function that makes the OptiStack structure into an S4 object, defined materials etc.

    ## Materials for the shapes need to be defined before you can do .SetRegion

    geom_list = [layer.geometry for layer in structure]
    geom_list.insert(0, {})  # incidence medium
    geom_list.append({})  # transmission medium

    shape_mats, geom_list_str = necessary_materials(geom_list)

    shapes_oc = np.zeros((len(wavelength), len(shape_mats)), dtype=complex)

    for i1, x in enumerate(shape_mats):
        shapes_oc[:, i1] = (x.n(wavelength/1e9) + 1j * x.k(wavelength/1e9)) ** 2

    stack_OS = OptiStack(structure, bo_back_reflection=False, substrate=substrate)
    widths = stack_OS.get_widths()
    layers_oc = np.zeros((len(wavelength/1e9), len(structure) + 2), dtype=complex)

    if incidence == None:
        layers_oc[:, 0] = 1
    else:
        layers_oc[:, 0] = (incidence.n(wavelength/1e9)) ** 2  # + 1j*incidence.k(wavelengths))**2


    if substrate == None:
        layers_oc[:, -1] = 1
    else:
        layers_oc[:, -1] = (substrate.n(wavelength/1e9) + 1j * substrate.k(wavelength/1e9)) ** 2

    for i1, x in enumerate(structure):
        layers_oc[:, i1 + 1] = (x.material.n(wavelength/1e9) + 1j * x.material.k(wavelength/1e9)) ** 2

    shapes_names = [str(x) for x in shape_mats]

    rcwa_options = default_options
    rcwa_options.update(user_options)

    if parallel:
        allres = Parallel(n_jobs=n_jobs)(delayed(RCWA_wl)
                                                         (wavelength[i1], geom_list_str,
                                                          layers_oc[i1], shapes_oc[i1],
                                                          shapes_names, pol,
                                                          theta, phi,
                                                          widths, size,
                                                          orders, rcwa_options)
                                                         for i1 in range(num_wl))

    else:
        allres = [
            RCWA_wl(wavelength[i1], geom_list_str, layers_oc[i1], shapes_oc[i1],
                         shapes_names, pol, theta, phi,
                         widths, size,
                         orders, rcwa_options)
            for i1 in range(num_wl)]

    R = np.stack([item[0] for item in allres])
    T = np.stack([item[1] for item in allres])
    A_mat = np.stack([item[2] for item in allres])

    output = {'R': R, 'A': np.sum(A_mat, 1), 'T': T,
              'A_per_layer': A_mat}

    return output


def rcwa_rat(S, n_layers):
    below = 'layer_' + str(n_layers)  # identify which layer is the transmission medium

    R = 1 - sum(S.GetPowerFlux('layer_2'))  # GetPowerFlux gives forward & backward Poynting vector, so sum to get power flux
    R_pfbo = -np.array(S.GetPowerFluxByOrder('layer_1'))[:,1] # real part of backwards power flow
    Nrm = np.real(np.sum(R_pfbo))
    R_pfbo = np.real((R/Nrm)*R_pfbo)

    T = sum(S.GetPowerFlux(below))

    T_pfbo = np.real(np.sum(np.array(S.GetPowerFluxByOrder(below)), 1))
    return {'R': np.real(R), 'T': np.real(T), 'R_pfbo': R_pfbo, 'T_pfbo': T_pfbo}


def initialise_S(size, orders, geom_list, mats_oc, shapes_oc, shape_mats, widths, options):
    # pass widths
    #print(widths)
    S = S4.New(size, orders)

    S.SetOptions(
        LatticeTruncation = options['LatticeTruncation'],
        DiscretizedEpsilon = options['DiscretizedEpsilon'],
        DiscretizationResolution = options['DiscretizationResolution'],
        PolarizationDecomposition = options['PolarizationDecomposition'],
        PolarizationBasis = options['PolarizationBasis'],
        LanczosSmoothing = options['LanczosSmoothing'],
        SubpixelSmoothing = options['SubpixelSmoothing'],
        ConserveMemory = options['ConserveMemory'],
        WeismannFormulation = options['WeismannFormulation']
    )


    for i1 in range(len(shapes_oc)):
        S.SetMaterial('shape_mat_' + str(i1 + 1), shapes_oc[i1])

    for i1 in range(len(widths)):
        S.SetMaterial('layer_' + str(i1 + 1), mats_oc[i1])

    for i1 in range(len(widths)):  # set base layers
        layer_name = 'layer_' + str(i1 + 1)
        #print(layer_name)
        if widths[i1] == float('Inf'):
            #print('zero width')
            S.AddLayer(layer_name, 0, layer_name)  # Solcore4 has incidence and transmission media widths set to Inf;
            # in S4 they have zero width
        else:
            S.AddLayer(layer_name, widths[i1], layer_name)

        geometry = geom_list[i1]

        if bool(geometry):
            for shape in geometry:

                mat_name = 'shape_mat_' + str(shape_mats.index(str(shape['mat'])) + 1)
                #print(str(shape['mat']), mat_name)
                if shape['type'] == 'circle':
                    S.SetRegionCircle(layer_name, mat_name, shape['center'], shape['radius'])
                elif shape['type'] == 'ellipse':
                    S.SetRegionEllipse(layer_name, mat_name, shape['center'], shape['angle'], shape['halfwidths'])
                elif shape['type'] == 'rectangle':
                    #print('rect')
                    S.SetRegionRectangle(layer_name, mat_name, shape['center'], shape['angle'], shape['halfwidths'])
                elif shape['type'] == 'polygon':
                    S.SetRegionPolygon(layer_name, mat_name, shape['center'], shape['angle'], shape['vertices'])

    return S


def necessary_materials(geom_list):
    shape_mats = []
    geom_list_str = deepcopy(geom_list)
    for i1, geom in enumerate(geom_list):
        if bool(geom):
            shape_mats.append([x['mat'] for x in geom])
            for i2, g in enumerate(geom):
                geom_list_str[i1][i2]['mat'] = str(g['mat'])
    return list(set([val for sublist in shape_mats for val in sublist])), geom_list_str


def update_epsilon(S, stack_OS, shape_mats_OS, wl):
    for i1 in range(len(stack_OS.get_widths())):
        S.SetMaterial('layer_' + str(i1 + 1), stack_OS.get_indices(wl)[i1] ** 2)
    for i1 in range(len(shape_mats_OS.widths)):  # initialise the materials needed for all the shapes in S4
        S.SetMaterial('shape_mat_' + str(i1 + 1), shape_mats_OS.get_indices(wl)[i1 + 1] ** 2)

    return S

def calculate_absorption_profile_rcwa(structure, size, orders, wavelength, rat_output_A,
                                      z_limit=None, steps_size=2, dist=None, theta=0, phi=0, pol='u', incidence=None,
                                      substrate=None,
                                      parallel=False, n_jobs=-1, user_options={}):
    """ It calculates the absorbed energy density within the material.
    Integrating this absorption profile in the whole stack gives the same result that the absorption obtained with
    calculate_rat as long as the spatial mesh is fine enough. If the structure is
    very thick and the mesh not thin enough, the calculation might diverge at short wavelengths.

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

    geom_list = [layer.geometry for layer in structure]
    geom_list.insert(0, {})  # incidence medium
    geom_list.append({})  # transmission medium
    # write a separate function that makes the OptiStack structure into an S4 object, defined materials etc.
    # write a separate function that makes the OptiStack structure into an S4 object, defined materials etc.

    ## Materials for the shapes need to be defined before you can do .SetRegion
    shape_mats, geom_list_str = necessary_materials(geom_list)

    shapes_oc = np.zeros((len(wavelength), len(shape_mats)), dtype=complex)

    for i1, x in enumerate(shape_mats):
        shapes_oc[:, i1] = (x.n(wavelength/1e9) + 1j * x.k(wavelength/1e9)) ** 2

    stack_OS = OptiStack(structure, bo_back_reflection=False, substrate=substrate)
    widths = stack_OS.get_widths()
    layers_oc = np.zeros((len(wavelength/1e9), len(structure) + 2), dtype=complex)

    if incidence == None:
        layers_oc[:, 0] = 1
    else:
        layers_oc[:, 0] = (incidence.n(wavelength/1e9)) ** 2  # + 1j*incidence.k(wavelengths))**2

    if substrate == None:
        layers_oc[:, -1] = 1
    else:
        layers_oc[:, -1] = (substrate.n(wavelength / 1e9) + 1j * substrate.k(wavelength / 1e9)) ** 2

    for i1, x in enumerate(structure):
        layers_oc[:, i1 + 1] = (x.material.n(wavelength/1e9) + 1j * x.material.k(wavelength/1e9)) ** 2

    shapes_names = [str(x) for x in shape_mats]
    rcwa_options = default_options
    rcwa_options.update(user_options)

    if dist is None:
        if z_limit is None:
            stack = OptiStack(structure)
            z_limit = np.sum(np.array(stack.widths))
        dist = np.arange(0, z_limit, steps_size)

    output = {'position': dist, 'absorption': np.zeros((num_wl, len(dist)))}

    if parallel:
        allres = Parallel(n_jobs=n_jobs)(delayed(RCWA_wl_prof)
                                                         (wavelength[i1], rat_output_A[i1],
                                                          dist,
                                                          geom_list_str,
                                                          layers_oc[i1], shapes_oc[i1],
                                                          shapes_names, pol,
                                                          theta, phi,
                                                          widths, size,
                                                          orders, rcwa_options)
                                                         for i1 in range(num_wl))

    else:
        allres = [
            RCWA_wl_prof(wavelength[i1], rat_output_A[i1],
                                                          dist,
                                                          geom_list_str,
                                                          layers_oc[i1], shapes_oc[i1],
                                                          shapes_names, pol,
                                                          theta, phi,
                                                          widths, size,
                                                          orders, rcwa_options)
            for i1 in range(num_wl)]

    profile = np.stack(allres)
    output['absorption'] = profile

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

    A = np.array([x if x > 0 else 0 for x in A])

    return A

def RCWA_wl(wl, geom_list, l_oc, s_oc, s_names, pol, theta, phi, widths, size, orders, rcwa_options):

    S = initialise_S(size, orders, geom_list, l_oc, s_oc, s_names, widths, rcwa_options)

    n_inc = np.real(np.sqrt(l_oc[0]))

    if len(pol) == 2:

        S.SetExcitationPlanewave((theta, phi), pol[0], pol[1], 0)
        S.SetFrequency(1 / wl)
        out = rcwa_rat(S, len(widths))
        R = out['R']
        T = out['T']
        A_layer = rcwa_absorption_per_layer(S, len(widths))

    else:
        if pol in 'sp':
            if pol == 's':
                s = 1
                p = 0
            elif pol == 'p':
                s = 0
                p = 1

            S.SetExcitationPlanewave((theta, phi), s, p, 0)
            S.SetFrequency(1 / wl)
            out = rcwa_rat(S, len(widths))
            R = out['R']
            T = out['T']
            A_layer = rcwa_absorption_per_layer(S, len(widths))

        else:

            S.SetFrequency(1 / wl)
            S.SetExcitationPlanewave((theta, phi), 0, 1, 0)  # p-polarization
            out_p = rcwa_rat(S, len(widths))
            A_layer_p = rcwa_absorption_per_layer(S, len(widths))


            S.SetExcitationPlanewave((theta, phi), 1, 0, 0)  # s-polarization
            out_s = rcwa_rat(S, len(widths))
            A_layer_s = rcwa_absorption_per_layer(S, len(widths))

            R = 0.5 * (out_p['R'] + out_s['R'])  # average
            T = 0.5 * (out_p['T'] + out_s['T'])
            A_layer = 0.5 * (A_layer_s + A_layer_p)

    T = n_inc * T / np.cos(theta * np.pi / 180)
    A_layer = n_inc * A_layer / np.cos(theta * np.pi / 180)
    return R, T, A_layer


def RCWA_wl_prof(wl, rat_output_A, dist, geom_list, layers_oc, shapes_oc, s_names, pol, theta, phi, widths, size, orders, rcwa_options):

    S = initialise_S(size, orders, geom_list, layers_oc, shapes_oc, s_names, widths, rcwa_options)
    profile_data = np.zeros(len(dist))


    A = rat_output_A

    if len(pol) == 2:

        S.SetExcitationPlanewave((theta, phi), pol[0], pol[1], 0)
        S.SetFrequency(1 / wl)
        for j, d in enumerate(dist):
            layer, d_in_layer = tmm.find_in_structure_with_inf(widths,
                                                               d)  # don't need to change this
            layer_name = 'layer_' + str(layer + 1)  # layer_1 is air above so need to add 1
            data = rcwa_position_resolved(S, layer_name, d_in_layer, A)
            profile_data[j] = data


    else:
        if pol in 'sp':
            if pol == 's':
                s = 1
                p = 0
            elif pol == 'p':
                s = 0
                p = 1

            S.SetExcitationPlanewave((theta, phi), s, p, 0)


            S.SetFrequency(1 / wl)

            for j, d in enumerate(dist):
                layer, d_in_layer = tmm.find_in_structure_with_inf(widths,
                                                                   d)  # don't need to change this
                layer_name = 'layer_' + str(layer + 1)  # layer_1 is air above so need to add 1
                data = rcwa_position_resolved(S, layer_name, d_in_layer, A)
                profile_data[j] = data

        else:


            S.SetFrequency(1 / wl)
            A = rat_output_A

            for j, d in enumerate(dist):
                layer, d_in_layer = tmm.find_in_structure_with_inf(widths,
                                                                   d)  # don't need to change this
                layer_name = 'layer_' + str(layer + 1)  # layer_1 is air above so need to add 1
                S.SetExcitationPlanewave((theta, phi), 0, 1, 0)  # p-polarization
                data_p = rcwa_position_resolved(S, layer_name, d_in_layer, A)
                S.SetExcitationPlanewave((theta, phi), 1, 0, 0)  # p-polarization
                data_s = rcwa_position_resolved(S, layer_name, d_in_layer, A)
                profile_data[j] = 0.5*(data_s + data_p)

    return profile_data

