import numpy as np
import tmm
from solcore.absorption_calculator import OptiStack, RCWASolverError
from joblib import Parallel, delayed

try:
    import S4
except ModuleNotFoundError:
    raise RCWASolverError


DEFAULT_OPTIONS = dict(LatticeTruncation='Circular',
                        DiscretizedEpsilon=False,
                        DiscretizationResolution=8,
                        PolarizationDecomposition=False,
                        PolarizationBasis='Default',
                        LanczosSmoothing=False,
                        SubpixelSmoothing=False,
                        ConserveMemory=False,
                        WeismannFormulation=False)


def calculate_rat_rcwa(structure, size, orders, wavelength, incidence, substrate, theta=0, phi=0, pol='u',
                       parallel=False, n_jobs=-1, user_options=None):
    """Calculates the reflected, absorbed and transmitted intensity of the structure for the wavelengths and angles
    defined using an RCWA method implemented using the S4 package.
    This function is analogous to calculate_rat from the transfer_matrix module.

    :param structure: A Solcore Structure object with layers and materials or a OptiStack object.
    :param size: a tuple of 2-D vectors in the format ((ux, uy), (vx, vy)) giving the x and y components of the lattice unit vectors in nm.
    :param orders: number of orders to retain in the RCWA calculations.
    :param wavelength: Wavelengths (in nm) in which calculate the data.
    :param incidence: a Solcore material describing the semi-infinite incidence medium
    :param substrate: a Solcore material describing the semi-infinite transmission medium
    :param theta: polar incidence angle (in degrees) of the incident light. Default: 0 (normal incidence)
    :param phi: azimuthal incidence angle (in degrees). Default: 0
    :param pol: Polarisation of the light: 's', 'p' or 'u'. Default: 'u' (unpolarised).
    :param parallel: whether or not to execute calculation in parallel (over wavelengths), True or False. Default is False
    :param n_jobs: the 'n_jobs' argument passed to Parallel from the joblib package. If set to -1, all available CPUs are used,
            if set to 1 no parallel computing is executed. The number of CPUs used is given by n_cpus + 1 + n_jobs. Default is -1.
    :param user_options: dictionary of options for S4. The list of possible entries and their values is:

            * LatticeTruncation: 'Circular' or 'Parallelogramic' (default 'Circular')
            * DiscretizedEpsilon: True or False (default False)
            * DiscretizationResolution: integer (default value 8)
            * PolarizationDecomposition: True or False (default False)
            * PolarizationBasis: 'Default' or 'Normal' or 'Jones' (default 'Default')
            * LanczosSmoothing: True or False (default False)
            * SubpixelSmoothing: True or False (default False)
            * ConserveMemory: True or False (default False)
            * WeismannFormulation: True or False (default False)

        Further information on the function of these options can be found in the S4 Python API documentation. If no
        options are provided, those from DEFAULT_OPTIONS are used.

    :return: A dictionary with the R, total A and T at the specified wavelengths and angle. A_pol lists total absorption
       at s and p polarizations if pol = 'u' was specified (this information is needed for the absorption profile calculation).
       Otherwise, A_pol is the same as A.
    """
    num_wl = len(wavelength)

    # write a separate function that makes the OptiStack structure into an S4 object, defined materials etc.
    # Materials for the shapes need to be defined before you can do .SetRegion

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

    rcwa_options = DEFAULT_OPTIONS
    if user_options is not None:
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
    if pol == 'u':
        A_per_layer = np.mean(A_mat, axis=1)
        A_pol = np.sum(A_mat, 2)
    else:
        A_per_layer = A_mat
        A_pol = np.sum(A_per_layer, 1)


    output = {'R': R, 'A': np.sum(A_per_layer, 1), 'T': T,
              'A_per_layer': A_per_layer, 'A_pol': A_pol}

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
    """Makes an S4 simulation object using S4.New() according to the user's specified structure and options.

    :param size: a tuple of 2-D vectors in the format ((ux, uy), (vx, vy)) giving the x and y components of the lattice unit vectors in nm.
    :param orders: number of Fourier orders to be retained in the RCWA/FMM calculation
    :param: geom_list: list containing shape information fpr each layer. Format is a list of lists; entries of inner lists are one dictionary
            per shape containing the shape information relevant for that shape.
    :param: mats_oc: complex dielectric constant for each of the base layers (including incidence and transmission media)
    :param: shapes_oc: complex dielectric constants for materials in the shapes inside layers
    :param: shape_mats: names of the shape materials
    :param: widths: widths of the layers, including the semi-infinite incidence and transmission medium. In nm.
    :param: options: S4 options (dictionary)

    :return: S4 simulation object
    """
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


    for i1, shape in enumerate(shapes_oc):
        S.SetMaterial('shape_mat_' + str(i1 + 1), shape)

    for i1 in range(len(widths)):
        S.SetMaterial('layer_' + str(i1 + 1), mats_oc[i1])

    for i1, width in enumerate(widths):  # set base layers
        layer_name = 'layer_' + str(i1 + 1)
        #print(layer_name)
        if width == float('Inf'):
            #print('zero width')
            S.AddLayer(layer_name, 0, layer_name)  # Solcore4 has incidence and transmission media widths set to Inf;
            # in S4 they have zero width
        else:
            S.AddLayer(layer_name, width, layer_name)

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
    geom_list_str = [None] * len(geom_list)
    for i1, geom in enumerate(geom_list):
        if bool(geom):
            shape_mats.append([x['mat'] for x in geom])
            geom_list_str[i1] = [{}] * len(geom)
            for i2, g in enumerate(geom):
                for item in g.keys():
                    if item != 'mat':
                        geom_list_str[i1][i2][item] = g[item]
                    else:
                        geom_list_str[i1][i2][item] = str(g[item])

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
                                      parallel=False, n_jobs=-1, user_options=None):
    """It calculates the absorbed energy density within the material. Integrating this absorption profile in the whole stack gives the same result as the absorption obtained with
    calculate_rat as long as the spatial mesh is fine enough. If the structure is
    very thick and the mesh not thin enough, the calculation might diverge at short wavelengths.
    This function is analogous to calculate_absorption_profile from the transfer_matrix module.

    :param structure: A Solcore structure with layers and materials.
    :param size: a tuple of 2-D vectors in the format ((ux, uy), (vx, vy)) giving the x and y components of the lattice unit vectors in nm.
    :param orders: number of orders to retain in the RCWA calculations.
    :param wavelength: Wavelengths (in nm) in which calculate the data.
    :param rat_output: A_pol' (polarization-dependent layer absorption) output from calculate_rat_rcwa
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
    rcwa_options = DEFAULT_OPTIONS

    if user_options is not None:
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

    """
    Calculates the reflection, absorption per layer and transmission using RCWA per wavelength. Called from calculate_rat_rcwa,
    can be used in parallel with joblib.

    :param wl: wavelength in nm
    :param geom_list: list containing shape information fpr each layer. Format is a list of lists; entries of inner lists are one dictionary
        per shape containing the shape information relevant for that shape.
    :param l_oc: list of optical constants of the main layers in the stack
    :param s_oc: list of optical constants for the shapes inside the layers
    :param s_names: list of names of the shapes
    :param pol: polarization ('s', 'p' or 'u')
    :param theta: polar incidence angle in degrees
    :param phi: azimuthal angle (clockwise from y-axis) in degrees)
    :param widths: layer widths in nm. Incidence and transmission medium have 0 thickness.
    :param size: a tuple of 2-D vectors in the format ((ux, uy), (vx, vy)) giving the x and y components of the lattice unit vectors in nm.
    :param orders: number of orders to retain in the RCWA calculations.
    :param rcwa_options: S4 options (dictionary)
    :return: R, T and A per layer at the specified wavelength
    """

    def vs_pol(s, p):
        S.SetExcitationPlanewave((theta, phi), s, p, 0)
        S.SetFrequency(1 / wl)
        out = rcwa_rat(S, len(widths))
        R = out['R']
        T = out['T']
        A_layer = rcwa_absorption_per_layer(S, len(widths))
        return R, T, A_layer

    S = initialise_S(size, orders, geom_list, l_oc, s_oc, s_names, widths, rcwa_options)

    n_inc = np.real(np.sqrt(l_oc[0]))

    if len(pol) == 2:

        R, T, A_layer = vs_pol(pol[0], pol[1])

    else:
        if pol in 'sp':
            R, T, A_layer = vs_pol(int(pol == "s"), int(pol == "p"))

        else:
            R_s, T_s, A_layer_s = vs_pol(1, 0)
            R_p, T_p, A_layer_p = vs_pol(0, 1)
            R = (R_s + R_p) / 2
            T = (T_s + T_p) / 2
            A_layer = np.stack([A_layer_s, A_layer_p])

    T = n_inc * T / np.cos(theta * np.pi / 180)
    A_layer = n_inc * A_layer / np.cos(theta * np.pi / 180)
    return R, T, A_layer


def RCWA_wl_prof(wl, A, dist, geom_list, layers_oc, shapes_oc, s_names, pol, theta, phi, widths, size, orders, rcwa_options):
    """
    Calculates the depth-dependent absorption profile using RCWA per wavelength. Called from calculate_absorption_profile_rcwa,
    can be used in parallel with joblib.

    :param wl: wavelength in nm
    :param A: 'A_pol' (polarization-dependent layer absorption) output from calculate_rat_rcwa
    :param geom_list: list containing shape information fpr each layer. Format is a list of lists; entries of inner lists are one dictionary
        per shape containing the shape information relevant for that shape.
    :param layers_oc: list of optical constants of the main layers in the stack
    :param shapes_oc: list of optical constants for the shapes inside the layers
    :param s_names: list of names of the shapes
    :param pol: polarization ('s', 'p' or 'u')
    :param theta: polar incidence angle in degrees
    :param phi: azimuthal angle (clockwise from y-axis) in degrees)
    :param widths: layer widths in nm. Incidence and transmission medium have 0 thickness.
    :param size: a tuple of 2-D vectors in the format ((ux, uy), (vx, vy)) giving the x and y components of the lattice unit vectors in nm.
    :param orders: number of orders to retain in the RCWA calculations.
    :param rcwa_options: S4 options (dictionary)
    :return: profile data (absorbed energy density at each depth)
    """
    def vs_pol_prof(s, p, A_total):
        profile = np.zeros(len(dist))
        S.SetExcitationPlanewave((theta, phi), s, p, 0)
        S.SetFrequency(1 / wl)
        for j, d in enumerate(dist):
            layer, d_in_layer = tmm.find_in_structure_with_inf(widths, d)
            layer_name = 'layer_' + str(layer + 1)  # layer_1 is air above so need to add 1
            data = rcwa_position_resolved(S, layer_name, d_in_layer, A_total)
            profile[j] = data

        return profile

    S = initialise_S(size, orders, geom_list, layers_oc, shapes_oc, s_names, widths, rcwa_options)
    n_inc = np.real(np.sqrt(layers_oc[0]))
    if len(pol) == 2:

        profile_data = vs_pol_prof(pol[0], pol[1], A)

    else:
        if pol in 'sp':
            profile_data = vs_pol_prof(int(pol == "s"), int(pol == "p"), A)

        else:
            profile_s = vs_pol_prof(1, 0, A[0])
            profile_p = vs_pol_prof(0, 1, A[1])
            profile_data = 0.5*(profile_s + profile_p)

    profile_data = n_inc * profile_data / np.cos(theta * np.pi / 180)
    return profile_data
