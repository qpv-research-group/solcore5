from sesame import Builder, IVcurve, Analyzer
import numpy as np
from scipy.optimize import root
from solcore.constants import q
from scipy.interpolate import interp1d
from solcore.state import State
import warnings

from solcore.registries import (
    register_iv_solver,
    register_equilibrium_solver
)

try:
    import sesame

    reason_to_exclude = None
except ImportError:
    reason_to_exclude = (
        "Sesame was not installed."
    )

def process_structure(junction, options):

    """
    Process the layers in the junction and convert it to a format which can be used by Sesame. This includes unit
    conversions from Solcore base SI units to the units used by Sesame (eV, cm). This function also sets the mesh,
    either from a user-supplied junction.mesh (in m) attribute, or by guessing an appropriate mesh depending on the
    doping profile and layer thicknesses. It will also extract the doping profile, which can be supplied in
    the following formats:

        - a constant value for each layer, set by the Nd or Na attribute of the material
        - a doping profile for the whole junction, passed as a function which accepts a depth in the junction in m
          and returns doping in m-3
        - a doping profile for each layer, passed as a function which accepts a depth in the layer in m and returns
          doping in m-3
        - Combinations of the above

    Note that since we can specify a doping profile for the whole junction, the Sesame solver can support a junction
    which is made up of a single layer, unlike the depletion approximation and Fortran-based PDD solvers.

    This function will add or update the following attributes in the Junction object:

        - mesh: the mesh used by the Sesame solver, in m
        - mesh_cm: the mesh used by the Sesame solver, in cm
        - sesame_sys: a Sesame Builder object defining the junction in Sesame format, used to solve the IV and QE

    :param junction: a Junction object
    :param options: a State object containing options for the solver

    """

    materials = [] # list of dictionaries for Sesame
    layer_widths = [] # layer widths in cm
    doping_profile_functions = [] # list of functions which take an argument of depth in junction in m and
    # return doping in cm-3
    offset = 0 # offset relative to front of junction of the current layer in cm

    for layer in junction:

        try:
            bulk_recombination_energy = layer.material.bulk_recombination_energy

        except:
            bulk_recombination_energy = 0 # should this be the degau

        layer_widths.append(layer.width*1e2) # m to cm

        # this will always be set, even if the layer has a doping profile
        constant_doping = layer.material.Nd / 1e6 if layer.material.Nd > layer.material.Na else -layer.material.Na /1e6

        if hasattr(layer, "doping_profile"):
            doping_profile_functions.append(lambda x: layer.doping_profile(x) / 1e6)
            # function which takes an argument of depth in junction in m and returns doping in cm-3

        else:
            doping_profile_functions.append(interp1d([offset/100, offset + layer.width],
                                                     [constant_doping, constant_doping]))
            # to be consistent with functions passed by user, this should be a function which takes an argument in m
            # and returns doping in cm-3

        try:
            auger_electron = layer.material.electron_auger_recombination * 1e12 # m6/s to cm6/s

        except:
            auger_electron = 0

        try:
            auger_hole = layer.material.hole_auger_recombination * 1e12 # m6/s to cm6/s

        except:
            auger_hole = 0

        new_mat = {
            'Nc': layer.material.Nc * 1e-6, # effective density of states at CB edge (cm-3)
            'Nv': layer.material.Nv * 1e-6, # effective density of states at VB edge (cm-3)
            'Eg': layer.material.band_gap / q, # material bandgap (eV)
            'affinity': layer.material.electron_affinity / q, # electron affinity (eV)
            'epsilon': layer.material.relative_permittivity, # relative permittivity (dimensionless)
            'mu_e': layer.material.electron_mobility*1e4, # electron mobility (cm2/(V s))
            'mu_h': layer.material.hole_mobility*1e4, # hole mobility (cm2/(V s))
            'tau_e': layer.material.electron_minority_lifetime, # electron bulk lifetime (s)
            'tau_h': layer.material.hole_minority_lifetime, # hole bulk lifetime (s)
            'Et': bulk_recombination_energy, # energy level of bulk recombination centres (eV)
            'B': layer.material.radiative_recombination*1e6, # radiative recombination constant (cm3/s)
            'Cn': auger_electron, # Auger recombination constant for electrons (cm6/s)
            'Cp': auger_hole, # Auger recombination constant for holes (cm6/s)
        }

        offset += layer.width*1e2
        # Note: in Sesame, all of these can be functions of position as well as constants - but this is untested!

        materials.append(new_mat)

    # see if we can determine the junction depth:

    edges = np.insert(np.cumsum(layer_widths), 0, 0) # locations of interfaces between layers
    edges[-1] = edges[-1] + 1e-10 # otherwise final point will not be assigned any values

    if hasattr(junction, "doping_profile"):
        junction_depth = 100*root(junction.doping_profile, np.sum(layer_widths) / (2*100)).x # in cm
        # this is where doping crosses zero - not necessarily the junction depth, but not sure how best to determine
        # this. Maybe want to find slope of doping profile and identify locations where it is steepest?

    else:
        # identify where the doping profile changes sign for constant doping layers or individual doping profiles
        # per layer
        doping_sign = [np.sign(doping_profile_functions[i1](edges[i1]+1e-9)) for i1 in range(len(edges) - 1)]
        junction_depth = edges[np.where(np.diff(doping_sign))[0] + 1]

    if not hasattr(junction, "mesh"):
        # calculate appropriate mesh (hopefully) if not supplied as input in junction.mesh
        make_mesh(junction, layer_widths, options, junction_depth)
        user_mesh = False

    else:
        # user-provided mesh
        junction.mesh_cm = junction.mesh*1e2 # m to cm
        user_mesh = True

    if hasattr(junction, "doping_profile"):
        doping_profile_x = junction.doping_profile(junction.mesh) / 1e6 # in cm-3

    else:
        doping_profile_x = np.zeros(len(junction.mesh_cm))
        for i1, doping_profile_fn in enumerate(doping_profile_functions):
            doping_profile_x[np.all((junction.mesh_cm >= edges[i1], junction.mesh_cm < edges[i1 + 1]), axis=0)] = (
                doping_profile_fn(junction.mesh[np.all((junction.mesh_cm >= edges[i1], junction.mesh_cm < edges[i1 + 1]), axis=0)]))

    # Make Sesame Builder object for simulations
    junction.sesame_sys = Builder(junction.mesh_cm, T=options.T) # Sesame system

    junction.sesame_sys.rho = doping_profile_x / junction.sesame_sys.scaling.density

    if not user_mesh:
        update_mesh(junction, layer_widths, options)

        if hasattr(junction, "doping_profile"):
            doping_profile_x = junction.doping_profile(junction.mesh) / 1e6  # in cm-3

        else:
            doping_profile_x = np.zeros(len(junction.mesh_cm))
            for i1, doping_profile_fn in enumerate(doping_profile_functions):
                doping_profile_x[np.all((junction.mesh_cm >= edges[i1], junction.mesh_cm < edges[i1 + 1]), axis=0)] = (
                    doping_profile_fn(junction.mesh[
                                          np.all((junction.mesh_cm >= edges[i1], junction.mesh_cm < edges[i1 + 1]),
                                                 axis=0)]))

        junction.sesame_sys = Builder(junction.mesh_cm, T=options.T)  # Sesame system

        junction.sesame_sys.rho = doping_profile_x / junction.sesame_sys.scaling.density

    for i1 in range(len(junction)):
        junction.sesame_sys.add_material(materials[i1],
                                         lambda x: np.all((x >= edges[i1], x < edges[i1 + 1]), axis=0))


    junction.sesame_sys.contact_type('Ohmic', 'Ohmic')
    # should also be able to choose other options, since Sesame can handle them

    # get surface recombination velocities
    junction.sesame_sys.contact_S(*get_srv(junction))
    print('mesh points:', len(junction.mesh_cm))


def make_mesh(junction, layer_width, options, junction_depth):
    # want the mesh to be dense close to the front surface, where the doping is changing rapidly/at the junction,
    # and where material changes happen.
    # (changing from n-type to p-type).
    # dense linear mesh near surface and junction
    # less dense linear mesh in rest of bulk

    total_width = np.sum(layer_width)

    # always do 1 nm spacing near front surface?

    if "minimum_spacing" in options.keys():
        minimum_spacing = options.minimum_spacing*1e2

    else:
        minimum_spacing = total_width / 1e4 # maximum of 10000 points

    if "maximum_spacing" in options.keys():
        maximum_spacing = options.maximum_spacing*1e2

    else:
        maximum_spacing = total_width / 500 # minimum of 500 points

    print(minimum_spacing, maximum_spacing)

    if maximum_spacing == minimum_spacing:
        mesh = np.arange(0, total_width, maximum_spacing)

    # identify location of junction. doping_profile is a function

    else:
        # edges = np.insert(np.cumsum(layer_width), 0, 0)

        front_spacing = np.min([minimum_spacing, 0.5e-7])

        dense_front_width = 50e-7

        end_front = np.min([dense_front_width - 1e-10, total_width])

        front = np.arange(0, end_front, front_spacing)

        # 1 nm in cm = 1e-7
        if end_front != total_width:

            if junction_depth[0] < total_width / 2:
                print('front junction')
                end_dense = junction_depth[0] + 0.1 * total_width if (
                        junction_depth[0] + 0.1 * total_width < total_width) else total_width

                include_end = True if end_dense == total_width else False

                n_points = int((end_dense - end_front) / minimum_spacing)
                mesh = np.linspace(end_front, end_dense, n_points, endpoint=include_end)

                mesh = np.concatenate((front, mesh))

                if not include_end:
                    n_points = int((total_width - end_dense) / maximum_spacing)
                    mesh_rest = np.linspace(end_dense, total_width, n_points, endpoint=True)

                    mesh = np.concatenate((mesh, mesh_rest))

            else:

                n_points = int((0.2 * total_width- end_front) / minimum_spacing)

                mesh_dense_front = np.linspace(end_front, 0.2 * total_width, n_points, endpoint=False)

                n_points = int((junction_depth[0] - 0.3 * total_width) / maximum_spacing)

                mesh_bulk = np.linspace(0.2 * total_width, junction_depth[0] - 0.1 * total_width, n_points, endpoint=False)

                end_dense = junction_depth[0] + 0.1 * total_width if (

                        junction_depth[0] + 0.1 * total_width < total_width) else total_width

                include_end = True if end_dense == total_width else False

                n_points = int((end_dense - (junction_depth[0] - 0.1 * total_width)) / minimum_spacing)

                mesh_dense = np.linspace(junction_depth[0] - 0.1 * total_width, end_dense, n_points, endpoint=include_end)

                if not include_end:
                    n_points = int((total_width - end_dense) / maximum_spacing)

                    mesh_rest = np.linspace(end_dense, total_width, n_points, endpoint=True)

                    mesh = np.unique(np.concatenate((front, mesh_dense_front, mesh_bulk, mesh_dense, mesh_rest)))

                else:
                    mesh = np.unique(np.concatenate((front, mesh_dense_front, mesh_bulk, mesh_dense)))

    junction.mesh_cm = mesh
    junction.mesh = junction.mesh_cm*1e-2 # cm to m

def update_mesh(junction, layer_width, options):

    total_width = np.sum(layer_width)

    # always do 1 nm spacing near front surface?

    if "minimum_spacing" in options.keys():
        minimum_spacing = options.minimum_spacing*1e2

    else:
        minimum_spacing = total_width / 1e4 # maximum of 10000 points

    current_mesh = junction.mesh_cm

    diff_rho = np.abs(np.gradient(junction.sesame_sys.rho, current_mesh))

    doping_change = np.where(diff_rho / np.max(diff_rho) > 1e-4)[0]

    # identify continuous stretches of high doping change, and where there are gaps:

    change_point = np.where(np.diff(doping_change) != 1)[0]

    last_points = doping_change[change_point]
    first_points = doping_change[change_point + 1]

    first_points = np.insert(first_points, 0, doping_change[0])
    last_points = np.append(last_points, doping_change[-1])

    first_x = current_mesh[first_points]
    last_x = current_mesh[last_points]

    for i1 in range(len(first_points)):
        # want to make mesh denser where there are large changes in doping
        # print('denser mesh between', first_x[i1], last_x[i1])

        mesh_above = current_mesh[current_mesh < first_x[i1]]
        mesh_below = current_mesh[current_mesh >= last_x[i1]]
        mesh_dense = np.arange(first_x[i1], last_x[i1]-1e-12, minimum_spacing) # last point not included

        current_mesh = np.unique(np.concatenate((mesh_above, mesh_dense, mesh_below)))

    junction.mesh_cm = current_mesh
    junction.mesh = junction.mesh_cm*1e-2 # cm to m

def get_srv(junction):
    # Define the surface recombination velocities for electrons and holes [cm/s]
    if hasattr(junction, "sn_front"):
        Sn_left = junction.sn_front * 100

    elif hasattr(junction, "sn"):
        Sn_left = junction.sn * 100

    else:
        Sn_left = 1e4 # is this the same as the Fortran PDD?

    if hasattr(junction, "sp_front"):
        Sp_left = junction.sp_front * 100

    elif hasattr(junction, "sp"):
        Sp_left = junction.sp * 100

    else:
        Sp_left = 1e4

    if hasattr(junction, "sn_rear"):
        Sn_right = junction.sn_rear * 100

    elif hasattr(junction, "sn"):
        Sn_right = junction.sn * 100

    else:
        Sn_right = 1e4 # is this the same as the Fortran PDD?

    if hasattr(junction, "sp_rear"):
        Sp_right = junction.sp_rear * 100

    elif hasattr(junction, "sp"):
        Sp_right = junction.sp * 100

    else:
        Sp_right = 1e4

    return Sn_left, Sp_left, Sn_right, Sp_right

@register_equilibrium_solver("sesame_PDD", reason_to_exclude=reason_to_exclude)
def equilibrium():
    pass

@register_iv_solver("sesame_PDD", reason_to_exclude=reason_to_exclude)
def iv_sesame(junction, options):

    if not hasattr(junction, "sesame_sys"):
        process_structure(junction, options)

    # gen_wl = profile_func

    if options.light_iv and np.all(junction.sesame_sys.g == 0):

        gen_wl = junction.absorbed(junction.mesh) / 100 # m-1 -> cm-1
        wls = options.wavelength

        gg = options.light_source.spectrum(wls, output_units="photon_flux_per_m")[1][:, None] * gen_wl.T

        g_vs_z = np.trapz(gg, wls, axis=0) / 1e4 # m^2 -> cm^2
        g_vs_z[np.isnan(g_vs_z)] = 0

        # can also pass a function to generation - more flexible?
        junction.sesame_sys.generation(g_vs_z)

    voltages = options.internal_voltages

    R_shunt = min(junction.R_shunt, 1e14) if hasattr(junction, "R_shunt") else 1e14

    # voltages need to go from 0 (or close) to highest applied +ve or -ve voltage. Split if necessary.

    if junction.sesame_sys.rho[junction.sesame_sys.nx-1] < 0:
        # this is necessary because Sesame will internally flip the sign for an n-p junction
        voltages_for_solve = -voltages

        if np.all(options.voltages >= 0):
            warnings.warn('All voltages are positive, but junction has been identified as n-p, so the '
                          'open-circuit voltage (Voc) of the junction will be negative.', UserWarning)

    else:
        voltages_for_solve = voltages


    # split +ve and -ve voltages if necessary:
    if np.any(voltages_for_solve < 0):
        if np.any(voltages_for_solve > 0):

            negative_voltages = voltages_for_solve[voltages_for_solve <= 0]
            positive_voltages = voltages_for_solve[voltages_for_solve >= 0]

            negative_voltages_order = np.argsort(negative_voltages)[::-1]
            positive_voltages_order = np.argsort(positive_voltages)

            negative_voltages = negative_voltages[negative_voltages_order]
            positive_voltages = positive_voltages[positive_voltages_order]

            j_positive, result_positive = sesame.IVcurve(junction.sesame_sys, positive_voltages)
            j_negative, result_negative = sesame.IVcurve(junction.sesame_sys, negative_voltages)

            j_negative = j_negative[::-1]

            result_positive = {key: result_positive[key][positive_voltages_order, :] for key in result_positive.keys()}
            result_negative = {key: result_negative[key][negative_voltages_order, :] for key in result_negative.keys()}

            negative_voltages = negative_voltages[::-1]

            if np.any(voltages_for_solve == 0):

                # V = 0 would have been included in both the +ve and -ve voltages, so
                # exclude it from the negative voltage results when concatenating

                j = np.concatenate((j_negative[:-1], j_positive))
                result = {key: np.concatenate((result_negative[key][:-1], result_positive[key])) for key in
                          result_positive.keys()}
                final_voltages = np.concatenate((negative_voltages[:-1], positive_voltages))


            else:

                j = np.concatenate((j_negative, j_positive))
                result = {key: np.concatenate((result_negative[key], result_positive[key])) for key in result_positive.keys()}
                final_voltages = np.concatenate((negative_voltages, positive_voltages))

            # this results in j and result in order of increasing values for voltages_for_solve.


    else:
        voltage_order = np.argsort(voltages_for_solve)

        final_voltages = voltages_for_solve[voltage_order]
        j, result = sesame.IVcurve(junction.sesame_sys, final_voltages)  # , verbose=False)

    # final_voltages are the voltages corresponding to the entries in j and result, USING
    # Sesame's sign convention. So if the voltage sign was flipped above, need to flip it back
    # for Solcore

    if junction.sesame_sys.rho[junction.sesame_sys.nx-1] < 0:
        print('flipping back')
        result_voltage = -final_voltages[::-1]
        j = j[::-1]
        result = {key: result[key][::-1, :] for key in result.keys()}

    else:
        result_voltage = final_voltages

    j = j * junction.sesame_sys.scaling.current * 1e4 # cm-2 -> m-2

    junction.current = j + result_voltage / R_shunt

    # current_increasing = j[voltage_increasing]

    junction.iv = interp1d(
        result_voltage,
        junction.current,
        kind="linear",
        bounds_error=False,
        assume_sorted=True,
        fill_value=(j[0], j[-1]),
    )

    junction.voltage = result_voltage

def qe_sesame(junction, options):

    if not hasattr(junction, "sesame_sys"):
        process_structure(junction, options)

    wls = options.wavelength

    eqe = np.zeros_like(wls)

    voltages = [0]

    # profile_func = interp1d(bulk_positions_cm, 1e7 * Si_profile, kind='linear', bounds_error=False, fill_value=0)
    profile_func = junction.absorbed # this returns an array of shape (mesh_points, wavelengths)

    A = np.trapz(junction.absorbed(junction.mesh), junction.mesh, axis=0)

    def make_gfcn_fun(wl_index, flux):
        def gcfn_fun(x, y):
            # print('internal', profile_func(x / 100).shape)
            return flux * profile_func(x/100)[wl_index]/100 # convert to cm-1 from m-1
            # NOTE: this ONLY works if x and y are scalars!
        return gcfn_fun

    import matplotlib.pyplot as plt

    for i1, wl in enumerate(wls):

        flux = 1e20

        junction.sesame_sys.generation(make_gfcn_fun(i1, flux))
        # print(make_gfcn_fun(i1, flux)(junction.mesh_cm, 0).shape)
        # plt.figure()
        # plt.semilogy(junction.mesh_cm, [make_gfcn_fun(i1, flux)(x, 0) for x in junction.mesh_cm])
        # plt.xlim(0, 350e-7)
        # plt.show()

        if i1 == 0:
            guess = None

        else:
            guess = {key: result[key][0, :] for key in result.keys()}

        j, result = sesame.IVcurve(junction.sesame_sys, voltages, guess=guess)

        eqe[i1] = np.abs(j) / (q * flux)

    eqe = eqe * junction.sesame_sys.scaling.current
    iqe = eqe / A


    # plt.figure()
    # plt.plot(wls, A)
    # plt.plot(wls, junction.layer_absorption, '--')
    # plt.plot(wls, eqe)
    # plt.show()

    # convert dimensionless current to dimension-ful current
    # iqe = iqe * junction.sesame_sys.scaling.current

    junction.iqe = interp1d(wls, iqe)

    junction.eqe = interp1d(
        wls,
        eqe,
        kind="linear",
        bounds_error=False,
        assume_sorted=True,
        fill_value=(eqe[0], eqe[-1]),
    )

    junction.qe = State(
        {
            "WL": wl,
            "IQE": junction.iqe(wl),
            "EQE": junction.eqe(wl),
        }
    )



