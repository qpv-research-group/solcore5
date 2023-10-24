from sesame import Builder, IVcurve, Analyzer
import numpy as np
from scipy.optimize import root
from solcore.constants import q
from solcore.light_source import LightSource
from scipy.interpolate import interp1d
from solcore.state import State

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

# sesame_options = State()
# sesame_options.sesame_minimum_spacing = 1e-9 # in m: 1 nm
# sesame_options.sesame_maximum_spacing = 1e-7 # in m: 100 nm

def process_structure(junction, options):

    # get material parameters from the junction and convert them to format required by Sesame (dictionary)
    # convert from Solcore units (i.e. base SI units like m) to Sesame's units (cm, eV).
    # Note that internally in sesame, and its outputs, many quantities are scaled/dimensionless.

    # acceptable formats would be:
    # single material with a doping profile
    # two materials (same or different underlying base material, i.e. homojunction or heterojunction) with different,
    # constant doping profiles n and p-type
    # two materials with a doping profile
    # Should be able to set a doping profile for the whole junction OR per layer in the junction

    materials = []
    layer_widths = []
    doping_profile_functions = []
    offset = 0

    for layer in junction:

        try:
            bulk_recombination_energy = layer.material.bulk_recombination_energy

        except:
            bulk_recombination_energy = 0

        layer_widths.append(layer.width*1e2) # m to cm

        constant_doping = layer.material.Nd / 1e6 if layer.material.Nd > layer.material.Na else -layer.material.Na /1e6

        if hasattr(layer, "doping_profile"):
            doping_profile_functions.append(lambda x: layer.doping_profile(x) / 1e6)
            # function which takes an argument of depth in junction in m and returns doping in cm-3

        else:
            doping_profile_functions.append(interp1d([offset/100, offset + layer.width],
                                                     [constant_doping, constant_doping]))
            # to be consistent with functions passed by user, this should be a function which takes an argument in m
            # and returns doping in cm-3


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
            'Cn': layer.material.electron_auger_recombination*1e12, # Auger recombination constant for electrons (cm6/s)
            'Cp': layer.material.hole_auger_recombination*1e12, # Auger recombination constant for holes (cm6/s)
        }

        offset += layer.width*1e2
        # Note: all of these can be functions of position as well as constants!

        materials.append(new_mat)

    # see if we can determine the junction depth:

    edges = np.insert(np.cumsum(layer_widths), 0, 0)
    edges[-1] = edges[-1] + 1e-10 # otherwise final point will not be assigned any values

    print('edges', edges)

    if hasattr(junction, "doping_profile"):
        junction_depth = 100*root(junction.doping_profile, np.sum(layer_widths) / (2*100)).x # in cm

    else:
        doping_sign = [np.sign(doping_profile_functions[i1](edges[i1]+1e-9)) for i1 in range(len(edges) - 1)]
        junction_depth = edges[np.where(np.diff(doping_sign))[0] + 1]

    if not hasattr(junction, "mesh"):
        print("set mesh")
        make_mesh(junction, layer_widths, options, junction_depth)

    else:
        junction.mesh_cm = junction.mesh*1e2 # m to cm

    print('junction depth', junction_depth)
    print('mesh points', len(junction.mesh))

    if hasattr(junction, "doping_profile"):
        doping_profile_x = junction.doping_profile(junction.mesh) / 1e6 # in cm-3

    else:

        doping_profile_x = np.zeros(len(junction.mesh_cm))
        for i1, doping_profile_fn in enumerate(doping_profile_functions):
            doping_profile_x[np.all((junction.mesh_cm >= edges[i1], junction.mesh_cm < edges[i1 + 1]), axis=0)] = (
                doping_profile_fn(junction.mesh[np.all((junction.mesh_cm >= edges[i1], junction.mesh_cm < edges[i1 + 1]), axis=0)]))

        # doping_profile = np.vectorize(doping_profile)
    # junction.mesh_cm = junction.mesh*1e2 # m to cm

    junction.sesame_sys = Builder(junction.mesh_cm) # Sesame system

    junction.sesame_sys.rho = doping_profile_x / junction.sesame_sys.scaling.density

    # import matplotlib.pyplot as plt
    # plt.figure()
    # plt.plot(junction.mesh_cm, doping_profile_x)
    # plt.show()
    # print(doping_profile_x)

    for i1 in range(len(junction)):
        junction.sesame_sys.add_material(materials[i1],
                                         lambda x: np.all((x >= edges[i1], x < edges[i1 + 1]), axis=0))

    # if constant_doping:
    #     doping_profile = np.zeros_like(junction.mesh_cm)
    #     for i1, dop in enumerate(doping_list):
    #         doping_profile[np.all((junction.mesh_cm >= edges[i1], junction.mesh_cm <= edges[i1 + 1]), axis=0)] = dop
    #
    # else:
    #     doping_profile = junction.doping_profile

    junction.sesame_sys.contact_type('Ohmic', 'Ohmic')
    # should also be able to choose other options, since Sesame can handle them

    # lower than 1e-6: diverges
    junction.sesame_sys.contact_S(*get_srv(junction))


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
        edges = np.insert(np.cumsum(layer_width), 0, 0)

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

        # identify other points of interest and make mesh dense around them
        #
        # for i1, layer_edge in enumerate(edges[1:]):
        #     start_dense = layer_edge - 0.1*layer_width[i1]
        #
        #     if i1 < len(layer_width) - 1:
        #         end_dense = layer_edge + 0.1*layer_width[i1+1]
        #
        #     else:
        #         end_dense = total_width
        #
        #     mesh_above = mesh[mesh <= start_dense]
        #     mesh_dense = np.arange(start_dense, end_dense, minimum_spacing)
        #     mesh_below = mesh[mesh >= end_dense]
        #     mesh = np.unique(np.concatenate((mesh_above, mesh_dense, mesh_below)))
        #
        #     # plt.figure()
        #     # plt.plot(mesh)
        #     # for layer_edge in np.cumsum(layer_width):
        #     #     plt.axhline(layer_edge, color='k')
        #     # plt.show()
        #
        # if len(junction_depth) > 1:
        #     for junc_depth in junction_depth[1:]:
        #
        #         start_dense = junc_depth - 0.1 * total_width
        #         start_dense = 0 if start_dense < 0 else start_dense
        #
        #         end_dense = junc_depth + 0.1*total_width
        #         end_dense = total_width if end_dense > total_width else end_dense
        #
        #         mesh_above = mesh[mesh <= start_dense]
        #         mesh_dense = np.arange(start_dense, end_dense, minimum_spacing)
        #         mesh_below = mesh[mesh >= end_dense]
        #         mesh = np.unique(np.concatenate((mesh_above, mesh_dense, mesh_below)))
        #

    junction.mesh_cm = mesh
    junction.mesh = junction.mesh_cm*1e-2 # cm to m

def get_srv(junction):
    # Define the surface recombination velocities for electrons and holes [cm/s]
    if hasattr(junction, "sn_front"):
        Sn_left = junction.sn_front / 100

    elif hasattr(junction, "sn"):
        Sn_left = junction.sn / 100

    else:
        Sn_left = 1e4 # is this the same as the Fortran PDD?

    if hasattr(junction, "sp_front"):
        Sp_left = junction.sp_front / 100

    elif hasattr(junction, "sp"):
        Sp_left = junction.sp / 100

    else:
        Sp_left = 1e4

    if hasattr(junction, "sn_rear"):
        Sn_right = junction.sn_front / 100

    elif hasattr(junction, "sn"):
        Sn_right = junction.sn / 100

    else:
        Sn_right = 1e4 # is this the same as the Fortran PDD?

    if hasattr(junction, "sp_rear"):
        Sp_right = junction.sp_front / 100

    elif hasattr(junction, "sp"):
        Sp_right = junction.sp / 100

    else:
        Sp_right = 1e4

    return Sn_left, Sp_left, Sn_right, Sp_right

@register_equilibrium_solver("sesame_PDD", reason_to_exclude=reason_to_exclude)
def equilibrium():
    pass

@register_iv_solver("sesame_PDD", reason_to_exclude=reason_to_exclude)
def iv_sesame(junction, options):

    # TODO: inclusion of shunt resistance
    if ~hasattr(junction, "sesame_sys"):
        process_structure(junction, options)

    # gen_wl = profile_func


    if options.light_iv:

        gen_wl = junction.absorbed(junction.mesh) / 100 # m-1 -> cm-1
        wls = options.wavelength

        gg = options.light_source.spectrum(wls, output_units="photon_flux_per_m")[1][:, None] * gen_wl.T

        g_vs_z = np.trapz(gg, wls, axis=0) / 1e4 # m^2 -> cm^2
        g_vs_z[np.isnan(g_vs_z)] = 0

        # can also pass a function to generation - more flexible?
        junction.sesame_sys.generation(g_vs_z)

    # equilibrium solution:
    # j0, result0 = sesame.IVcurve(junction.sesame_sys, [0])

    voltages = options.voltages

    j, result = sesame.IVcurve(junction.sesame_sys, voltages)  # , verbose=False)
    j = j * junction.sesame_sys.scaling.current * 1e4 # cm-2 -> m-2

    # Jsc = j[0] # units?
    junction.current = j

    junction.iv = interp1d(
        voltages,
        j,
        kind="linear",
        bounds_error=False,
        assume_sorted=True,
        fill_value=(junction.current[0], junction.current[-1]),
    )

    junction.voltage = voltages
    # jsc = q*np.trapz(eqe*light_source.spectrum()[1], wls) # A/m2

    # this is all done externally, same for all junction types
    # zero_crossing = np.where(np.diff(np.sign(j)))[0][0]
    # j_above = j[zero_crossing]
    # j_below = j[zero_crossing + 1]
    #
    # Voc = voltages[zero_crossing] + (voltages[zero_crossing + 1] - voltages[zero_crossing]) * j_above / (
    #             j_above - j_below)
    #
    # voltages_for_mpp = voltages * (np.abs(voltages) <= Voc)
    # Vmpp = voltages[np.nanargmax(np.abs(j * voltages_for_mpp))]
    # Jmpp = j[np.nanargmax(np.abs(j * voltages_for_mpp))] * 1e4
    #
    # FF = Vmpp * Jmpp / (Jsc * Voc)


def qe_sesame(junction, options):

    if ~hasattr(junction, "sesame_sys"):
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
        print('mesh shape', junction.mesh_cm.shape)

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


    plt.figure()
    plt.plot(wls, A)
    plt.plot(wls, junction.layer_absorption, '--')
    plt.plot(wls, eqe)
    plt.show()

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



