from sesame import Builder, IVcurve, Analyzer
import numpy as np
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

def process_structure(junction):

    # get material parameters from the junction and convert them to format required by Sesame (dictionary)
    # convert from Solcore units (i.e. base SI units like m) to Sesame's units (cm, eV).
    # Note that internally in sesame, many quantities are scaled/dimensionless.

    if hasattr(junction, "constant_doping"):
        constant_doping = junction.constant_doping
        # if false, expect only a single material, rather than "emitter" and "base", with a
        # pre-defined doping profile

    else:
        constant_doping = True

    material_list = []
    layer_width = []
    doping_list = []

    for layer in junction:

        try:
            bulk_recombination_energy = layer.material.bulk_recombination_energy

        except:
            bulk_recombination_energy = 0

        layer_width.append(layer.width*1e2) # m to cm

        doping_list.append(layer.material.Nd / 1e6 if layer.material.Nd > layer.material.Na else -layer.material.Na /1e6)

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
        # Note: all of these can be functions of position; should eventually give users a way to use this
        # functionality

        material_list.append(new_mat)

    if not hasattr(junction, "mesh"):
        print("set mesh")
        junction.mesh = np.linspace(0, np.sum(layer_width), 1001) # keep this in m
        # needs to be improved

    junction.mesh_cm = junction.mesh*1e2 # m to cm

    junction.sesame_sys = sesame.Builder(junction.mesh_cm) # Sesame system

    edges = np.insert(np.cumsum(layer_width), 0, 0)
    edges[-1] = edges[-1] + 1e-10 # otherwise final point will not be assigned any values

    for i1 in range(len(junction)):
        junction.sesame_sys.add_material(material_list[i1],
                                         lambda x: np.all((x >= edges[i1], x < edges[i1 + 1]), axis=0))

    if constant_doping:
        doping_profile = np.zeros_like(junction.mesh_cm)
        for i1, dop in enumerate(doping_list):
            doping_profile[np.all((junction.mesh_cm >= edges[i1], junction.mesh_cm <= edges[i1 + 1]), axis=0)] = dop

    else:
        doping_profile = junction.doping_profile

    junction.sesame_sys.rho = doping_profile / junction.sesame_sys.scaling.density

    junction.sesame_sys.contact_type('Ohmic', 'Ohmic')
    # should also be able to choose other options, since Sesame can handle them

    # lower than 1e-6: diverges
    junction.sesame_sys.contact_S(*get_srv(junction))


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
    print('hello')

    # TODO: inclusion of shunt resistance
    if ~hasattr(junction, "sesame_sys"):
        process_structure(junction)

    # gen_wl = profile_func


    if options.light_iv:

        gen_wl = junction.absorbed(junction.mesh) / 100 # m-1 -> cm-1
        wls = options.wavelength

        gg = options.light_source.spectrum(wls, output_units="photon_flux_per_m")[1][:, None] * gen_wl.T

        g_vs_z = np.trapz(gg, wls, axis=0) / 1e4 # m^2 -> cm^2
        g_vs_z[np.isnan(g_vs_z)] = 0
        print(g_vs_z.shape)
        junction.sesame_sys.generation(g_vs_z)


    # equilibrium solution:
    # j0, result0 = sesame.IVcurve(junction.sesame_sys, [0])

    voltages = options.voltages

    j, result = sesame.IVcurve(junction.sesame_sys, voltages)  # , verbose=False)
    j = j * junction.sesame_sys.scaling.current * 1e4 # cm-2 -> m-2
    print(result['v'])

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
        process_structure(junction)

    wls = options.wavelength

    eqe = np.zeros_like(wls)

    voltages = [0]

    # profile_func = interp1d(bulk_positions_cm, 1e7 * Si_profile, kind='linear', bounds_error=False, fill_value=0)
    profile_func = junction.absorbed

    A = np.trapz(junction.absorbed, junction.mesh, axis=0)
    print(A)
    print(junction.layer_absorption)

    def make_gfcn_fun(wl_index, flux):
        def gcfn_fun(x, y):
            return flux * profile_func(x/100)[wl_index]/100 # convert to cm-1 from m-1

        return gcfn_fun

    for i1, wl in enumerate(wls):

        flux = 1e20
        print(wl)
        junction.sesame_sys.generation(make_gfcn_fun(i1, flux))
        print(make_gfcn_fun(i1, flux)(0, 0))
        print(make_gfcn_fun(i1, flux)(1e-6, 0))

        if i1 == 0:
            guess = None

        else:
            guess = {key: result[key][0, :] for key in result.keys()}

        j, result = sesame.IVcurve(junction.sesame_sys, voltages, guess=guess)

        eqe[i1] = np.abs(j) / (q * flux)
        print('j, q, flux', j, q, flux)

    eqe = eqe * junction.sesame_sys.scaling.current
    iqe = eqe / A

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



