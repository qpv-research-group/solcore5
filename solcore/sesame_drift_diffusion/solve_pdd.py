from sesame import Builder, IVcurve, Analyzer
import numpy as np
from solcore.constants import q
from solcore.light_source import LightSource
from solcore.interpolate import interp1d
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
    j = j * junction.sesame_sys.scaling.current

    Jsc = j[0]*1e4 # units?

    # jsc = q*np.trapz(eqe*light_source.spectrum()[1], wls) # A/m2

    zero_crossing = np.where(np.diff(np.sign(j)))[0][0]
    j_above = j[zero_crossing]
    j_below = j[zero_crossing + 1]

    Voc = voltages[zero_crossing] + (voltages[zero_crossing + 1] - voltages[zero_crossing]) * j_above / (
                j_above - j_below)

    voltages_for_mpp = voltages * (np.abs(voltages) <= Voc)
    Vmpp = voltages[np.nanargmax(np.abs(j * voltages_for_mpp))]
    Jmpp = j[np.nanargmax(np.abs(j * voltages_for_mpp))] * 1e4

    FF = Vmpp * Jmpp / (Jsc * Voc)


def qe_sesame(junction, options):

    if ~hasattr(junction, "sesame_sys"):
        process_structure(junction)

    wls = options.wavelength

    eqe = np.zeros_like(wls)

    voltages = [0]

    for i1, wl in enumerate(wls):

        flux = 1e20
        print(wl)
        junction.sesame_sys.generation(gfcn_fun(i1, flux))

        if i1 == 0:
            guess = None

        else:
            guess = {key: result[key][0, :] for key in result.keys()}

        j, result = sesame.IVcurve(junction.sesame_sys, voltages, guess=guess)

        eqe[i1] = j / (q * flux)

    A = junction.layer_absorption
    iqe = eqe / A

    # convert dimensionless current to dimension-ful current
    eqe = eqe * junction.sesame_sys.scaling.current
    iqe = iqe * junction.sesame_sys.scaling.current

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



from solcore import material, si
from solcore.solar_cell import SolarCell, Junction, Layer
from solcore.state import State
from solcore.solar_cell_solver import solar_cell_solver
from solcore.light_source import LightSource

options = State()
options.wavelength = np.linspace(280, 600, 20)*1e-9
options.optics_method = 'TMM'
options.light_iv = True
options.light_source = LightSource(source_type="standard",
                           version="AM1.5g", x=options.wavelength, output_units="photon_flux_per_m")
options.voltages = np.linspace(0, 2, 100)

T = 293

add_args = {'relative_permittivity': 10, 'electron_minority_lifetime': 5e-6,
            'hole_minority_lifetime': 5e-6,
            'electron_auger_recombination': 1e-45,
            'hole_auger_recombination': 1e-45}

ARC = material('Si3N4')()
window = material('AlGaAs')(T=T, Na=5e24, Al=0.8, **add_args)
p_AlGaAs = material('AlGaAs')(T=T, Na=1e24, Al=0.4, **add_args)
n_AlGaAs = material('AlGaAs')(T=T, Nd=8e22, Al=0.4, **add_args)
bsf = material('AlGaAs')(T=T, Nd=2e24, Al=0.6, **add_args)

junction = Junction([Layer(width=si('30nm'), material=window, role="Window"),
                   Layer(width=si('150nm'), material=p_AlGaAs, role="Emitter"),
                   Layer(width=si('1000nm'), material=n_AlGaAs, role="Base"),
                   Layer(width=si('200nm'), material=bsf, role="BSF")], sn=1e6, sp=1e6, T=T, kind='PDD')

widths = [layer.width for layer in junction]
junction.mesh = np.linspace(0, np.sum(widths), 1000)

solar_cell = SolarCell(
    [Layer(60e-0, ARC), junction]
)

solar_cell_solver(solar_cell, 'optics', options)

# qe_sesame(solar_cell[1], options)

iv_sesame(solar_cell[1], options)