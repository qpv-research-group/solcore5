from solcore.solar_cell import SolarCell, Junction, Layer
from solcore.state import State
from solcore.solar_cell_solver import solar_cell_solver
from solcore.sesame_drift_diffusion.solve_pdd import process_structure
import numpy as np
import os
from solcore.light_source import LightSource
from solcore.absorption_calculator import search_db
from solcore.material_system import create_new_material
from scipy.special import erfc
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import seaborn as sns
from cycler import cycler
from solcore.constants import q
from scipy.interpolate import interp1d, RectBivariateSpline

from solcore.absorption_calculator.dielectric_constant_models import DielectricConstantModel, Drude, Cauchy

from rayflare.textures import regular_pyramids #, rough_pyramids
from solcore import material, si
from rayflare.options import default_options
from rayflare.ray_tracing import rt_structure
from rayflare.utilities import make_absorption_function

plt.rc('axes', axisbelow=True)
sns.set_style("whitegrid")

def base_SHJ(additional_generation=None,
             light_source=None,
             voltages=None,
             plot=False,
             surface_recomb=None,
             ):

    rayflare_options = default_options()
    d_bulk = 130e-6
    # shading = 0.02
    # Note Sesame assumes all quantities of length are input in units of cm. Other assumed
    # input units are: time in s, energy in eV.

    # LONGi paper:
    # To further increase JSC, a 150-nm-thick MgF2 film was evaporated on the front TCO layer
    # as a second antireflective coating. For 26.81% cell, an additional 120-nm-thick MgF2/150-nm-thick
    # Ag stack was evaporated on the rear TCO layer, which means this cell is a monofacial solar cell.

    # 130 um Si, n-type 1.2–1.5 Ω cm (~ 1e15-1e16 bulk doping)

    # Jsc = 41.16 mA/cm2
    # optics first (ray-tracing)

    aSi_i = material("aSi_i")()
    aSi_p = material("aSi_p")()
    aSi_n = material("aSi_n")()
    Air = material("Air")()
    MgF2_pageid = str(search_db(os.path.join("MgF2", "Rodriguez-de Marcos"))[0][0])
    MgF2 = material(MgF2_pageid, nk_db=True)()
    Ag = material("Ag_Jiang")()

    # An = 1.6
    # Bn = 0.07
    # Ak = 0.2
    # Ck = 4
    # Bk = 3
    # Cn = 0

    An = 2.1
    Bn = 0.006
    Ak = 0.19
    Ck = 4
    Bk = 3.5
    Cn = 0

    osc = [
        Cauchy(An=An, Bn=Bn, Cn=Cn,
               Ak=Ak, Ck=Ck, Bk=Bk)
    ]

    Output = DielectricConstantModel(oscillators=osc, e_inf=0)

    Si_pn = material("Si")(Nc=3e25, Nv=1e25, electron_mobility=si("1e4cm2"), hole_mobility=si("1e3cm2"),
                           electron_minority_lifetime=0.01, hole_minority_lifetime=0.01,
                           radiative_recombination=si("1.89e-15cm3"),
                           electron_auger_recombination=si("3e-31cm6"), hole_auger_recombination=si("1e-31cm6"))

    wavelengths = np.linspace(280, 1200, 160)*1e-9

    wavelengths_iv = np.arange(280, 1200.01, 0.5)*1e-9

    eps = Output.dielectric_constants(wavelengths*1e9)

    nk = np.sqrt(eps)

    np.savetxt('results/TCO_n.txt', np.array([wavelengths, nk.real]).T)
    np.savetxt('results/TCO_k.txt', np.array([wavelengths, nk.imag]).T)

    create_new_material('TCO_SHJ', 'results/TCO_n.txt', 'results/TCO_k.txt',
                        overwrite=True)

    TCO = material('TCO_SHJ')()

    front_materials = [Layer(150e-9, MgF2), Layer(55e-9, TCO),
                       Layer(5e-9, aSi_i), Layer(1e-9, aSi_p)]
    back_materials = [Layer(2e-9, aSi_n), Layer(2e-9, aSi_i), Layer(55e-9, TCO),
                      Layer(120e-9, MgF2), Layer(150e-9, Ag)]

    triangle_surf = regular_pyramids(elevation_angle=52, upright=True, size=1, interface_layers=front_materials)

    triangle_surf_rear = regular_pyramids(52, upright=False, size=1, interface_layers=back_materials)

    rayflare_options.wavelengths = wavelengths
    rayflare_options.depth_spacing_bulk = 10e-9
    rayflare_options.project_name = 'SF_whitepaper'
    rayflare_options.n_rays = 7000
    rayflare_options.randomize_surface = True

    # set up ray-tracing options
    rtstr = rt_structure(
        textures=[triangle_surf, triangle_surf_rear],
        materials=[Si_pn],
        widths=[d_bulk],
        incidence=Air,
        transmission=Air,
        use_TMM=True,
        options=rayflare_options,
        save_location="current",
        overwrite=True,
    )

    rt_filename = f'results/SHJ_optics{rayflare_options.n_rays}.npy'
    if os.path.exists(rt_filename):
        result_rt = np.load(rt_filename, allow_pickle=True)[()]

    else:
        result_rt = rtstr.calculate_profile(rayflare_options)
        np.save(rt_filename, result_rt)

    pal = sns.color_palette('pastel6', 5)
    cols = cycler("color", pal)

    params = {"axes.prop_cycle": cols}

    plt.rcParams.update(params)

    # fig, ax = plt.subplots(1, 1, figsize=(5, 4))
    #
    result_stack = [
        result_rt['A_per_layer'][:,0],
        np.sum(result_rt['A_per_interface'][0], 1),
        np.sum(result_rt['A_per_interface'][1], 1),
        result_rt['R0'],
        result_rt['R'] - result_rt['R0'],
    ]

    positions_Si = np.arange(0, d_bulk, rayflare_options.depth_spacing_bulk)

    rayflare_options.wavelengths = wavelengths_iv

    profile_data = result_rt['profile']

    interp_spline = RectBivariateSpline(wavelengths, positions_Si, profile_data)

    interp_profile = interp_spline(wavelengths_iv, positions_Si)

    result_rt_interp = {'profile': interp_profile}

    R_interp_f = interp1d(wavelengths, result_rt['R'], bounds_error=False, fill_value=0)

    R_interp = R_interp_f(wavelengths_iv)

    positions, diff_absorb_fn = make_absorption_function(
        result_rt_interp, rtstr, rayflare_options,
    )

    options = State()
    options.wavelength = wavelengths_iv
    options.optics_method = 'external'
    options.light_iv = True
    options.T = 298

    if light_source is None:
        options.light_source = LightSource(source_type="standard",
                               version="AM1.5g", x=options.wavelength, output_units="photon_flux_per_m")

    else:
        options.light_source = light_source


    if voltages is None:
        options.voltages = np.linspace(0, 0.8, 50)

    else:
        options.voltages = voltages

    options.mpp = True
    options.recalculate_absorption = True

    # ray-tracing for optics first

    limiting_current = q*np.trapz(result_rt['A_per_layer'][:,0]*
                                  options.light_source.spectrum(x=wavelengths,
                                                                output_units="photon_flux_per_m")[1], wavelengths)

    print('Limiting current: {:.2f} mA/cm2'.format(limiting_current/10))

    front_width = np.sum([layer.width for layer in front_materials])
    back_width = np.sum([layer.width for layer in back_materials])

    # L = 130e-6

    # rear junction

    nD = si("1e20cm-3")
    nA = si("1e19cm-3")
    bulk_doping = si("5e15cm-3") # n type bulk

    # rear junction (n-type)
    def doping_profile_func(x):

        L = d_bulk

        doping_profile = - nA * erfc(x/150e-9) # characteristic depth of 150 nm

        doping_profile_rear = nD * erfc((L - x)/200e-9) # characteristic depth of 200 nm

        return doping_profile + doping_profile_rear + bulk_doping

    L_cm = d_bulk*100 # length of the system in the x-direction [cm]

    # Mesh
    x = np.concatenate((np.arange(0, 200e-7-1e-10, 0.25e-7),
                        np.linspace(200e-7,2e-4, 1000, endpoint=False),
                        np.linspace(2e-4, L_cm - 2e-4, 2000, endpoint=False),
                        np.linspace(L_cm - 1e-4, L_cm, 500)))

    options.position = np.concatenate((np.linspace(0, front_width, 100), front_width + x/100,
                                      np.linspace(front_width + d_bulk,
                                                  front_width + d_bulk + back_width, 100)))

    options.minimum_spacing = 5e-9

    if surface_recomb is None:
        surface_recomb = dict(sn_front=si('5e3 cm s-1'),  # important
                              sp_front=si('1e4 cm s-1'),
                              sn_rear=si('1e4 cm s-1'),
                              sp_rear=si('1e4 cm s-1'))  # important

    Si_junction = [Junction([Layer(d_bulk, Si_pn)],
                            doping_profile=doping_profile_func, kind='sesame_PDD',
                            mesh=x/100,
                            **surface_recomb,
                            )
                   ]

    process_structure(Si_junction[0], options)

    if additional_generation is not None:

        interp_profile_gen = interp_spline(wavelengths_iv, x/100)

        gen_wl = interp_profile_gen * 1e9 / 100  # nm-1 -> m-1 -> cm-1
        wls = options.wavelength

        gg = options.light_source.spectrum(wls, output_units="photon_flux_per_m")[1][:, None] * gen_wl

        g_vs_z = np.trapz(gg, wls, axis=0) / 1e4  # m^2 -> cm^2
        g_vs_z[np.isnan(g_vs_z)] = 0

        extra_generation = 1e-6*additional_generation(x/100)
        # this returns a result in terms of photons m-2 (area) m-1 (depth). Need to multiply by 1e-6 to convert to cm

        g_vs_z += extra_generation

        # can also pass a function to generation - more flexible?
        Si_junction[0].sesame_sys.generation(g_vs_z)

    Si_cell = SolarCell(front_materials +
                         Si_junction +
                         back_materials,
                        external_reflected=R_interp,
                        external_absorbed=diff_absorb_fn,
                        )

    Si_ind = len(front_materials)

    solar_cell_solver(Si_cell, 'iv', options)

    if plot:

        rayflare_options.wavelengths = wavelengths
        options.wavelength = wavelengths

        positions_2, diff_absorb_fn_2 = make_absorption_function(
            result_rt, rtstr, rayflare_options,
        )

        Si_cell_EQE = SolarCell(front_materials +
                            Si_junction +
                            back_materials,
                            external_reflected=result_rt['R'],
                            external_absorbed=diff_absorb_fn_2,
                            )

        solar_cell_solver(Si_cell_EQE, 'qe', options)

        shading = 0.02

        sns.set_style('ticks')
        fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10, 3.5))
        ax.stackplot(wavelengths * 1e9, 100*(1-shading)*np.array(result_stack), linewidth=0.5)
        ax.fill_between(wavelengths*1e9, 100*(1-shading), 100, color='white')
        ax.plot(wavelengths* 1e9, 100*(1-shading)*Si_cell_EQE[Si_ind].eqe(wavelengths), '-r', linewidth=2)
        ax.legend(['Si absorption', 'Front stack absorption', 'Rear stack aborption',
                   'Front reflection', 'Front escape', 'Electrode shading', 'Simulated EQE'],
                  loc='lower center')
        ax.set_xlim(280, 1200)
        ax.set_ylim(0, 100)
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("R / A / EQE (%)")
        ax.set_title('a) EQE and cell optics', loc='left')
        # plt.show()


        jsc = Si_cell.iv.Isc / 10

        ax2.plot(options.voltages, -(1 - shading) * Si_cell[Si_ind].iv(options.voltages) / 10, '-', label='IV',
                 linewidth=2, color='k')
        # plt.plot(options.voltages, solar_cell_old[1].iv(options.voltages), 'o-', label='IV')
        # plt.xlim(-2, 1.8)

        ax2.set_ylim(0, 1.03 * jsc)
        ax2.set_xlim(np.min(options.voltages), np.max(options.voltages))
        ax2.set_xlabel('Voltage (V)')
        ax2.set_ylabel('Current density (mA/cm$^2$)')
        ax2.set_title('b) IV characteristics and power output', loc='left')

        ax3 = ax2.twinx()
        ax3.plot(options.voltages, Si_cell.iv['IV'][0] * Si_cell.iv['IV'][1],
                 '-', color=pal[2], label='Power', linewidth=2)
        ax3.set_ylabel('Power density (W m$^{-2}$)')
        ax3.set_ylim(0, 1.03 * jsc * 10)


        ax3.spines['right'].set_color(pal[2])
        ax3.yaxis.label.set_color(pal[2])
        ax3.tick_params(axis='y', colors=pal[2])

        ax2.set_axisbelow(True)
        ax3.set_axisbelow(True)

        ax2.text(0.02, 0.9 * jsc, r'$J_{SC}$', zorder=5)
        ax2.text(0.02, 0.8 * jsc, r'$V_{OC}$')
        ax2.text(0.02, 0.7 * jsc, 'FF')
        ax2.text(0.02, 0.6 * jsc, r'$\eta$')
        ax2.text(0.02, 0.5 * jsc, r'$J_{MPP}$')
        ax2.text(0.02, 0.4 * jsc, r'$V_{MPP}$')

        ax2.text(0.1, 0.9 * jsc, r'= {:.2f} mA/cm$^2$'.format(jsc))
        ax2.text(0.1, 0.8 * jsc, r'= {:.3f} V'.format(Si_cell.iv.Voc))
        ax2.text(0.1, 0.7 * jsc, '= {:.2f} %'.format(Si_cell.iv.FF * 100))
        ax2.text(0.1, 0.6 * jsc, r'= {:.2f} %'.format((1 - shading) * Si_cell.iv.Pmpp/10))
        ax2.text(0.1, 0.5 * jsc, r'= {:.2f} mA/cm$^2$'.format((1 - shading) * Si_cell.iv.Impp / 10))
        ax2.text(0.1, 0.4 * jsc, r'= {:.3f} V'.format(Si_cell.iv.Vmpp))
        ax2.grid(which='major', alpha=0.35)
        ax2.add_patch(Rectangle((0.01, 0.38 * jsc), 0.31, 0.58 * jsc, facecolor='white', edgecolor='lightgray',
                                alpha=0.8))

        ax3.grid(False)
        plt.tight_layout()
        fig.savefig('SHJ_plot.pdf', bbox_inches='tight')
        plt.show()

    return Si_cell, Si_ind, result_rt

if __name__ == "__main__":
    surface_recomb = dict(sn_front=si('1e4 cm s-1'), # important
                          sp_front=si('1e4 cm s-1'),
                          sn_rear=si('1e4 cm s-1'),
                          sp_rear=si('1e4 cm s-1')) # important

    Si_cell, Si_ind, result_rt = base_SHJ(plot=True,
                                          surface_recomb=surface_recomb,)