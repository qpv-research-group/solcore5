import numpy as np

import solcore.poisson_drift_diffusion as PDD
import solcore.analytic_solar_cells as ASC
from solcore.light_source import LightSource
from solcore.state import State
from solcore.optics import solve_beer_lambert

default_options = State()

# General
default_options.T_ambient = 298
default_options.T = 298

# Illumination spectrum
default_options.wavelength = np.linspace(300, 1300, 251)
default_options.light_source = LightSource(source_type='standard', version='AM1.5g', x=default_options.wavelength,
                                           output_units='photon_flux_per_nm')

# IV control
default_options.voltages = np.linspace(0, 1.2, 100)
default_options.mpp = False
default_options.light_iv = False
default_options.db_mode = 'top_hat'
default_options.E_range = 0.3
default_options.internal_voltages = np.linspace(-2.5, 2.5, 300)
default_options.DA_print_parameters = False

# Optics control
default_options.optics_method = 'BL'


def merge_dicts(*dict_args):
    """
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    """
    result = State()
    for dictionary in dict_args:
        result.update(dictionary)
    return result


def solar_cell_solver(solar_cell, task, user_options=None):
    """ Solves the properties of a solar cell object, either calculating its optical properties (R, A and T), its quantum efficiency or its current voltage characteristics in the dark or under illumination. The general options for the solvers are passed as dicionaries.

    :param solar_cell: A solar_cell object
    :param task: Task to perform. It has to be "optics", "dark_iv" or "qe"
    :param solver_options: General solver options
    :param optics_options: Options specific for the optics solver
    :param iv_options: Options specific for the IV solver.
    :param qe_options: Options specific for the QE solver.
    :return:
    """
    if type(user_options) in [State, dict]:
        options = merge_dicts(default_options, user_options)
    else:
        options = merge_dicts(default_options)

    if task == 'optics':
        solve_optics(solar_cell, options)
    elif task == 'iv':
        solve_iv(solar_cell, options)
    elif task == 'qe':
        solve_qe(solar_cell, options)
    else:
        raise ValueError(
            'ERROR in "solar_cell_solver":\n\tValid values for "task" are: "optics", "iv" and "qe".')


def solve_optics(solar_cell, options):
    """ Solves the optical properties of the structure, calculating the reflectance, absorptance and transmitance. The "optics_method" option controls which method is used to calculate the optical properties of the solar cell:

    - None: The calculation is skipped. Only useful for solar cells involving just "2-diode" kind of junctions.
    - BL: Uses the Beer-Lambert law to calculate the absorption in each layer. Front surface reflexion has to provided externally. It is the default method and the most flexible one.
    - TMM: Uses a transfer matrix calculation to obtain the RAT.
    - RWCA: Uses the rigorous wave coupled analysisto obtain the RAT. This allows to include 2D photonic crystals in the structure, for example.

    :param solar_cell: A solar_cell object
    :param options: Options for the optics solver
    :return: None
    """
    if options.optics_method is None:
        print('Warning: Not solving the optics of the solar cell.')
    elif options.optics_method == 'BL':
        solve_beer_lambert(solar_cell, options)
    else:
        raise ValueError(
            'ERROR in "solar_cell_solver":\n\tOptics solver method must be None, "BL", "TMM" or "RWCA".')


def solve_iv(solar_cell, options):
    """ Calculates the IV at a given voltage range, providing the IVs of the individual junctions in addition to the total QE

    :param solar_cell: A solar_cell object
    :param iv_options: Options for the iv solver
    :param optics_options: Options for the optics solver
    :return: None
    """

    solve_optics(solar_cell, options)

    for j in solar_cell.junction_indices:

        if solar_cell[j].kind == 'PDD':
            PDD.iv_pdd(solar_cell[j], options)
        elif solar_cell[j].kind == 'DA':
            ASC.iv_depletion(solar_cell[j], options)
        elif solar_cell[j].kind == '2D':
            ASC.iv_2diode(solar_cell[j], options)
        elif solar_cell[j].kind == 'DB':
            ASC.iv_detailed_balance(solar_cell[j], options)
        else:
            raise ValueError(
                'ERROR in "solar_cell_solver":\n\tJunction {} has an invalid "type". It must be "PDD", "DA", "2D" or "DB".'.format(
                    j))

    ASC.iv_multijunction(solar_cell, options)


def solve_qe(solar_cell, options):
    """ Calculates the QE of all the junctions

    :param solar_cell: A solar_cell object
    :param qe_options: Options for the QE solver
    :param optics_options: Options for the optics solver
    :return: None
    """

    solve_optics(solar_cell, options)

    for j in solar_cell.junction_indices:

        if solar_cell[j].kind is 'PDD':
            PDD.qe_pdd(solar_cell[j], options)
        elif solar_cell[j].kind is 'DA':
            ASC.qe_depletion(solar_cell[j], options)
        elif solar_cell[j].kind is '2D':
            pass
        elif solar_cell[j].kind is 'DB':
            ASC.qe_detailed_balance(solar_cell[j])
        else:
            raise ValueError(
                'ERROR in "solar_cell_solver":\n\tJunction {} has an invalid "type". It must be "PDD", "DA", "2D" or "DB".'.format(
                    j))

    ASC.qe_multijunction(solar_cell, options)


if __name__ == '__main__':
    from solcore.structure import Junction
    from solcore.solar_cell import SolarCell, default_GaAs
    import matplotlib.pyplot as plt
    from solcore.structure import Layer
    from solcore import si

    T = 298

    ##### Input data for the 2D kind of junction
    # From the QE object we get the short circuit currents
    Isc_array = [250, 150, 200]
    I01_array = [4.93e-24, 1.0e-21, 4.93e-6]
    I02_array = [3.28e-15, 2.7e-10, 1.0e-5]

    # This is the structure to calculate.
    # db_junction = Junction(kind='2D', Eg=0.66, j01=I01_array[2], j02=I02_array[2], R_shunt=115,
    #                        n1=1.00, n2=2.0, jsc=Isc_array[2])
    # db_junction2 = Junction(kind='2D', Eg=1.4, j01=I01_array[1], j02=I02_array[1], R_shunt=1.5e6,
    #                         n1=1.00, n2=2.0, jsc=Isc_array[1])
    # db_junction3 = Junction(kind='2D', Eg=1.9, j01=I01_array[0], j02=I02_array[0], R_shunt=3e6,
    #                         n1=1.00, n2=2.0, jsc=Isc_array[0])

    ##### Input data for the DB kind of junction
    db_junction = Junction(kind='DB', T=T, Eg=0.77, A=1, R_shunt=np.inf, n=3.5)
    db_junction2 = Junction(kind='DB', T=T, Eg=1.4, A=1, R_shunt=np.inf, n=3.5)
    db_junction3 = Junction(kind='DB', T=T, Eg=1.9, A=1, R_shunt=np.inf, n=3.5)

    ##### Input data for the DA kind of junction
    from solcore.examples.ASC_examples.ASC_example_1 import window_material, top_cell_n_material, top_cell_p_material, \
        mid_cell_n_material, mid_cell_p_material, bot_cell_n_material, bot_cell_p_material, ref
    #
    # window = Layer(material=window_material, width=si("25nm"))
    # tJ = Layer(material=bot_cell_n_material, width=si("25nm"))
    # db_junction3 = Junction([Layer(si("100nm"), material=top_cell_n_material, role='emitter'),
    #                          Layer(si("600nm"), material=top_cell_p_material, role='base')], kind='DA', sn=1, sp=1)
    # db_junction2 = Junction([Layer(si("100nm"), material=mid_cell_n_material, role='emitter'),
    #                          Layer(si("3.5um"), material=mid_cell_p_material, role='base')], kind='DA', sn=1, sp=1)
    # db_junction = Junction([Layer(si("400nm"), material=bot_cell_n_material, role='emitter'),
    #                         Layer(si("100um"), material=bot_cell_p_material, role='base')], kind='DA', sn=1, sp=1)

    #
    # ##### Input data for PDD kind of junction
    # from solcore import material
    #
    # def default_GaAsP(T):
    #     # We create the other materials we need for the device
    #     window = material('AlGaAs')(T=T, Na=5e24, Al=0.8)
    #     p_GaAs = material('AlGaAs')(T=T, Na=1e24, Al=0.4)
    #     n_GaAs = material('AlGaAs')(T=T, Nd=8e22, Al=0.4)
    #     bsf = material('AlGaAs')(T=T, Nd=2e24, Al=0.4)
    #
    #     output = Junction([Layer(width=si('30nm'), material=window, role="Window"),
    #                        Layer(width=si('150nm'), material=p_GaAs, role="Emitter"),
    #                        Layer(width=si('600nm'), material=n_GaAs, role="Base"),
    #                        Layer(width=si('200nm'), material=bsf, role="BSF")], sn=1e6, sp=1e6, T=T, kind='PDD')
    #
    #     return output

    ##### Common comands
    my_solar_cell = SolarCell([db_junction3, db_junction2, db_junction], T=T, R_series=0)
    my_solar_cell.reflectivity = ref
    # my_solar_cell = SolarCell([db_junction2], T=T, R_series=0)

    # my_solar_cell = SolarCell([default_GaAsP(T), default_GaAs(T)], T=T, R_series=0)


    V = np.linspace(0, 3.5, 100)
    wl = np.linspace(300, 1900, 801)
    z = np.linspace(0, 1, 1001)

    solar_cell_solver(my_solar_cell, 'iv',
                      user_options={'T_ambient': T, 'db_mode': 'top_hat', 'voltages': V, 'light_iv': True,
                                    'internal_voltages': np.linspace(-6, 2.5, 400), 'wavelength': wl})

    plt.figure()
    # for j in my_solar_cell.junction_indices:
    #     plt.plot(wl, my_solar_cell(j).absorptance(wl))

    for data in my_solar_cell.iv['junction IV']:
        plt.plot(data[0], abs(data[1]), 'o')
    # #
    # # for j in my_solar_cell.junction_indices:
    # #     plt.semilogy(V, abs(my_solar_cell(j).iv(V)))
    #
    plt.plot(my_solar_cell.iv['IV'][0], abs(my_solar_cell.iv['IV'][1]), 'k', linewidth=2)
    # plt.semilogy(abs(my_solar_cell.iv['Rseries IV'][0]), abs(my_solar_cell.iv['Rseries IV'][1]))

    # for j in range(len(my_solar_cell)):
    #     absorbed = my_solar_cell[j].absorbed(z*my_solar_cell[j].width)
    #     abs_vs_wl = np.trapz(absorbed, z*my_solar_cell[j].width, axis=0)
    #     # abs_vs_z = np.trapz(absorbed, wl, axis=1)
    #
    #     # plt.contourf(wl*1e9, z, absorbed, 50, cmap='jet')
    #     # plt.plot(z, abs_vs_z)
    #     plt.plot(wl, abs_vs_wl, label='junction {}'.format(j))

    # absorbed = my_solar_cell.absorbed(z * 1e-9)
    # # plt.contourf(wl, z, np.log(absorbed), 50, cmap='jet')
    # #
    # abs_vs_wl = np.trapz(absorbed, z * 1e-9, axis=0)
    # plt.plot(wl, abs_vs_wl, 'k', linewidth=2)

    # abs_vs_z = np.trapz(absorbed, wl, axis=1)
    # plt.plot(z, abs_vs_z)

    # plt.plot(wl, my_solar_cell.absorbed, label='Absorbed')
    # plt.plot(wl, my_solar_cell.transmitted, label='Transmitted')
    # plt.plot(wl, my_solar_cell.reflected, label='Reflected')
    plt.legend()

    plt.show()
