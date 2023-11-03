

def test_constant_doping():



    pass

def test_layer_doping_profile():
    pass

def test_junction_doping_profile():
    pass

def test_combination_doping():
    pass


def test_parameter_extraction():
    pass

def test_dark_IV():
    pass

def test_light_IV():
    pass


def test_multijunction_IV():
    pass


def test_compare_DA():

    from solcore import material, si
    from solcore.structure import Junction, Layer
    from solcore.solar_cell_solver import solar_cell_solver, SolarCell
    from solcore.analytic_solar_cells import iv_depletion
    from solcore.sesame_drift_diffusion import iv_sesame
    from solcore.state import State
    from solcore.optics import solve_tmm
    from solcore.light_source import LightSource
    import numpy as np

    GaAs_n = material('GaAs')(T=300, Nd=1e24, hole_minority_lifetime=1e-6, electron_minority_lifetime=1e-6)
    GaAs_p = material('GaAs')(T=300, Na=1e24, hole_minority_lifetime=1e-6, electron_minority_lifetime=1e-6)

    options = State()
    options.wavelength = np.linspace(300, 950, 100)*1e-9
    options.voltages = np.linspace(0, 1.3, 10)
    options.internal_voltages = np.linspace(0, 1.3, 10)
    options.T = 300
    options.light_iv = True
    options.light_source = LightSource(source_type='standard', version='AM1.5g', x=options.wavelength, output_units='photon_flux_per_m')
    options.da_mode = 'green'
    options.optics_method = 'TMM'

    pn_junction = Junction([Layer(si('200nm'), GaAs_n, role='emitter'), Layer(si('2000nm'), GaAs_p, role='base')],
                           kind='DA')

    solar_cell_solver(SolarCell([pn_junction]), 'optics', options)

    import matplotlib.pyplot as plt
    plt.figure()
    iv_depletion(pn_junction, options)
    plt.plot(options.voltages, pn_junction.iv(options.voltages), label='DA')
    iv_sesame(pn_junction, options)
    plt.plot(options.voltages, pn_junction.iv(options.voltages), label='Sesame')
    plt.ylim(-250, 250)
    plt.show()

