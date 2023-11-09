import numpy as np
from pytest import approx

def test_constant_doping():
    # test sesame PDD process_structure when constant doping levels are passed
    from solcore import material, si
    from solcore.structure import Junction, Layer
    from solcore.sesame_drift_diffusion.solve_pdd import process_structure
    from solcore.state import State

    Si_n = material('Si')(Nd=1e24, electron_minority_lifetime=1e-6, hole_minority_lifetime=1e-7)
    Si_i = material('Si')(Na=0, Nd=0, electron_minority_lifetime=2e-6, hole_minority_lifetime=2e-7)
    Si_p = material('Si')(Na=2e24, electron_minority_lifetime=1.5e-6, hole_minority_lifetime=1.5e-7)

    junction = Junction([Layer(si('200nm'), Si_n),
                            Layer(si('2000nm'), Si_i),
                            Layer(si('2000nm'), Si_p)])

    options = State()
    options.T = 300

    process_structure(junction, options)

    sesame_obj = junction.sesame_sys

    assert len(np.unique(sesame_obj.rho)) == 3
    assert sesame_obj.rho[0] == 1e24*1e-6/sesame_obj.scaling.density
    assert sesame_obj.rho[-1] == -2e24*1e-6/sesame_obj.scaling.density

    pass

def test_layer_doping_profile():
    # test process_structure when a doping profile is passed for one of the layers
    pass

def test_junction_doping_profile():
    # test process_structure when a doping profile is passed for the whole junction
    pass

def test_parameter_extraction():
    # test all material parameters are correctly extracted by process_structure
    pass

def test_dark_IV():
    pass

def test_light_IV():
    pass

def test_multijunction_IV():
    pass

def test_compare_DA_np():

    from solcore import material, si
    from solcore.structure import Junction, Layer
    from solcore.solar_cell_solver import solar_cell_solver, SolarCell
    from solcore.analytic_solar_cells import iv_depletion
    from solcore.sesame_drift_diffusion import iv_sesame
    from solcore.state import State
    from solcore.light_source import LightSource

    GaAs_n = material('GaAs')(T=300, Nd=1e24,
                              hole_minority_lifetime=1e-6, electron_minority_lifetime=1e-6)
    GaAs_p = material('GaAs')(T=300, Na=1e24,
                              hole_minority_lifetime=1e-6, electron_minority_lifetime=1e-6)

    options = State()
    options.wavelength = np.linspace(300, 950, 100)*1e-9
    options.voltages = np.linspace(-1.3, 0.5, 30)
    options.internal_voltages = np.linspace(-1.3, 0.5, 30)
    options.T = 300
    options.light_iv = True
    options.light_source = LightSource(source_type='standard', version='AM1.5g', x=options.wavelength, output_units='photon_flux_per_m')
    options.da_mode = 'green'
    options.optics_method = 'TMM'

    mesh = np.linspace(0, 2200, 500)*1e-9

    pn_junction = Junction([Layer(si('200nm'), GaAs_n, role='emitter'), Layer(si('2000nm'), GaAs_p, role='base')],
                           kind='DA', R_shunt=0.1, mesh=mesh)

    solar_cell_solver(SolarCell([pn_junction]), 'optics', options)

    iv_depletion(pn_junction, options)
    depletion_current = pn_junction.current
    depletion_current_interp = pn_junction.iv(options.voltages)

    iv_sesame(pn_junction, options)
    sesame_current = pn_junction.current
    sesame_current_interp = pn_junction.iv(options.voltages)

    assert depletion_current[-1] == approx(sesame_current[-1], rel=0.05)
    assert np.sign(depletion_current[0]) == np.sign(sesame_current[0])

    assert depletion_current_interp[-1] == approx(sesame_current_interp[-1], rel=0.05)
    assert np.sign(depletion_current_interp[0]) == np.sign(sesame_current_interp[0])

    assert np.all(sesame_current == sesame_current_interp)


def test_compare_DA_pn():

    from solcore import material, si
    from solcore.structure import Junction, Layer
    from solcore.solar_cell_solver import solar_cell_solver, SolarCell
    from solcore.analytic_solar_cells import iv_depletion
    from solcore.sesame_drift_diffusion import iv_sesame
    from solcore.state import State
    from solcore.light_source import LightSource

    GaAs_n = material('GaAs')(T=300, Na=1e24,
                              hole_minority_lifetime=1e-6, electron_minority_lifetime=1e-6)
    GaAs_p = material('GaAs')(T=300, Nd=1e24,
                              hole_minority_lifetime=1e-6, electron_minority_lifetime=1e-6)

    options = State()
    options.wavelength = np.linspace(300, 950, 100)*1e-9
    options.voltages = np.linspace(-0.5, 1.3, 30)
    options.internal_voltages = np.linspace(-0.5, 1.3, 30)
    options.T = 300
    options.light_iv = True
    options.light_source = LightSource(source_type='standard', version='AM1.5g', x=options.wavelength, output_units='photon_flux_per_m')
    options.da_mode = 'green'
    options.optics_method = 'TMM'

    mesh = np.linspace(0, 2200, 500)*1e-9

    pn_junction = Junction([Layer(si('200nm'), GaAs_n, role='emitter'), Layer(si('2000nm'), GaAs_p, role='base')],
                           kind='DA', R_shunt=0.1, mesh=mesh)

    solar_cell_solver(SolarCell([pn_junction]), 'optics', options)

    iv_depletion(pn_junction, options)
    depletion_current = pn_junction.current
    depletion_current_interp = pn_junction.iv(options.voltages)

    iv_sesame(pn_junction, options)
    sesame_current = pn_junction.current
    sesame_current_interp = pn_junction.iv(options.voltages)

    assert depletion_current[0] == approx(sesame_current[0], rel=0.05)
    assert np.sign(depletion_current[-1]) == np.sign(sesame_current[-1])

    assert depletion_current_interp[0] == approx(sesame_current_interp[0], rel=0.05)
    assert np.sign(depletion_current_interp[-1]) == np.sign(sesame_current_interp[-1])

    assert np.all(sesame_current == sesame_current_interp)
