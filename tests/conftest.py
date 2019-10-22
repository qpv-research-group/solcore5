from pytest import fixture
import numpy as np


@fixture
def wavelength():
    return np.linspace(300, 1200)


@fixture
def gauss_spectrum(wavelength):
    from scipy.interpolate import interp1d

    centre = wavelength.mean()
    width = (wavelength.max() - wavelength.min()) / 6
    sp = np.exp(-(wavelength - centre) ** 2 / width)
    return sp, interp1d(wavelength, sp, bounds_error=False, fill_value=0)


@fixture
def np_junction():
    from solcore.structure import Junction, Layer
    from solcore import si, material
    from solcore.solar_cell import SolarCell
    from solcore.constants import vacuum_permittivity
    from solcore.solar_cell_solver import prepare_solar_cell
    from solcore.state import State
    from solcore.light_source import LightSource
    from solcore.optics import solve_beer_lambert

    Nd = np.power(10, np.random.uniform(24, 26))
    Na = np.power(10, np.random.uniform(23, 25))

    Lp = np.power(10, np.random.uniform(-9, -5))
    Ln = np.power(10, np.random.uniform(-9, -5))

    AlInP = material("AlInP")
    InGaP = material("GaInP")
    window_material = AlInP(Al=0.52)
    top_cell_n_material = InGaP(In=0.48, Nd=Nd,
                                hole_diffusion_length=Lp)
    top_cell_p_material = InGaP(In=0.48, Na=Na,
                                electron_diffusion_length=Ln)

    rel_perm = np.random.uniform(1,20)
    for mat in [top_cell_n_material, top_cell_p_material]:
        mat.permittivity = rel_perm * vacuum_permittivity

    n_width = np.random.uniform(500,1000)*1e-9
    p_width = np.random.uniform(3000, 5000)*1e-9

    test_junc = SolarCell([Junction([Layer(si("25nm"), material=window_material, role='window'),
                  Layer(n_width, material=top_cell_n_material, role='emitter'),
                  Layer(p_width, material=top_cell_p_material, role='base'),
                 ], sn=1, sp=1, kind='DA')])

    light_source = LightSource(source_type="standard", version="AM1.5g")

    options = State()
    wl = np.linspace(290, 700, 150) * 1e-9
    options.T = np.random.uniform(10, 350)
    options.wavelength = wl
    options.light_source = light_source
    options.position = None
    options.internal_voltages = np.linspace(-6, 4, 20)
    prepare_solar_cell(test_junc, options)

    solve_beer_lambert(test_junc, options)
    return test_junc, light_source, options



@fixture
def pn_junction():
    from solcore.structure import Junction, Layer
    from solcore import si, material
    from solcore.solar_cell import SolarCell
    from solcore.constants import vacuum_permittivity
    from solcore.solar_cell_solver import prepare_solar_cell
    from solcore.state import State
    from solcore.light_source import LightSource
    from solcore.optics import solve_beer_lambert

    Na = np.power(10, np.random.uniform(24, 26))
    Nd = np.power(10, np.random.uniform(23, 25))

    Lp = np.power(10, np.random.uniform(-9, -5))
    Ln = np.power(10, np.random.uniform(-9, -5))

    AlInP = material("AlInP")
    InGaP = material("GaInP")
    window_material = AlInP(Al=0.52)
    top_cell_n_material = InGaP(In=0.48, Nd=Nd,
                                hole_diffusion_length=Lp)
    top_cell_p_material = InGaP(In=0.48, Na=Na,
                                electron_diffusion_length=Ln)

    rel_perm = np.random.uniform(1,20)
    for mat in [top_cell_n_material, top_cell_p_material]:
        mat.permittivity = rel_perm * vacuum_permittivity

    p_width = np.random.uniform(500,1000)*1e-9
    n_width = np.random.uniform(3000, 5000)*1e-9

    test_junc = SolarCell([Junction([Layer(si("25nm"), material=window_material, role='window'),
                  Layer(p_width, material=top_cell_p_material, role='emitter'),
                  Layer(n_width, material=top_cell_n_material, role='base'),
                 ], sn=1, sp=1, kind='DA')])

    light_source = LightSource(source_type="standard", version="AM1.5g")

    options = State()
    wl = np.linspace(290, 700, 150) * 1e-9
    options.T = np.random.uniform(10, 350)
    options.wavelength = wl
    options.light_source = light_source
    options.position = None
    options.internal_voltages = np.linspace(-6, 4, 20)
    prepare_solar_cell(test_junc, options)

    solve_beer_lambert(test_junc, options)
    return test_junc, light_source, options
