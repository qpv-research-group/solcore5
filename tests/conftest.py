from pytest import fixture
import numpy as np
from unittest.mock import patch
import os
import tempfile

temp_folder = tempfile.TemporaryDirectory()
SOLCORE_USER_DATA = temp_folder.name
os.environ["SOLCORE_USER_DATA"] = SOLCORE_USER_DATA


def patch_plots(function):
    from functools import wraps

    @wraps(function)
    def decorated(*args, **kwargs):

        with patch("matplotlib.pyplot.show", lambda *x, **y: None):
            import matplotlib
            matplotlib.use("Agg")
            return function(*args, **kwargs)

    return decorated


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


def da_light_source():
    from solcore.light_source import LightSource

    light_source = LightSource(source_type="standard", version="AM1.5g")
    return light_source


def da_options():
    from solcore.state import State

    options = State()
    wl = np.linspace(290, 700, 150) * 1e-9
    options.T = np.random.uniform(250, 350)
    options.wavelength = wl
    options.light_source = da_light_source()
    options.position = None
    options.internal_voltages = np.linspace(-6, 4, 20)

    return options


def junction(nd_top, na_top, nd_bottom, na_bottom):
    from solcore.structure import Junction, Layer
    from solcore import si, material
    from solcore.solar_cell import SolarCell
    from solcore.constants import vacuum_permittivity
    from solcore.solar_cell_solver import prepare_solar_cell
    from solcore.optics import solve_beer_lambert

    Lp = np.power(10, np.random.uniform(-8, -6))
    Ln = np.power(10, np.random.uniform(-8, -6))

    AlInP = material("AlInP")
    InGaP = material("GaInP")
    window_material = AlInP(Al=0.52)
    top_cell_n_material = InGaP(
        In=0.48,
        Nd=nd_top,
        Na=na_top,
        hole_diffusion_length=Lp,
        electron_diffusion_length=Ln,
    )
    top_cell_p_material = InGaP(
        In=0.48,
        Nd=nd_bottom,
        Na=na_bottom,
        hole_diffusion_length=Lp,
        electron_diffusion_length=Ln,
    )

    rel_perm = np.random.uniform(1, 20)
    for mat in [top_cell_n_material, top_cell_p_material]:
        mat.permittivity = rel_perm * vacuum_permittivity

    n_width = np.random.uniform(500, 1000) * 1e-9
    p_width = np.random.uniform(3000, 5000) * 1e-9

    test_junc = SolarCell(
        [
            Junction(
                [
                    Layer(si("25nm"), material=window_material, role="window"),
                    Layer(n_width, material=top_cell_n_material, role="emitter"),
                    Layer(p_width, material=top_cell_p_material, role="base"),
                ],
                sn=1,
                sp=1,
                kind="DA",
            )
        ]
    )

    options = da_options()
    options.light_source = da_light_source()
    prepare_solar_cell(test_junc, options)
    solve_beer_lambert(test_junc, options)

    return test_junc, options


@fixture
def np_junction():

    Nd_top = np.power(10, np.random.uniform(24, 26))
    Na_bottom = np.power(10, np.random.uniform(23, 25))

    Na_top = 0
    Nd_bottom = 0

    test_junc, options = junction(Nd_top, Na_top, Nd_bottom, Na_bottom)
    return test_junc, options


@fixture
def pn_junction():

    Na_top = np.power(10, np.random.uniform(24, 26))
    Nd_bottom = np.power(10, np.random.uniform(23, 25))

    Na_bottom = 0
    Nd_top = 0

    test_junc, options = junction(Nd_top, Na_top, Nd_bottom, Na_bottom)
    return test_junc, options


@fixture
def prepare_test_cell():

    from solcore import material, si
    from solcore.structure import Junction, Layer, TunnelJunction
    from solcore.solar_cell import SolarCell

    GaAs = material("GaAs")()
    MgF2 = material("MgF2")()
    TiO2 = material("TiO2")()
    Ge = material("Ge")()

    widths = np.random.rand(9) * 200

    solar_cell = SolarCell(
        [
            Layer(si(widths[0], "nm"), material=MgF2),
            Layer(si(widths[1], "nm"), material=TiO2),
            Junction(
                [
                    Layer(si(widths[2], "nm"), material=GaAs, role="window"),
                    Layer(si(widths[3], "nm"), material=GaAs, role="emitter"),
                    Layer(si(widths[4], "nm"), material=GaAs, role="base"),
                ],
                kind="DA",
            ),
            TunnelJunction(
                [
                    Layer(si(widths[5], "nm"), material=GaAs),
                    Layer(si(widths[6], "nm"), material=GaAs),
                ]
            ),
            Junction(
                [
                    Layer(si(widths[7], "nm"), material=Ge, role="emitter"),
                    Layer(si(widths[8], "nm"), material=Ge, role="base"),
                ],
                kind="PDD",
            ),
        ]
    )

    return solar_cell, widths


@fixture
def AlGaAs():
    from solcore import material
    from solcore.structure import Layer, Junction
    from solcore import si

    T = 298
    # We create the other materials we need for the device
    window = material("AlGaAs")(T=T, Na=5e24, Al=0.8)
    p_AlGaAs = material("AlGaAs")(T=T, Na=1e24, Al=0.4)
    n_AlGaAs = material("AlGaAs")(T=T, Nd=8e22, Al=0.4)
    bsf = material("AlGaAs")(T=T, Nd=2e24, Al=0.6)

    return Junction(
        [
            Layer(width=si("30nm"), material=window, role="Window"),
            Layer(width=si("150nm"), material=p_AlGaAs, role="Emitter"),
            Layer(width=si("1000nm"), material=n_AlGaAs, role="Base"),
            Layer(width=si("200nm"), material=bsf, role="BSF"),
        ],
        sn=1e6,
        sp=1e6,
        T=T,
        kind="PDD",
    )


@fixture
def light_source(wavelength):
    from solcore.light_source import LightSource

    return LightSource(
        source_type="standard",
        version="AM1.5g",
        x=np.linspace(350, 1000, 301) * 1e-9,
        output_units="photon_flux_per_m",
        concentration=1,
    )
