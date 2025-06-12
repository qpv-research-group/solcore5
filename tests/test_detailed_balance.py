import pytest
import numpy as np
from solcore.structure import Junction
from pytest import approx
from solcore.state import State
from solcore.analytic_solar_cells.detailed_balance import (
    iv_detailed_balance,
    qe_detailed_balance,
    absorptance_detailed_balance,
    surface_integral,
)


def test_absorptance_detailed_balance():
    # Test the quantum efficiency calculation
    junc = Junction(Eg=1.3, A=0.9)

    absorptance_detailed_balance(junc)

    assert junc.absorptance((1240 / junc.Eg - 10) * 1e-9) == junc.A
    assert junc.absorptance((1240 / junc.Eg + 10) * 1e-9) == 0


def test_absorptance_detailed_balance_error():
    # Test the quantum efficiency calculation
    junc = Junction(Eg=1.3)

    with pytest.raises(AttributeError):

        absorptance_detailed_balance(junc)


def test_surface_integral_perfectA():

    junc = Junction(Eg=1240 / 900, reflected=lambda x: 0 * np.ones_like(x), A=1)

    absorptance_detailed_balance(junc)
    wavelength = np.linspace(300, 1000, 50) * 1e-9
    As = surface_integral(junc, wavelength)

    assert np.all(As[wavelength * 1e9 > 1240 / junc.Eg] == 0)
    assert np.all(As[wavelength * 1e9 < 1240 / junc.Eg] == 2 * np.pi)


def test_surface_integral_nonperfectA():

    junc = Junction(Eg=1240 / 900, reflected=lambda x: 0.05 * np.ones_like(x), A=0.8)

    absorptance_detailed_balance(junc)
    wavelength = np.linspace(300, 1000, 50) * 1e-9
    As = surface_integral(junc, wavelength)

    assert np.all(As[wavelength * 1e9 > 1240 / junc.Eg] == 0)
    assert np.all(As[wavelength * 1e9 < 1240 / junc.Eg] < 2 * np.pi)


def test_qe_detailed_balance():
    # Test the quantum efficiency calculation
    junc = Junction(Eg=1.3, A=0.9)
    wavelength = np.linspace(300, 1000, 50) * 1e-9

    qe_detailed_balance(junc, wavelength)

    assert np.all(junc.qe["EQE"][wavelength * 1e9 > 1240 / junc.Eg] == 0)
    assert np.all(junc.qe["IQE"][wavelength * 1e9 > 1240 / junc.Eg] == 0)
    assert np.all(junc.qe["EQE"][wavelength * 1e9 < 1240 / junc.Eg] == junc.A)
    assert np.all(junc.qe["IQE"][wavelength * 1e9 < 1240 / junc.Eg] == junc.A)


def test_qe_detailed_balance_absorbed():
    # Test the quantum efficiency calculation with an absorption profile from an optical calculation
    wavelength = np.linspace(300, 1000, 50) * 1e-9
    alpha = 1e5 * np.ones_like(wavelength)

    alpha[wavelength * 1e9 > 1240 / 1.3] = (
        0  # Absorption coefficient is zero above the bandgap
    )

    def absorption_profile(z):
        return alpha[None, :] * np.exp(-alpha[None, :] * z[:, None])

    junc = Junction(Eg=1.3, A=1, width=10e-6, absorbed=absorption_profile)
    wavelength = np.linspace(300, 1000, 50) * 1e-9

    qe_detailed_balance(junc, wavelength)

    B = 1 - np.exp(-alpha[0] * junc.width)

    assert np.all(junc.qe["EQE"][wavelength * 1e9 > 1240 / junc.Eg] == 0)
    assert np.all(junc.qe["IQE"][wavelength * 1e9 > 1240 / junc.Eg] == 0)
    assert junc.qe["EQE"][wavelength * 1e9 < 1240 / junc.Eg] == approx(B)
    assert np.all(junc.qe["IQE"][wavelength * 1e9 < 1240 / junc.Eg] == junc.A)


def test_iv_detailed_balance_dark():
    options = State()
    options.T = 300
    options.light_iv = False
    options.db_mode = "boltzmann"
    options.T_ambient = 300
    options.wavelength = np.linspace(300, 1000, 50) * 1e-9
    options.internal_voltages = np.linspace(0, 1.2, 20)

    junc = Junction(Eg=1.3, A=0.9)

    with pytest.raises(AttributeError):
        iv_detailed_balance(junc, options)

    junc.n = 3  # refractive index
    junc.reflected = lambda x: 0.05 * np.ones_like(x)  # some reflection
    iv_detailed_balance(junc, options)

    iv_nobackref = junc.iv(options.internal_voltages)

    junc.back_reflector = True

    iv_detailed_balance(junc, options)
    iv_backref = junc.iv(options.internal_voltages)

    assert iv_nobackref / 2 == approx(iv_backref, rel=0.01)

    options.db_mode = "planck"
    iv_detailed_balance(junc, options)

    iv_planck = junc.iv(options.internal_voltages)

    assert iv_planck == approx(iv_backref, rel=0.01)





def test_iv_detailed_balance_light():
    from solcore.light_source import LightSource
    options = State()
    options.T = 300
    options.light_iv = True
    options.db_mode = "boltzmann"
    options.T_ambient = 300
    options.wavelength = np.linspace(300, 1000, 50) * 1e-9
    options.internal_voltages = np.linspace(0, 1.3, 20)
    options.light_source = LightSource(source_type='standard', x=options.wavelength,
                                       version='AM1.5g',
                           output_units="photon_flux_per_m")

    junc = Junction(Eg=1.3, A=0.9)

    junc.n = 3  # refractive index
    junc.reflected = lambda x: 0.05 * np.ones_like(x)  # some reflection
    iv_detailed_balance(junc, options)

    iv_nobackref = junc.iv(options.internal_voltages)

    junc.back_reflector = True

    iv_detailed_balance(junc, options)
    iv_backref = junc.iv(options.internal_voltages)

    assert iv_nobackref[0] == approx(iv_backref[0], rel=0.001)

    assert np.all(iv_nobackref >= iv_nobackref)
