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
