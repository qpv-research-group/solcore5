import numpy as np
from typing import Callable, Union, Optional
from functools import wraps
from abc import ABC, abstractmethod
from copy import deepcopy
from pytest import approx

from solcore import (
    spectral_conversion_nm_ev,
    spectral_conversion_nm_hz,
    eVnm,
    nmHz,
    nmJ,
)
from solcore.constants import q, h, c


REGISTERED_CONVERTERS: dict = {}
""" Registered spectrum conversion functions."""


SPECTRUM_SIGNATURE = Callable[[Union[np.ndarray, float]], Union[np.ndarray, float]]
"""Signature for the spectrum functions."""


class LightSourceBase(ABC):
    def __init__(
        self, output_units="power_density_per_nm", concentration=1, **kwargs,
    ):
        """ Base class for all light sources

        :param x: Array with the spectral range in which to calculate the
        spectrum. It must be in the "units" defined by the output_units parameter.
        :param output_units: Units of the output spectrum.
        :param concentration: Concentration of the light source.
        :param kwargs: Options to update the light source.
        """
        msg = (
            f"Unknown output units {output_units}. "
            f"Valid options are: {tuple(REGISTERED_CONVERTERS.keys())}"
        )
        assert output_units in REGISTERED_CONVERTERS, msg

        self.power_density = 0
        self.output_units = output_units
        self.concentration = concentration
        self.options = deepcopy(kwargs)
        self.cache = None

    def spectrum(
        self,
        x: Union[np.ndarray, float],
        output_units: str = None,
        concentration: Union[float, int] = None,
        **kwargs,
    ):
        """ Returns the spectrum of the light in the requested units. Internally,
        the spectrum is always managed in power density per nanometers, but the
        output can be provided in other compatible units, such as power density per
        Hertz or photon flux per eV.

        :param x: (Default=None) If "x" is provided, it must be an array with the
        spectral range in which to calculate the spectrum. It must be in the "units"
        defined when creating the light source, or in those required by the
        output_units parameter, if provided.
        :param output_units: Units of the output spectrum.
        :param concentration: Concentration of the light source.
        :param kwargs: Options to update the light source. It can be any of the
        options specific to the chosen type of light already defined source.
        :return: Array with the spectrum in the requested units.
        """
        if (
            self.cache
            and not output_units
            and not concentration
            and not kwargs
            and x == approx(self.cache[0])
        ):
            return self.cache

        output_units = output_units if output_units is not None else self.output_units
        con = concentration if concentration is not None else self.concentration

        msg = f"Valid units are: {list(REGISTERED_CONVERTERS.keys())}."
        assert output_units in REGISTERED_CONVERTERS, msg

        self.options.update({k: v for k, v in kwargs.items() if k in self.options})
        self.cache = (x, REGISTERED_CONVERTERS[output_units](self._spectrum, x) * con)
        return self.cache

    @abstractmethod
    def _spectrum(self) -> SPECTRUM_SIGNATURE:
        """Provides the callable that calculates the spectrum at the given inputs.

        The callable will be passed to the conversion function corresponding
        to the requested output units and, therefore, it MUST accept one single input
        parameter, the wavelengths in nm, and provide as output the spectrum in
        power_density_per_nm.
        """
        pass


def register_conversion_function(fun: Callable):
    """Registers a view method that will trigger an event. """

    @wraps(fun)
    def wrapper(*args, **kwargs):
        return fun(*args, **kwargs)

    REGISTERED_CONVERTERS[fun.__name__] = wrapper

    return wrapper


@register_conversion_function
def power_density_per_nm(spectrum: Callable[[np.ndarray], np.ndarray], x: np.ndarray):
    """ Function that returns the spectrum in power density per nanometer.

    The input spectrum is assumed to be in power density per nanometer.

    :param spectrum: The spectrum to interpolate.
    :param x: Array with the wavelengths (in nm)
    :return: The spectrum in the chosen units.
    """
    return spectrum(x)


@register_conversion_function
def photon_flux_per_nm(spectrum: Callable[[np.ndarray], np.ndarray], x: np.ndarray):
    """ Function that returns the spectrum in photon flux per nanometer.

    The input spectrum is assumed to be in power density per nanometer.

    :param spectrum: The spectrum to interpolate.
    :param x: Array with the wavelengths (in nm)
    :return: The spectrum in the chosen units.
    """
    return spectrum(x) / (c * h * 1e9 / x)


@register_conversion_function
def power_density_per_m(spectrum: Callable[[np.ndarray], np.ndarray], x: np.ndarray):
    """ Function that returns the spectrum in power density per meter.

    The input spectrum is assumed to be in power density per nanometer.

    :param spectrum: The spectrum to interpolate.
    :param x: Array with the wavelengths (in m)
    :return: The spectrum in the chosen units.
    """
    return spectrum(x * 1e9) * 1e9


@register_conversion_function
def photon_flux_per_m(spectrum: Callable[[np.ndarray], np.ndarray], x: np.ndarray):
    """ Function that returns the spectrum in photon flux per meter.

    The input spectrum is assumed to be in power density per nanometer.

    :param spectrum: The spectrum to interpolate.
    :param x: Array with the wavelengths (in m)
    :return: The spectrum in the chosen units.
    """
    return spectrum(x * 1e9) / (c * h / x) * 1e9


@register_conversion_function
def power_density_per_ev(spectrum: Callable[[np.ndarray], np.ndarray], x: np.ndarray):
    """ Function that returns the spectrum in power density per eV.

    The input spectrum is assumed to be in power density per nanometer.

    :param spectrum: The spectrum to interpolate.
    :param x: Array with the energies (in eV)
    :return: The spectrum in the chosen units.
    """
    wavelength = eVnm(x)[::-1]
    output = spectrum(wavelength)
    _, output = spectral_conversion_nm_ev(wavelength, output)
    return output


@register_conversion_function
def photon_flux_per_ev(spectrum: Callable[[np.ndarray], np.ndarray], x: np.ndarray):
    """ Function that returns the spectrum in photon flux per eV.

    The input spectrum is assumed to be in power density per nanometer.

    :param spectrum: The spectrum to interpolate.
    :param x: Array with the energies (in eV)
    :return: The spectrum in the chosen units.
    """
    wavelength = eVnm(x)[::-1]
    output = spectrum(wavelength)
    _, output = spectral_conversion_nm_ev(wavelength, output)
    return output / (q * x)


@register_conversion_function
def power_density_per_joule(
    spectrum: Callable[[np.ndarray], np.ndarray], x: np.ndarray
):
    """ Function that returns the spectrum in power density per Joule.

    The input spectrum is assumed to be in power density per nanometer.

    :param spectrum: The spectrum to interpolate.
    :param x: Array with the energies (in J)
    :return: The spectrum in the chosen units.
    """
    wavelength = nmJ(x)[::-1]
    output = spectrum(wavelength)
    _, output = spectral_conversion_nm_ev(wavelength, output)
    return output / q


@register_conversion_function
def photon_flux_per_joule(spectrum: Callable[[np.ndarray], np.ndarray], x: np.ndarray):
    """ Function that returns the spectrum in photon flux per Joule.

    The input spectrum is assumed to be in power density per nanometer.

    :param spectrum: The spectrum to interpolate.
    :param x: Array with the energies (in J)
    :return: The spectrum in the chosen units.
    """
    wavelength = nmJ(x)[::-1]
    output = spectrum(wavelength)
    _, output = spectral_conversion_nm_ev(wavelength, output)
    return output / (q * x)


@register_conversion_function
def power_density_per_hz(spectrum: Callable[[np.ndarray], np.ndarray], x: np.ndarray):
    """ Function that returns the spectrum in power density per hertz.

    The input spectrum is assumed to be in power density per nanometer.

    :param spectrum: The spectrum to interpolate.
    :param x: Array with the frequencies (in hz)
    :return: The spectrum in the chosen units.
    """
    wavelength = nmHz(x)[::-1]
    output = spectrum(wavelength)
    _, output = spectral_conversion_nm_hz(wavelength, output)
    return output


@register_conversion_function
def photon_flux_per_hz(spectrum: Callable[[np.ndarray], np.ndarray], x: np.ndarray):
    """ Function that returns the spectrum in photon flux per hertz.

    The input spectrum is assumed to be in power density per nanometer.

    :param spectrum: The spectrum to interpolate.
    :param x: Array with the frequencies (in hz)
    :return: The spectrum in the chosen units.
    """
    wavelength = nmHz(x)[::-1]
    output = spectrum(wavelength)
    frequency, output = spectral_conversion_nm_hz(wavelength, output)
    return output / (h * frequency)
