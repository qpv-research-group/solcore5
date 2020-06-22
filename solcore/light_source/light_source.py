import os
from scipy.interpolate import interp1d
import numpy as np
from typing import Callable, Optional
from functools import wraps

from solcore.science_tracker import science_reference
from solcore import (
    spectral_conversion_nm_ev,
    spectral_conversion_nm_hz,
    eVnm,
    nmHz,
    nmJ,
)
from solcore.constants import q, h, c, kb

from solcore.light_source.spectral2 import (
    get_default_spectral2_object,
    calculate_spectrum_spectral2,
)
from solcore.light_source.smarts import (
    get_default_smarts_object,
    calculate_spectrum_smarts,
)


REGISTERED_CONVERTERS: dict = {}
""" Registered spectrum conversion functions."""


def reference_spectra():
    """ Function providing the standard reference spectra: AM0, AM1.5g and AM1.5d.

    :return: A 2D array with 4 columns representing the wavelength, AM0, AM1.5g and
    AM1.5d standard spectra.
    """

    science_reference(
        "Standard solar spectra",
        "ASTM G173-03(2012), Standard Tables for Reference Solar Spectral Irradiances: "
        "Direct Normal and Hemispherical on 37° Tilted Surface, ASTM International, "
        "West Conshohocken, PA, 2012, www.astm.org",
    )

    this_dir = os.path.split(__file__)[0]
    output = np.loadtxt(
        os.path.join(this_dir, "astmg173.csv"), dtype=float, delimiter=",", skiprows=2
    )

    return output


class LightSource:
    """ Common interface to access all types of light sources supported by Solcore.

    It includes standard solar spectra (AM0, AM1.5g, and AM1.5d), blackbody radiation,
    laser light or spectra created from atmospheric data using SPECTRAL2 or SMARTS.
    Additionally, it can also use experimentally measured spectra.
    """

    type_of_source = [
        "laser",
        "black body",
        "standard",
        "SMARTS",
        "SPECTRAL2",
        "custom",
    ]

    def __init__(self, source_type, x=None, output_units = "power_density_per_nm", concentration=1, **kwargs):
        """

        :param source_type:
        :param kwargs:
        """
        msg = f"Unknown source {source_type}. " \
              f"Valid options are: {self.type_of_source}"
        assert source_type in self.type_of_source, msg

        msg = f"Unknown output units {output_units}. " \
              f"Valid options are: {tuple(REGISTERED_CONVERTERS.keys())}"
        assert output_units in REGISTERED_CONVERTERS, msg

        self.source_type = source_type
        self.x = x
        self.x_internal = x
        self.power_density = 0

        self.options = {}
        self.output_units = output_units
        self.concentration = concentration
        self.options.update(kwargs)
        self._spectrum = None

        self._update_get_spectrum(self.output_units)
        self._update_spectrum_function()

        self.ready = False
        self.cache_spectrum = None

    def spectrum(self, x=None, output_units=None, concentration=None, **kwargs):
        """ Returns the spectrum of the light in the requested units. Internally,
        the spectrum is always managed in power density per nanometers, but the
        output can be provided in other compatible units, such as power density per
        Hertz or photon flux per eV.

        :param x: (Default=None) If "x" is provided, it must be an array with the
        spectral range in which to calculate the spectrum. Depending on the "units"
        defined when creating the light source, this array must be in nm, m, eV,
        J or hz.
        :param output_units: Units of the output spectrum
        :param concentration: Concentration of the light source
        :param kwargs: Options to update the light source. It can be "units",
        "concentration" or any of the options specific to the chosen type of light
        source.
        :return: Array with the spectrum in the requested units
        """
        self.x = x if x is not None else self.x_internal

        output_units = output_units if output_units is not None else self.output_units
        con = concentration if concentration is not None else self.concentration

        self._update_get_spectrum(output_units)
        if kwargs:
            self._update(**kwargs)
            self._update_spectrum_function()
        return self.x, self._get_spectrum(self._spectrum, self.x) * con

    def _update(self, **kwargs):
        """ Updates the options of the light source with new values.

        It only updates existing options. No new options are added.

        :param kwargs: A dictionary with the options to update and their new values.
        :return: None
        """

        for opt in kwargs:
            if opt in self.options:
                self.options[opt] = kwargs[opt]

        self.ready = False
        self.cache_spectrum = None

    def _update_get_spectrum(self, output_units):
        """ Updates the function to get the spectrum, depending on the chosen output units.

        :return: None
        """
        msg = f"Valid units are: {list(REGISTERED_CONVERTERS.keys())}."
        assert output_units in REGISTERED_CONVERTERS, msg

        self._get_spectrum = REGISTERED_CONVERTERS[output_units]

    def _update_spectrum_function(self):
        """ Updates the spectrum function during the light source creation or just after updating one or more of the options. It also updates the "options" property with any default options available to the chosen light source, if any.

        :return: True
        """
        try:
            if self.source_type == "standard":
                self._spectrum = self._get_standard_spectrum(self.options)
            elif self.source_type == "laser":
                self._spectrum = self._get_laser_spectrum(self.options)
            elif self.source_type == "black body":
                self._spectrum = self._get_blackbody_spectrum(self.options)
            elif self.source_type == "SPECTRAL2":
                self._spectrum = self._get_spectral2_spectrum(self.options)
            elif self.source_type == "SMARTS":
                self._spectrum = self._get_smarts_spectrum(self.options)
            elif self.source_type == "custom":
                self._spectrum = self._get_custom_spectrum(self.options)
            else:
                raise ValueError(
                    "Unknown type of light source: {0}.\nValid light sources are: {1}".format(
                        self.source_type, self.type_of_source
                    )
                )

            self.ready = True

        except ValueError as err:
            print(err)

    def _get_standard_spectrum(self, options):
        """ Gets one of the reference standard spectra: AM0, AM1.5g or AM1.5d.

        :param options: A dictionary that contains the 'version' of the standard spectrum: 'AM0', 'AM1.5g' or 'AM1.5d'
        :return: A function that takes as input the wavelengths and return the standard spectrum at those wavelengths.
        """

        try:
            version = options["version"]

            spectra = reference_spectra()
            wl = spectra[:, 0]

            if version == "AM0":
                spectrum = spectra[:, 1]
            elif version == "AM1.5g":
                spectrum = spectra[:, 2]
            elif version == "AM1.5d":
                spectrum = spectra[:, 3]
            else:
                raise KeyError(
                    'ERROR when creating a standard light source. Input parameters must include "version" '
                    'which can be equal to "AM0", "AM1.5g" or "AM1.5d" only.'
                )

            self.x_internal = wl
            self.power_density = (
                np.trapz(y=spectrum, x=wl) * self.concentration
            )
            output = interp1d(
                x=wl, y=spectrum, bounds_error=False, fill_value=0, assume_sorted=True
            )
            return output

        except KeyError as err:
            print(err)

    def _get_laser_spectrum(self, options):
        """ Creates a gaussian light source with a given total power, linewidth and central wavelength. These three parameters must be provided in the "options" diccionary.

        :param options: A dictionary that must contain the 'power', the 'linewidth' and the 'center' of the laser emission.
        :return: A function that takes as input the wavelengths and return the laser spectrum at those wavelengths.
        """
        try:
            power = options["power"]
            sigma2 = options["linewidth"] ** 2
            center = options["center"]

            def output(x):
                out = (
                    power
                    / np.sqrt(2 * np.pi * sigma2)
                    * np.exp(-(x - center) ** 2 / 2 / sigma2)
                )
                return out

            self.x_internal = np.arange(
                center - 5 * options["linewidth"],
                center + 5 * options["linewidth"],
                options["linewidth"] / 20,
            )
            self.power_density = power * self.concentration
            return output

        except KeyError:
            print(
                'ERROR when creating a laser light source. Input parameters must include "power", "linewidth"'
                ' and "center".'
            )

    def _get_blackbody_spectrum(self, options):
        """ Gets the expontaneous emission in W/m2/sr/nm from a black body source chemical potential = 0

        :param options: A dictionary that must contain the temperature of the blackbody, in kelvin 'T' and the 'entendue' in sr. If not provided, the entendue will be taken as 1 sr. Possible values for entendue are:
            - 'Sun': The entendue will be taken as 6.8e-5 sr.
            - 'hemispheric': The entendue wil be taken as pi/2 sr.
            - A numeric value
        :return: A function that takes as input the wavelengths and return the black body spectrum at those wavelengths.
        """

        try:
            T = options["T"]
            if "entendue" in options:
                if options["entendue"] == "Sun":
                    entendue = 6.8e-5
                elif options["entendue"] == "hemispheric":
                    entendue = np.pi / 2
                else:
                    entendue = options["entendue"]
            else:
                entendue = 1
                options["entendue"] = 1

            def BB(x):
                x = x * 1e-9
                out = (
                    2
                    * entendue
                    * h
                    * c ** 2
                    / x ** 5
                    / (np.exp(h * c / (x * kb * T)) - 1)
                )
                return out * 1e-9

            wl_max = 2.897_772_9e6 / T
            self.x_internal = np.arange(0, wl_max * 10, wl_max / 100)
            sigma = 5.670_367e-8
            self.power_density = (
                sigma * T ** 4 * entendue / np.pi * self.concentration
            )

            return BB

        except KeyError:
            print(
                'ERROR when creating a blackbody light source. Input parameters must include "T" and, '
                'optionally, an "entendue", whose values can be "Sun" "hemispheric" or a number. '
                "Equal to 1 if omitted."
            )

    def _get_spectral2_spectrum(self, options):
        """ Get the solar spectrum calculated with the SPECTRAL2 calculator. The options dictionary is updated with the default options of all parameters if they are not provided.

        :param options: A dictionary that contain all the options for the calculator.
        :return: A function that takes as input the wavelengths and return the SPECTRAL2 calculated spectrum at those wavelengths.
        """
        default = get_default_spectral2_object()

        for opt in options:
            if opt in default:
                default[opt] = options[opt]

        options.update(default)

        wl, irradiance = calculate_spectrum_spectral2(options, power_density_in_nm=True)

        self.x_internal = wl
        self.power_density = (
            np.trapz(y=irradiance, x=wl) * self.concentration
        )
        output = interp1d(
            x=wl, y=irradiance, bounds_error=False, fill_value=0, assume_sorted=True
        )
        return output

    def _get_smarts_spectrum(self, options):
        """ Get the solar spectrum calculated with the SMARTS calculator. The options dictionary is updated with the default options of all parameters if they are not provided.

        :param options: A dictionary that contain all the options for the calculator.
        :return: A function that takes as input the wavelengths and return the SMARTS calculated spectrum at those wavelengths.
        """
        outputs = {
            "Extraterrestial": 2,
            "True direct": 2,
            "Experimental direct": 3,
            "Global horizontal": 4,
            "Global tilted": 5,
        }

        default = get_default_smarts_object()

        for opt in options:
            if opt in default:
                default[opt] = options[opt]

        options.update(default)

        if "output" not in options:
            options["output"] = "Global horizontal"

        try:
            output = outputs[options["output"]]

            out = calculate_spectrum_smarts(options)

            self.x_internal = out[0]
            self.power_density = (
                np.trapz(y=out[output], x=out[0]) * self.concentration
            )
            output = interp1d(
                x=out[0],
                y=out[output],
                bounds_error=False,
                fill_value=0,
                assume_sorted=True,
            )
            return output

        except KeyError:
            print(
                "ERROR: Output option no recognized. Avaliable options are: {}".format(
                    outputs
                )
            )
        except RuntimeError as err:
            print("ERROR in SMARTS: {}".format(err))

    def _get_custom_spectrum(self, options):
        """ Convert an custom spectrum into a solcore light source object.

        :param options: A dictionary that contains the following information:
            - 'x_data' and 'y_data' of the custom spectrum.
            - 'input_units', the units of the spectrum, such as 'photon_flux_per_nm' or 'power_density_per_eV'
        :return: A function that takes as input the wavelengths and return the custom spectrum at those wavelengths.
        """

        try:
            x_data = options["x_data"]
            y_data = options["y_data"]
            units = options["input_units"]

            # We check the units type by type.
            # Regardless of the input, we want power_density_per_nm
            if units == "power_density_per_nm":
                wl = x_data
                spectrum = y_data
            elif units == "photon_flux_per_nm":
                wl = x_data
                spectrum = y_data * (c * h * 1e9 / wl)
            elif units == "power_density_per_m":
                wl = x_data * 1e-9
                spectrum = y_data * 1e-9
            elif units == "photon_flux_per_m":
                wl = x_data * 1e-9
                spectrum = y_data * (c * h / wl)
            elif units == "power_density_per_eV":
                wl, spectrum = spectral_conversion_nm_ev(x_data, y_data)
            elif units == "photon_flux_per_eV":
                wl, spectrum = spectral_conversion_nm_ev(x_data, y_data)
                spectrum = spectrum * (c * h * 1e9 / wl)
            elif units == "power_density_per_J":
                wl, spectrum = spectral_conversion_nm_ev(x_data / q, y_data * q)
            elif units == "photon_flux_per_J":
                wl, spectrum = spectral_conversion_nm_ev(x_data / q, y_data * q)
                spectrum = spectrum * (c * h * 1e9 / wl)
            elif units == "power_density_per_hz":
                wl, spectrum = spectral_conversion_nm_hz(x_data, y_data)
            elif units == "photon_flux_per_hz":
                wl, spectrum = spectral_conversion_nm_hz(x_data, y_data)
                spectrum = spectrum * (h * x_data)
            else:
                raise ValueError(
                    "Unknown units: {0}.\nValid units are: {1}.".format(
                        units, list(REGISTERED_CONVERTERS.keys())
                    )
                )

            self.x_internal = wl
            self.power_density = (
                np.trapz(y=spectrum, x=wl) * self.concentration
            )
            output = interp1d(
                x=wl, y=spectrum, bounds_error=False, fill_value=0, assume_sorted=True
            )
            return output

        except KeyError as err:
            print(err)


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
