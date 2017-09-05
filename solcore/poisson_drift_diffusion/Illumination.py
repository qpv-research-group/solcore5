import numpy as np
import os

from solcore import constants
from solcore.absorption_calculator import adachi_alpha
from solcore.light_source import reference_spectra

q = constants.q
pi = constants.pi
Gauss = lambda x, centre, hwhm: 1 / np.sqrt(2 * pi) / (0.5 * hwhm) * np.exp(
    -0.5 * (x - centre) ** 2 / (0.5 * hwhm) ** 2)  # Gaussian (area normalised to 1)


class Illumination(object):
    """ Creates the illumination spectrum in photons m-2 m-1, so it can be used by the drift diffusion module.

    :param spectrum: Either a string with the name of the standard spectrum (AM1.5g, AM1.5d, AM0) or an array containing the power density (in W m-2 nm-1). Default: AM1.5d
    :param irradiance: A multiplicative factor of the spectrum. If the spectrum is in absolutine units, it represents the concentration factor. Default: 1.
    :param wavelengths: An array of the same length that spectrum and containing the corresponding wavelenghts (in m). If a standard spectrum is used, the photon flux is interpolated at these values. Default: None.
    :param central_WL: Central wavelength on an LED type illumination source, in nm.
    :param width: Half-width of the gaussian illumination profile in the LED type source, in nm.
    :param power_density: Power density in the LED type illumination profile, in W/m2.
    :param resolution: Resolution of the spectrum in the LED type illumination profile, in nm.

    """

    def __init__(self, spectrum="AM1.5d", irradiance=1, wavelengths=None, central_WL=532, width=10, power_density=1000,
                 resolution=1):

        if spectrum in ["AM1.5g", "AM1.5d", "AM0"]:
            # DEFAULT. We are using a standard spectrum from the ASTM
            spectrumfile = reference_spectra()

            if spectrum == "AM1.5g":
                self.photon_flux = spectrumfile[:, 2]
                self.total_power = np.trapz(self.photon_flux, spectrumfile[:, 0])
            elif spectrum == "AM0":
                self.photon_flux = spectrumfile[:, 1]
                self.total_power = np.trapz(self.photon_flux, spectrumfile[:, 0])
            elif spectrum == "AM1.5d":
                # DEFAULT
                self.photon_flux = spectrumfile[:, 3]
                self.total_power = np.trapz(self.photon_flux, spectrumfile[:, 0])

            # We interpolate the standard spectrum if we choose a different wavelength range
            if wavelengths is not None:
                self.wl = wavelengths
                self.photon_flux = np.interp(self.wl, spectrumfile[:, 0] * 1e-9, self.photon_flux)
            else:
                self.wl = spectrumfile[:, 0] * 1e-9

        elif spectrum == "LED":
            # We are using an LED type source (LED, laser, etc). We assume a gaussian profile
            # The spread of the spectrum is 5xwidth on each side
            central_WL = central_WL * 1e-9
            width = width * 1e-9
            resolution = resolution * 1e-9

            self.wl = np.arange(central_WL - 5 * width, central_WL + 5 * width, resolution)
            self.photon_flux = power_density * Gauss(self.wl, central_WL, width) * 1e-9
            self.total_power = np.trapz(self.photon_flux, self.wl) * 1e9

        else:
            # We are using a custom spectrum
            self.wl = wavelengths
            self.photon_flux = spectrum
            self.total_power = np.trapz(self.photon_flux, self.wl)

        # The above is not the actual photon flux but power density. We calculate the photon flux
        self.power_density = self.photon_flux
        self.photon_flux = self.photon_flux * self.wl / 1240 / q * 1e18

        # We calculate some extra things    
        self.OD = np.zeros(len(self.photon_flux))
        self.filter_transmission = np.exp(-self.OD)

        if irradiance == 0:
            self.irradiance = 1e-9
        else:
            self.irradiance = irradiance

        output = self.photon_flux * self.filter_transmission * self.irradiance
        self.total_current = np.trapz(output * q, self.wl)

    def filter(self, edge=None, OD=2, material=None, thickness=0, T=293, func=None):
        """Creates a filter to be applied to the spectrum.

        It can be either a square, long pass filter, defined by an edge wavelength (edge) and an optical density (OD)
        or assuming a layer of certain thickness (thickness) made of certain solcore material (material). It can also
        be provided externally as a function (func) that calculates a wavelength dependent optical density.

        The optical density is used to calculate the transmission of the filter as:

        .. math:: T = exp(-OD)

        :param edge: Wavelength edge
        :param OD: Optcial density above the edge
        :param material: Solcore material
        :param thickness: Thickness of the solcore material
        :param T: Temperature
        :param func: Function that takes wavelengths as input and returns the optical density at those wavelengths
        :return: None
        """
        if func is not None:
            # We are using a function that provides the optical density
            assert callable(func)
            self.OD = func(self.wl)
        elif edge is not None:
            # DEFAULT. We are using a square filter with some "nice", 10nm wide, exponential drop at the edge.
            self.OD[self.wl <= edge] = OD
            self.OD[self.wl > edge] = OD * np.exp(- (self.wl[self.wl > edge] - edge) ** 2 / 2 / 10e-9 ** 2)
        elif material is not None:
            # We are using a layer made of a certain material and thickness as filter
            try:
                self.OD = material.alpha(self.wl) * thickness
            except:
                self.OD = adachi_alpha.create_adachi_alpha(material, T=T, wl=self.wl) * thickness
        else:
            print("Error: filter settings are empty. No filter will be used ")

        # We transform the optical density in actual absorptivity    
        self.filter_transmission = np.exp(-self.OD)

    def spectrum(self, irradiance=None):
        """ Returns the filtered spectrum optionally scaled by an irradiance.

        :param irradiance: Irradiance multiplicative factor (in Suns).
        :return: spectrum weighted by the irradiance and any filter defined previously.
        """

        if irradiance is not None:
            if irradiance == 0:
                self.irradiance = 1e-9
            else:
                self.irradiance = irradiance

        output = self.photon_flux * self.filter_transmission * self.irradiance
        self.total_current = np.trapz(output * q, self.wl)
        self.total_power = np.trapz(self.power_density, self.wl) * 1e9 * self.irradiance

        return output
