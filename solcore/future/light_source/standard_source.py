import os
import numpy as np
from scipy import interpolate
from solcore.science_tracker import science_reference

from .light_source_base import LightSourceBase, SPECTRUM_SIGNATURE


def reference_spectra():
    """Function providing the standard reference spectra: AM0, AM1.5g and AM1.5d.

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


class StandardSource(LightSourceBase):
    @classmethod
    def am0(cls, output_units="power_density_per_nm", concentration=1):
        spectra = reference_spectra()
        return cls(
            version="AM0",
            x=spectra[:, 0],
            y=spectra[:, 1],
            output_units=output_units,
            concentration=concentration,
        )

    @classmethod
    def am15g(cls, output_units="power_density_per_nm", concentration=1):
        spectra = reference_spectra()
        return cls(
            version="AM1.5g",
            x=spectra[:, 0],
            y=spectra[:, 2],
            output_units=output_units,
            concentration=concentration,
        )

    @classmethod
    def am15d(cls, output_units="power_density_per_nm", concentration=1):
        spectra = reference_spectra()
        return cls(
            version="AM1.5d",
            x=spectra[:, 0],
            y=spectra[:, 3],
            output_units=output_units,
            concentration=concentration,
        )

    def __init__(
        self,
        version: str,
        x: np.ndarray,
        y: np.ndarray,
        output_units="power_density_per_nm",
        concentration=1,
    ):
        super().__init__(output_units, concentration)
        self.version = version
        self._x = x
        self._y = y
        self.power_density = np.trapz(y=self._y, x=self._x) * self.concentration

    def _spectrum(self) -> SPECTRUM_SIGNATURE:
        return interpolate.interp1d(
            x=self._x, y=self._y, bounds_error=False, fill_value=0, assume_sorted=True
        )
