from typing import Optional
from abc import ABC, abstractmethod

import numpy as np
import xarray as xr


class JunctionBase(ABC):
    @property
    @abstractmethod
    def total_width(self) -> Optional[float]:
        """ Provides the total width of the junction in meters.

        If this quantity is meaningless for the junction, it should return None.

        Returns:
            The total width of the junction or None.
        """

    @property
    @abstractmethod
    def widths(self) -> Optional[xr.DataArray]:
        """ Provides the widths of all layers the junction contains in meters.

        If this quantity is meaningless for the junction, it should return None.

        Returns:
            An xr.DataArray with the widths or None. The only coordinate must be
            'layer'.
        """

    def nk(self, wavelength: np.ndarray) -> Optional[xr.DataArray]:
        """ Provides the complex refractive index of all layers of the junction.

        This function collects the n and k data of each material of each layer and
        returns a xr.DataArray with the refractive index in complex form interpolated
        at the requested wavelengths. The output DataArray must have two coordinates,
        'wavelength' and 'layer'.

        If this quantity is meaningless for the junction, it should return None.

        Args:
            wavelength (np.ndarray): Array with the wavelengths in meters.

        Returns:
            A xr.DataArray with the complex refractive index as a function of two
            coordinates, 'wavelength' and 'layer'.
        """

    def absorptivity(
        self, wavelength: np.ndarray, angle: Optional[np.ndarray] = None
    ) -> Optional[xr.DataArray]:
        """ Provides the absorptivity of the junction.

        For those junctions not based on layers and materials, the absorptivity
        represents the fraction of light that is absorbed in the junction as a function
        of its wavelength and the angle of incidence. If the input angle is None, then
        it is assumed normal incidence.

        If this quantity is meaningless for the junction, it should return None.

        Args:
            wavelength (np.ndarray): Array with the wavelengths in meters.
            angle (Optional[np.ndarray]): Array with the angles in radians.

        Returns:
            A xr.DataArray with the absorptivity as a function of two coordinates,
            'wavelength' and 'angle'.
        """

    def solve_iv(self):
        raise NotImplementedError

    def solve_qe(self):
        raise NotImplementedError

    def solve_equilibrium(self):
        raise NotImplementedError

    def solve_short_circuit(self):
        raise NotImplementedError
