from __future__ import annotations

from typing import Optional, Type, Dict, TypeVar, Union, Sequence
from abc import ABC, abstractmethod

import numpy as np
import xarray as xr


class JunctionBase(ABC):
    def __init_subclass__(cls, **kwargs):
        """Registers any subclass of JunctionBase in the Junctions registry."""
        if cls not in JUNCTIONS_REGISTRY:
            JUNCTIONS_REGISTRY[cls.__name__] = cls
        super().__init_subclass__(**kwargs)

    @property
    def total_width(self) -> Optional[float]:
        """Provides the total width of the junction in meters.

        If this quantity is meaningless for the junction, it should return None.

        Returns:
            The total width of the junction or None."""
        return None

    @property
    def widths(self) -> Optional[xr.DataArray]:
        """Provides the widths of all layers the junction contains in meters.

        If this quantity is meaningless for the junction, it should return None.

        Returns:
            An xr.DataArray with the widths or None. The only coordinate must be
            'layer'."""
        return None

    def nk(self, wavelength: np.ndarray) -> Optional[xr.DataArray]:
        """Provides the complex refractive index of all layers of the junction.

        This function collects the n and k data of each material of each layer and
        returns a xr.DataArray with the refractive index in complex form interpolated
        at the requested wavelengths. The output DataArray must have two coordinates,
        'wavelength' and 'layer'.

        If this quantity is meaningless for the junction, it should return None.

        Args:
            wavelength (np.ndarray): Array with the wavelengths in meters.

        Returns:
            A xr.DataArray with the complex refractive index as a function of two
            coordinates, 'wavelength' and 'layer'."""
        return None

    def absorptivity(
        self, wavelength: np.ndarray, angle: Optional[np.ndarray] = None
    ) -> Optional[xr.DataArray]:
        """Provides the absorptivity of the junction.

        For those junctions not based on layers and materials, the absorptivity
        represents the fraction of light that is absorbed in the junction as a function
        of its wavelength and the angle of incidence. If the input angle is None, then
        it is assumed normal incidence.

        If this quantity is meaningless for the junction, it should return None.

        Args:
            wavelength (np.ndarray): Array with the wavelengths in meters.
            angle (np.ndarray, optional): Array with the angles in radians.

        Returns:
            A xr.DataArray with the absorptivity as a function of two coordinates,
            'wavelength' and 'angle'."""
        return None

    @abstractmethod
    def solve_iv(
        self,
        voltage: np.ndarray,
        absorption: Optional[xr.DataArray] = None,
        source: Optional[xr.DataArray] = None,
    ) -> xr.Dataset:
        """Calculates the IV curve of the junction.

        If absorption is provided, then light_source must also be provided and the
        light IV curve should be calculated instead. In this case, parameters like
        Voc, Isc, fill factor, etc. are also calculated.

        Args:
            voltage (np.ndarray): Array of voltages at which to calculate the IV curve.
            absorption (xr.DataArray, optional): Array with the fraction of absorbed
                light as a function of 'wavelength' and 'position'.
            source (xr.DataArray, optional): Light source to use providing the number
                of photons as a junction of 'wavelength'.

        Returns:
            A xr.Dataset with the output of the calculation. At a minimum, it must
            contain a 'current' DataArray giving the current in amps as a function of
            the input 'voltage'. If light IV is calculated, the curve parameters (Voc,
            Isc, FF, Vmpp, Impp and Pmpp) must be provided as attributes of the
            Dataset.
            Other DataArrays containing extra information resulting from the calculation
            might be available, depending on the junction."""

    @abstractmethod
    def solve_qe(
        self, absorption: xr.DataArray, source: Optional[xr.DataArray] = None
    ) -> xr.Dataset:
        """Calculates the external and internal quantum efficiency of the junction.

        Args:
            absorption (xr.DataArray, optional): Array with the fraction of absorbed
                light as a function of 'wavelength' and 'position'.
            source (xr.DataArray, optional): Light source to use providing the number
                of photons as a junction of 'wavelength'.

        Returns:
            A xr.Dataset with the output of the calculation. At a minimum, it must
            contain 'eqe' and 'iqe' DataArrays giving the external and internal quantum
            efficiencies, respectively, as a function of 'wavelength'.
            Other DataArrays containing extra information resulting from the calculation
            might be available, depending on the junction."""

    @abstractmethod
    def solve_equilibrium(self):
        """Calculates the junction band structure at equilibrium.

        Returns:
            A xr.Dataset with the output of the calculation. At a minimum, it must
            contain 'conduction_band' and 'valence_band' DataArrays giving the
            conduction and valence band profiles, respectively, as a function of
            'position'.
            Other DataArrays containing extra information resulting from the calculation
            might be available, depending on the junction."""

    @abstractmethod
    def solve_short_circuit(
        self, absorption: xr.DataArray, source: xr.DataArray,
    ) -> xr.Dataset:
        """Calculates the junction band structure at short circuit.

        Args:
            absorption (xr.DataArray): Array with the fraction of absorbed
                light as a function of 'wavelength' and 'position'.
            source (xr.DataArray): Light source to use providing the number
                of photons as a junction of 'wavelength'.

        Returns:
            A xr.Dataset with the output of the calculation. At a minimum, it must
            contain 'conduction_band' and 'valence_band' DataArrays giving the
            conduction and valence band profiles, respectively, as a function of
            'position'.
            Other DataArrays containing extra information resulting from the calculation
            might be available, depending on the junction."""

    @staticmethod
    def iv_parameters(voltage: np.ndarray, current: np.ndarray) -> dict:
        """Calculates Voc, Jsc, FF, Vmpp, Jmpp and Pmpp.

        Args:
            voltage (np.ndarray):
            current (np.ndarray):

        Returns:
            A dictionary with the calculated parameters.
        """
        return iv_parameters(voltage, current)

    def create_mesh(self, points: Union[int, Sequence[int], np.ndarray]) -> np.ndarray:
        """Creates a 1D array of points, serving as spacial mesh for the calculations.

        Args:
            points (int, Sequence[int], np.ndarray): One of the following:
                - Total number of points for the junction.
                - A sequence with the number of points per layer.
                - An array with the mesh already created.

        Raises:
            ValueError: If 'total_width' is None.
            ValueError: If 'widths' is None when points is a Sequence.
            ValueError: If 'widths' and 'points' have different length when points
                is a Sequence.
            ValueError: If points is and array and the first and last points are not
                0 and 'total_width' respectively.

        Returns:
            An array with the mesh created according to the rules above.
        """
        if self.total_width is None:
            msg = "Total junction width cannot be None when creating a mesh."
            raise ValueError(msg)

        if isinstance(points, int):
            return np.linspace(0, self.total_width, points, endpoint=True)

        elif isinstance(points, np.ndarray):
            if points[0] != 0.0 or not np.isclose(points[-1], self.total_width):
                msg = (
                    "First and last elements of the points array must be 0 and "
                    "the total junction width, respectively."
                )
                raise ValueError(msg)
            return points

        elif isinstance(points, Sequence):
            if self.widths is None:
                msg = (
                    "Layer widths cannot be None when creating a mesh from a "
                    "sequence."
                )
                raise ValueError(msg)
            if len(self.widths) != len(points):
                msg = (
                    f"'widths' and 'points' have different length: "
                    f"{len(self.widths)} and {len(points)}."
                )
                raise ValueError(msg)
            w = np.cumsum([0] + list(self.widths))
            return np.append(
                np.concatenate(
                    [
                        np.linspace(start, stop, p, endpoint=False)
                        for start, stop, p in zip(w[:-1], w[1:], points)
                    ]
                ),
                self.total_width,
            )


Junction = TypeVar("Junction", bound=JunctionBase)

JUNCTIONS_REGISTRY: Dict[str, Type[Junction]] = {}
"""Registry of all junctions available in Solcore."""


def iv_parameters(voltage: np.ndarray, current: np.ndarray) -> Dict[str, float]:
    """Calculates Voc, Jsc, FF, Vmpp, Jmpp and Pmpp.

    Args:
        voltage (np.ndarray):
        current (np.ndarray):

    Returns:
        A dictionary with the calculated parameters.
    """
    i0 = np.abs(voltage).argmin()
    iVoc = np.abs(current).argmin()
    step = 1 if iVoc > i0 else -1
    q = slice(i0, iVoc, step)
    power = np.abs(voltage * current)[q]

    Jsc = current[i0]
    Voc = voltage[iVoc]
    Pmmp = power.max()
    Jmpp = current[q][power.argmax()]
    Vmpp = voltage[q][power.argmax()]
    FF = Pmmp / (Jsc * Voc)

    return dict(Jsc=Jsc, Voc=Voc, Pmmp=Pmmp, Jmpp=Jmpp, Vmpp=Vmpp, FF=FF)
