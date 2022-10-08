from __future__ import annotations

from typing import NamedTuple, Optional, Dict, Sequence, Callable, Union
from dataclasses import dataclass
from warnings import warn

import xarray as xr
import numpy as np
import pandas as pd

from ...junction_base import JunctionBase
from ...structure import Structure, Layer
from ...constants import kb, q, vacuum_permittivity as e0
from ...analytic_solar_cells.depletion_approximation import get_Jsrh


NK_DATA_SIGNATURE = Callable[[np.ndarray], xr.DataArray]


class LayerData(NamedTuple):
    """Parameters per layer required for a PDD simulation.

    width (float): Layer width, in m.
    bandgap (float): Bandgap of the material, in eV.
    Xi (float): Electron affinity, in eV.
    mun (float): Mobility of electrons, in m^2 / Vs.
    muh (float): Mobility of holes, in m^2 / Vs.
    Nc (float): Effective density of state in the conduction band, in m^3.
    Nv (float): Effective density of state in the valence band, in m^3.
    tn (float): Electrons minority lifetime, in s.
    th (float): Holes minority lifetime, in s.
    es (float): Relative permittivity.
    Br (float): Radiative recombination coefficient, in m^3/s.
    Cn (float): Auger recombination coefficient for electrons, in m^6/s.
    Ch (float): Auger recombination coefficient for holes, in m^6/s.
    Nd (float): Doping of the N region, in m^-3.
    Na (float): Doping of the P region, in m^-3.
    """

    width: float = 5e-6
    bandgap: float = 1.42
    Xi: float = 4.13
    mun: float = 0.47
    muh: float = 0.049
    Nc: float = 3.727e23
    Nv: float = 5.6e24
    te: float = 3e-6
    th: float = 2.5e-7
    es: float = 12.4
    Br: float = 3.53e-17
    Cn: float = 1e-42
    Ch: float = 1e-42
    Nd: float = 1.0
    Na: float = 8e22

    @classmethod
    def from_layer(cls, layer: Layer) -> LayerData:
        """Extracts the parameters required for a PDD simulation from a Layer object.

        Args:
            layer (Layer): Layer object defined by a width and a material.

        Raises:
            AttributeError: If one essential attribute is not found in the layer or the
                material.

        Returns:
            An instance of LayerData with the extracted parameters.
        """


@dataclass(frozen=True)
class PDDJunction(JunctionBase):
    """Junction class solving self-consistent PDD equations.

    Attributes:
        data (AbruptHomojunctionData): Object storing the parameters of the junction.
        nk_data (xr.DataArray, Callable[[np.ndarray], xr.DataArray]): The refractive
            index data for all layers. Can be provided as a DataArray with dimensions
            'wavelength' and 'layer' or as a function that takes as input an array of
            wavelengths and return the required DataArray.
        points (np.ndarray): Array with the discrete points across the junction. When
            creating the junction, either an integer with the total number of points,
            a list with the points per layer or the complete array can be provided. In
            the first two cases, the array of points will be created.
        layer_widths (xr.DataArray): Widths of each of the layers of the junction. If
            not provided, the information is extracted from the data attribute.
    """

    data: pd.DataFrame
    nk_data: Union[xr.DataArray, NK_DATA_SIGNATURE, None] = None
    points: Union[int, Sequence[int], np.ndarray] = 1001
    structure: Optional[Structure] = None

    def __post_init__(self):
        """Performs some checks on inputs and postprocessing."""

        if isinstance(self.nk_data, xr.DataArray) and len(self.nk_data.layer) != len(
            self.widths
        ):
            msg = (
                f"nk data must have the same length than the number of layers "
                f"along the 'layer' dimension."
            )
            raise ValueError(msg)

        object.__setattr__(self, "points", self.create_mesh(self.points))