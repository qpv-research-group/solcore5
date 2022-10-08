from __future__ import annotations

import json
from pathlib import Path
from typing import Tuple, Union, Dict

import numpy as np

from pint import Quantity
from ..constants import kb_
from ..parameter import (
    InputArgumentMissing,
    Parameter,
    ParameterMissing,
    ParameterSourceBase,
    alloy_parameter,
)

ForQ = Union[float, Quantity]
"""Annotation shortcut when both a float and Quantity are valid types."""


def mobility_low_field(
    T: Quantity,
    N: Quantity,
    mu_min: Quantity,
    mu_max: Quantity,
    Nref: Quantity,
    ll: ForQ,
    t1: ForQ,
    t2: ForQ,
) -> Quantity:
    """Low field mobility model.

    This implements Eq. 4 of Sotoodeh et al.:

        mu = mu_min +   (mu_max * (300 / T) ** t1 - mu_min) / (
                        1 + (N / (Nref * (T / 300) ** t2)) ** ll
                        )

    Some trends:
        - At a given temperature, the lower the doping, the higher the mobility.
        - At a given doping:
            - If ll*t2 > t1, there's a maximum in the mobility at a certain temperature
            - Else, mobility decreases with increasing temperature

    Args:
        T: Temperature (in K)
        N: Impurity concentration (in m-3)
        mu_min: Minimum mobility (in m2/V/s)
        mu_max: Maximum mobility (in m2/V/s)
        Nref: Reference doping concentration (m-3)
        ll: Lambda exponent
        t1: Theta 1 exponent
        t2: Theta 2 exponent

    Returns:
        The low field mobility (in m2/V/s)
    """
    T = T.to("K")
    Tref = Quantity(300, "K")

    m = mu_min + (mu_max * (Tref / T) ** t1 - mu_min) / (
        1 + (N / (Nref * (T / Tref) ** t2)) ** ll
    )

    return m.to("m**2/V/s")


class SotoodehMobilitySource(ParameterSourceBase):

    name: str = "SotoodehJAP2000"
    _path: Path = (
        Path(__file__).parent.parent / "material_data" / "sotoodehJAP2000_mobility.json"
    )
    _instance = None

    def __new__(cls, *args, **kwargs):
        if cls._instance is None:
            cls._instance = ParameterSourceBase.__new__(cls)

            with cls._path.open("r") as f:
                data = json.load(f)

            cls._instance.reference = data.pop("reference", "")
            cls._instance._descriptions = data.pop("descriptions", {})
            cls._instance._data = data

        return cls._instance

    def __init__(self):
        self.reference: str
        self._descriptions: Dict[str, str]
        self._data: Dict

    @classmethod
    def load_source(cls, source_name: str = "") -> SotoodehMobilitySource:
        """Factory method to initialise the source.

        Args:
            source_name: The name of the source, needed when a general base source
                might have several concrete sources.

        Returns:
            An instance of the source class
        """
        return cls()

    @property
    def materials(self) -> Tuple[str, ...]:
        """Materials this source provides parameters for.

        Todo: At the moment, the support for alloys is a work in progress.
            Skipping those...

        Returns:
            A tuple with the list of materials.
        """
        return tuple((m for m, p in self._data.items() if "x" not in p))

    def parameters(self, material: str) -> Tuple[str, ...]:
        """Parameters available in this source for the requested material.

        Args:
            material (str): The material whose parameters are of interests.

        Returns:
            A tuple with the parameters for this material that this source provides.
        """
        if material not in self.materials:
            return ()

        return "electron_mobility", "hole_mobility"

    def get_parameter(self, material: str, parameter: str, **kwargs) -> Parameter:
        """Retrieve the parameter for the material.

        Any arguments that obtaining this parameter requires must be included as
        keyword arguments in the call.

        Args:
            material (str): Material the enquiry is about.
            parameter (str): The parameter of interest.
            **kwargs: Any other argument needed to calculate the requested parameter.

        Raises
            MaterialMissing: if the material does not exist in the source
            ParameterMissing: if the parameter does not exist for that material
            InputArgumentMissing: if there is a problem when retrieving the parameter.

        Returns:
            A Parameter object with the requested parameter.
        """
        if parameter not in self.parameters(material):
            raise ParameterMissing(self.name, material, parameter)

        carrier = parameter.split("_")[0]
        if "T" not in kwargs:
            raise InputArgumentMissing("T")

        if "Nd" not in kwargs and "Na" not in kwargs:
            raise InputArgumentMissing("Nd or Na")

        if "x" in self._data[material]:
            return self._get_parameter_alloy(material, parameter, **kwargs)

        N = max(kwargs.get("Nd", 0.0), kwargs.get("Na", 0.0))
        params = {p: Quantity(v) for p, v in self._data[material][carrier].items()}
        out = mobility_low_field(T=kwargs["T"], N=N, **params)

        return Parameter(
            out, description=self._descriptions[parameter], reference=(self.name,)
        )

    def get_nk(self, material: str, **kwargs):
        raise ParameterMissing(self.name, material, "nk")

    def _get_parameter_alloy(
        self, material: str, parameter: str, **kwargs
    ) -> Parameter:
        pass


def fraction_electrons_direct_valley(
    EgG: Quantity,
    EgX: Quantity,
    EgL: Quantity,
    mnG: Quantity,
    mnX: Quantity,
    mnL: Quantity,
    T: Quantity,
) -> Quantity:
    """Fraction of electrons in the direct conduction band valley.

    Implements equation A13 of Sutherland et al.

    References:
        "A computer analysis of heterojunction and graded composition solar cells,"
        J. E. Sutherland and J. R. Hauser, in IEEE Transactions on Electron Devices,
        24, 363 (1977), https://doi.org/10.1109/T-ED.1977.18742.

    Args:
        EgG: Gamma-valley bandgap (in eV)
        EgX: X-valley bandgap (in eV)
        EgL: L-valley bandgap (in eV)
        mnG: Gamma-valley effective mass (in kg or m0)
        mnX: X-valley effective mass (in kg or m0)
        mnL: L-valley effective mass (in kg or m0)
        T: The temperature (in K)

    Returns
        The fraction.
    """
    T = T.to("K")
    fx = (mnX / mnG) ** 1.5 * np.exp((EgG - EgX) / (kb_ * T))
    fy = (mnL / mnG) ** 1.5 * np.exp((EgG - EgL) / (kb_ * T))
    fraction = 1.0 / (1 + fx + fy)
    return fraction


def effective_indirect_mass(
    EgX: Quantity, EgL: Quantity, mnX: Quantity, mnL: Quantity, T: Quantity
) -> Quantity:
    """Effective "effective" mass considering the X and L valleys.

    It is based on equation A13 of Sutherland et al. considering only the X and L
    valleys.

    References:
        "A computer analysis of heterojunction and graded composition solar cells,"
        J. E. Sutherland and J. R. Hauser, in IEEE Transactions on Electron Devices,
        24, 363 (1977), https://doi.org/10.1109/T-ED.1977.18742.

    Args:
        EgX: X-valley bandgap (in eV)
        EgL: L-valley bandgap (in eV)
        mnX: X-valley effective mass (in kg or m0)
        mnL: L-valley effective mass (in kg or m0)
        T: The temperature (in K)

    Returns:
        The effective "effective" mass
    """
    T = T.to("K")
    fx = np.exp((EgL - EgX) / (kb_ * T))
    fraction = 1.0 / (1 + fx)
    return mnL * fraction + mnX * (1 - fraction)


def interpolate_epsilon(e1: ForQ, e2: ForQ, x: float) -> ForQ:
    """Interpolate the permittivity.

    This follows Eq. A3 of Sotoodeh et al. for epsilon

        (eo - 1) / (eo + 2) = x * (e1 - 1) / (e1 + 2) + (1 - x) * (e2 - 1) / (e2 + 2)

    References:
        "Empirical low-field mobility model for III-V compounds applicable in device
        simulation codes", M. Sotoodeh, A. H. Khalid, and A. A. Rezazadeh, Journal of
        Applied Physics 87, 2890 (2000), https://doi.org/10.1063/1.372274",


    Args:
        e1: Permittivity of material 1
        e2: Permitivity of material 2
        x: Fraction of material 1 in the alloy

    Returns:
        The interpolated permittivity
    """
    ratio1 = (e1 - 1) / (e1 + 2)
    ratio2 = (e2 - 1) / (e2 + 2)

    k = alloy_parameter(ratio2, ratio1, x)

    return (2 * k + 1) / (1 - k)


def interpolate_concentration(n1: ForQ, n2: ForQ, x: float) -> ForQ:
    """Interpolate the carrier/dopant concentrations.

    As they often vay over several orders of magnitude, their logarithm is linearly
    interpolated

        log_10(n_out) = x * log_10(n1) + (1 - x) * log_10(n2)


    Args:
        n1: Carrier/dopant concentration of material 1
        n2: Carrier/dopant concentration of material 1
        x: Fraction of material 1 in the alloy

    Returns:
        The interpolated concentration
    """
    return 10 ** alloy_parameter(np.log10(n2), np.log10(n1), x)
