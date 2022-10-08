from __future__ import annotations

from typing import Optional, Tuple, Dict, Any

import xarray as xr
import numpy as np
from pint import Quantity as Q_
from frozendict import frozendict
from deprecated import deprecated
from warnings import warn
from functools import partial

from .parameter import ParameterManager, Parameter, validate_nk


class Material:

    __slots__ = ("name", "comp", "sources", "_nk", "_params")

    def __init__(
        self,
        name: str,
        comp: Optional[dict] = None,
        sources: Tuple[str, ...] = (),
        nk: xr.DataArray = xr.DataArray(),
        params: Optional[dict] = None,
    ):
        # Define the types
        self.name: str
        self.comp: frozendict
        self.sources: Tuple[str, ...]
        self._nk: xr.DataArray
        self._params: Dict[str, Q_]

        # Actually create the attributes. Needed this way since it is inmutable
        composition = frozendict(comp if isinstance(comp, dict) else {})
        parameters = params if isinstance(params, dict) else {}
        object.__setattr__(self, "name", name)
        object.__setattr__(self, "comp", composition)
        object.__setattr__(self, "sources", tuple(sources))
        object.__setattr__(self, "_nk", nk)
        object.__setattr__(self, "_params", parameters)

    def __getattr__(self, item: str) -> Q_:
        """Retrieve attributes for the material.

        If the requested attribute is not already in the params dictionary, then it will
        be retrieved from the available sources and the result stored in the sources.

        Raises:
            ParameterMissing: If the requested attribute does not exists in the
                parameters dictionary not in any of the available databases.

        Return:
            The value of the parameter.
        """
        if item not in self._params:
            self._params[item] = ParameterManager().get_parameter(
                material=self.name,
                parameter=item,
                source=self.sources,
                comp=self.comp,
                **self._params,
            )
        return self._params[item]

    def __setattr__(self, item, value) -> Any:
        raise AttributeError("Attributes of a Material object cannot be change.")

    @property
    def params(self) -> Tuple[str, ...]:
        """List of parameters already stored in the material."""
        return tuple(self._params.keys())

    @property
    def nk(self) -> xr.DataArray:
        """DataArray with the complex refractive index of the material.

        Raises:
            ParameterMissing: If the requested attribute does not exists in the
                parameters dictionary not in any of the available databases.

        Return:
            The refractive index
        """
        if self._nk.shape == ():
            nk = ParameterManager().get_nk(
                material=self.name, source=self.sources, comp=self.comp, **self._params,
            )
            object.__setattr__(self, "_nk", nk)
        return self._nk

    @property
    def material_str(self) -> str:
        """Return the material name embedding the composition information."""
        result = self.name
        for k, v in self.comp.items():
            result = result.replace(k, f"{k}{v:.2}")
        return result

    @classmethod
    def factory(
        cls,
        name: str,
        comp: Optional[dict] = None,
        include: Tuple[str, ...] = (),
        sources: Tuple[str, ...] = (),
        **kwargs,
    ) -> Material:
        """Create a material object out of the existing databases.

        Args:
            name: Name of the material.
            comp: Composition of the material, eg. {"In": 0.17}. If more
                than one element is given, then the first one will become x, the
                 second y, etc.
            include: Parameters to retrieve form the sources during the creation of
                the material. Any parameter can be retrieved later on, as well.
            sources: Sources in which to look for parameters. By default, all of the
                available sources are used. If provided, the sources will be scanned
                in the the order given.
            **kwargs: Any extra argument will be incorporated to the parameters
                dictionary. These have preference over the parameters built into
                Solcore. They should be Quantities, otherwise the results and errors
                might be unexpected. Common arguments are the temperature (T) and the
                doping density (Na and Nd).

        Returns:
            A new Material object.
        """
        comp = comp if comp else {}
        nk = kwargs.pop("nk", xr.DataArray())

        with_units = cls._validate_args(**kwargs)

        to_retrieve = tuple((p for p in include if p != "nk"))
        params: Dict[str, Parameter] = (
            ParameterManager().get_multiple_parameters(
                material=name,
                include=to_retrieve,
                source=sources,
                comp=comp,
                **with_units,
            )
            if to_retrieve != ()
            else {}
        )
        params.update(with_units)

        if nk.shape != ():
            validate_nk(nk)
        elif "nk" in include:
            nk = ParameterManager().get_nk(
                material=name, source=sources, comp=comp, **params
            )

        return cls(name=name, comp=comp, sources=sources, nk=nk, params=params)

    @staticmethod
    def _validate_args(**kwargs) -> Dict[str, Q_]:
        """Provide units to those arguments without them.

        T, Nd and Na are assumed to be in S.I. units (Kelvin and 1/m3). Any other is
        assumed dimensionless, which might have unexpected consequences. A warning is
        provided in this case.

        Args:
            - kwargs: Magnitudes to give units to.

        Return:
            A dictionary of the inputs including units.
        """
        result = dict()
        for k, v in kwargs.items():
            if isinstance(v, (Q_, Parameter)):
                result[k] = v
            elif k == "T":
                result[k] = Q_(v, "K")
            elif k in ("Nd", "Na"):
                result[k] = Q_(v, "1/m**3")
            else:
                result[k] = Q_(v, "dimensionless")
                warn(
                    f"Input argument '{k}' asigned 'dimensionless' units.",
                    category=UserWarning,
                )
        return result

    def __repr__(self) -> str:
        nk = self._nk.__repr__().replace("\n", " ")
        return (
            f"<Material(name={self.name}, comp={self.comp._dict}, "
            f"sources={self.sources}, nk={nk}, params={self._params})>"
        )

    @deprecated("Use 'Material.nk.real.interp' instead.", version="6.0.0")
    def n(self, wavelength: np.ndarray) -> np.ndarray:
        """Real part of the refractive index.

        Wavelength is assumed to be in meters.
        """
        return self.nk.real.interp(wavelength=Q_(wavelength, "m")).data

    @deprecated("Use 'Material.nk.real.interp' instead.", version="6.0.0")
    def k(self, wavelength: np.ndarray) -> np.ndarray:
        """Imaginary part of the refractive index or extinction coefficient.

        Wavelength is assumed to be in meters.
        """
        return self.nk.real.interp(wavelength=Q_(wavelength, "m")).data

    @deprecated("Use 'Material.nk.alpha().interp' instead.", version="6.0.0")
    def alpha(self, wavelength: np.ndarray) -> np.ndarray:
        """Absorption coefficient at the given wavelength.

        Wavelength is assumed to be in meters.
        """
        return self.nk.alpha().interp(wavelength=Q_(wavelength, "m")).data

    @deprecated("Use 'Material.nk.alpha().interp' instead.", version="6.0.0")
    def alphaE(self, energy: np.ndarray) -> np.ndarray:
        """Absorption coefficient at the given energies, given in J."""
        from .constants import h, c

        return self.alpha(wavelength=(h * c / Q_(energy, "J")))


@deprecated(version="6.0.0", reason="Use 'Material.factory' instead.")
def material(name: str):
    """Firts step of the old, 2-step interface to get materials. Deprecated.

    Used for compatibility with Solcore v5, where a 2-step material initialization
    is used. It will be removed in future releases of Solcore 6.
    """
    return partial(_material, name)


def _material(name: str, **kwargs):
    """Second step of the old, 2-step interface to get materials. Deprecated.

    Used for compatibility with Solcore v5, where a 2-step material initialization
    is used. It will be removed in future releases of Solcore 6.
    """
    comp = {k: v for k, v in kwargs.items() if k in ("In", "Al", "P", "Sb", "As", "N")}
    for k in comp.keys():
        kwargs.pop(k)
    return Material.factory(name=name, comp=comp, **kwargs)
