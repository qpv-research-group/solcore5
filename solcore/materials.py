from __future__ import annotations

from typing import Optional, Any, Union
from dataclasses import dataclass, field

import xarray as xr


class MaterialParameterError(Exception):
    pass


class ReadOnlyDict(dict):
    def __readonly__(self, *args, **kwargs):
        raise RuntimeError("Cannot modify ReadOnlyDict")

    __setitem__ = __readonly__
    __delitem__ = __readonly__
    pop = __readonly__
    popitem = __readonly__
    clear = __readonly__
    update = __readonly__
    setdefault = __readonly__
    del __readonly__

    def __getattr__(self, item):
        return self[item]


@dataclass(frozen=True)
class Material:
    name: str = "No name"
    composition: ReadOnlyDict = field(default=ReadOnlyDict)
    T: float = 273.0
    Na: float = 0.0
    Nd: float = 0.0
    params: ReadOnlyDict = field(default=ReadOnlyDict)
    nk: Optional[xr.DataArray] = None

    @classmethod
    def factory(
        cls,
        name: str = "No name",
        composition: Optional[dict] = None,
        T: float = 273.0,
        Na: float = 0.0,
        Nd: float = 0.0,
        nk: Union[xr.DataArray, str, None] = None,
        parametric: bool = True,
        **kwargs,
    ) -> Material:
        """Create a material object out of the existing databases.

        Args:
            name (str): Name of the material.
            composition (dict): Composition of the material, eg. {"In": 0.17}. If more
                than one element is given, then the first one will become x, the
                 second y, etc.
            T (float): Temperature, in K.
            Na (float): Density of acceptors, in m^-3
            Nd (float): Density of acceptors, in m^-3
            nk (xr.DataArray, str): Either a DataArray with the complex
                refractive index as function of wavelength, in m; or the name of the
                database from where to retrieve the data.
            parametric (bool): If true (default) material parameters from the parameters
                DB will be retrieved. If this flag is set to True and the material
                is not in the parameters DB, the creation of the material will fail.
            **kwargs: Any extra argument will be incorporated to the parameters
                dictionary. If a parameter with the same name already exist, provided
                by the chosen parameters database, it will be overwritten.

        Returns:
            A new Material object.
        """
        from . import get_all_parameters, NK

        composition = composition if composition else {}

        nk_data: Optional[xr.DataArray] = None
        if isinstance(nk, str):
            nk_data = NK.get_data(nk, name, composition, T)
        elif isinstance(nk, xr.DataArray):
            nk_data = nk

        if parametric:
            params = get_all_parameters(name, composition, T, Na, Nd)
            params.update(kwargs)
        else:
            params = kwargs

        return cls(
            name=name,
            composition=ReadOnlyDict(**composition),
            T=T,
            Na=Na,
            Nd=Nd,
            params=ReadOnlyDict(**params),
            nk=nk_data,
        )

    def __getattr__(self, item: str) -> Any:
        """Retrieve attributes stored as parameters.

        Raises:
            MaterialParameterError: If the requested attribute does not exists in the
            parameters dictionary.
        """
        if item not in self.param:
            raise MaterialParameterError(
                f"Parameter '{item}' does not exist in "
                f"material '{self.material_str}'."
            )
        return self.param[item]

    @property
    def material_str(self) -> str:
        """Return the material name embedding the composition information."""
        result = self.name
        for k, v in self.composition.items():
            result = result.replace(k, f"{k}{v:.2}")
        return result
