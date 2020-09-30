from __future__ import annotations

from typing import Optional, Any, Union
from dataclasses import dataclass, field

import xarray as xr


class MaterialParameterError(Exception):
    pass


class MaterialNKDatabaseError(Exception):
    pass


NK_DB: dict = {}
"""Dictionary with the NK data databases."""


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
            nk (Optional, xr.DataArray, str): Either a DataArray with the complex
                refractive index as function of wavelength, in m; or the name of the
                database from where to retrieve the data.
            **kwargs: Any extra argument will be incorporated to the parameters
                dictionary. If a parameter with the same name already exist, provided
                by the chosen parameters database, it will be overwritten.

        Returns:
            A new Material object.
        """
        composition = composition if composition else {}

        nk_data: Optional[xr.DataArray] = None
        if isinstance(nk, str):
            nk_data = get_nk_data(nk, name, composition, T)
        elif isinstance(nk, xr.DataArray):
            nk_data = nk

        params = get_parameters(name, composition, T, Na, Nd)
        params.update(kwargs)

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


def get_parameters(
    name: str,
    composition: Optional[dict] = None,
    T: float = 273.0,
    Na: float = 0.0,
    Nd: float = 0.0,
) -> dict:
    """Extracts all the available parameters of the material from a database.

    Args:
        name (str): Name of the material.
        composition (dict): Composition of the material, eg. {"In": 0.17}.
        T (float): Temperature, in K.
        Na (float): Density of acceptors, in m^-3
        Nd (float): Density of acceptors, in m^-3

    Returns:
        A ProxyDict with the parameters extracted from the chosen database.
    """
    from solcore import ParameterSystem as PS

    # First we retrieve all the available parameters
    db = PS().database[name]
    param_list = set(db.keys())
    if "x" in db:
        for key in (k for k in db if "parent" in k):
            param_list.update(set(PS().database[db[key]].keys()))

    # We include the immediate and final calculables
    param_list.update(set(PS().database["Immediate Calculables"].keys()))
    param_list.update(set(PS().database["Final Calculables"].keys()))

    # Finally, we get the parameters
    composition = composition if composition else {}
    return {
        p: PS().get_parameter(name, p, composition=composition, T=T, Na=Na, Nd=Nd)
        for p in param_list
    }


def get_nk_data(
    nk: str, name: str, composition: Optional[dict] = None, T: float = 273.0,
) -> xr.DataArray:
    """Gets the complex refractive index from the database.

    Args:
        nk (str): the name of the database from where to retrieve the data.
        name (str): Name of the material.
        composition (dict): Composition of the material, eg. {"In": 0.17}.
        T (float): Temperature, in K.

    Returns:
        DataArray with the complex refractive index as function of wavelength, in m.
    """
    composition = composition if composition else {}
    if nk not in NK_DB:
        msg = (
            f"Unknown nk database '{nk}'. "
            f"Available databases are '{list(NK_DB.keys())}'."
        )
        raise MaterialNKDatabaseError(msg)

    return NK_DB[nk](name, composition, T)


if __name__ == "__main__":
    from pprint import pprint

    print(get_parameters("GaAs"))
