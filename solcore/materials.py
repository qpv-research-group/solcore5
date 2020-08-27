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


@dataclass(frozen=True)
class Material:
    name: str = "No name"
    composition: dict = field(default_factory=dict)
    T: float = 273.0
    Na: float = 0.0
    Nd: float = 0.0
    params: dict = field(default_factory=dict)
    nk: Optional[xr.DataArray] = None
    param_db: Optional[str] = None
    nk_db: Optional[str] = None

    @classmethod
    def factory(
        cls,
        name: str = "No name",
        composition: Optional[dict] = None,
        T: float = 273.0,
        Na: float = 0.0,
        Nd: float = 0.0,
        nk: Union[xr.DataArray, str, None] = None,
        param_db: Optional[str] = None,
        **kwargs,
    ) -> Material:
        """Create a material object out of the existing databases.

        Args:
            name (str): Name of the material.
            composition (dict): Composition of the material, eg. {"In": 0.17}.
            T (float): Temperature, in K.
            Na (float): Density of acceptors, in m^-3
            Nd (float): Density of acceptors, in m^-3
            nk (Optional, xr.DataArray, str): Either an DataArray with the complex
                refractive index as function of wavelength, in m; or the name of the
                database from where to retrieve the data.
            param_db (Optional, str): Name of the parameters database.
            **kwargs: Any extra argument will be incorporated to the parameters
                dictionary. If a parameter with the same name already exist, provided
                by the chosen parameters database, it will be overwritten.

        Returns:
            A new Material object.
        """
        composition = composition if composition else {}
        nk_db: Optional[str]
        nk_data: Optional[xr.DataArray]
        if isinstance(nk, str):
            nk_db = nk
            nk_data = get_nk_data(nk, name, composition, T)
        elif isinstance(nk, xr.DataArray):
            nk_db = None
            nk_data = nk
        else:
            nk_db = None
            nk_data = None

        params = (
            get_parameters(param_db, name, composition, T, Na, Nd)
            if param_db is not None
            else {}
        )
        params.update(kwargs)
        return cls(name, composition, T, Na, Nd, params, nk_data, param_db, nk_db,)

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
    param_db: str,
    name: str,
    composition: Optional[dict] = None,
    T: float = 273.0,
    Na: float = 0.0,
    Nd: float = 0.0,
) -> dict:
    """Extracts all the available parameters of the material from a database.

    Args:
        param_db (Optional, str): Name of the parameters database.
        name (str): Name of the material.
        composition (dict): Composition of the material, eg. {"In": 0.17}.
        T (float): Temperature, in K.
        Na (float): Density of acceptors, in m^-3
        Nd (float): Density of acceptors, in m^-3

    Returns:
        A dictionary with the parameters extracted from the chosen database.
    """
    composition = composition if composition else {}
    return {}


def get_nk_data(
    nk: str,
    name: str,
    composition: Optional[dict] = None,
    T: float = 273.0,
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
    InGaAs = Material(name="InGaAs", composition={"In": 0.17})
    print(InGaAs.material_str)
    print(InGaAs.cat)
