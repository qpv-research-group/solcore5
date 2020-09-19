from __future__ import annotations

from typing import Optional, Any, Union, NamedTuple, Tuple
from dataclasses import dataclass, field
from pathlib import Path
from collections import ChainMap

import xarray as xr
import toml


class MaterialParameterError(Exception):
    pass


class MaterialNKDatabaseError(Exception):
    pass


NK_DB: dict = {}
"""Dictionary with the NK data databases."""


class Composition(NamedTuple):
    """Composition of the material.

    In each field, the first component gives the element name and the second one gives
    the composition as a float number between 0 and 1. For example, if we are talking
    about InGaAs and composition has x = ("In", 0.52), it will be In0.52Ga0.48As.
    """
    x: Optional[Tuple[str, float]] = None
    y: Optional[Tuple[str, float]] = None
    z: Optional[Tuple[str, float]] = None


class ParametersDB:
    __instance = None

    def __new__(cls, data: dict):
        if ParametersDB.__instance is None:
            ParametersDB.__instance = object.__new__(cls)
            ParametersDB.__instance._data = data
        return ParametersDB.__instance

    @classmethod
    def load_databases(cls):
        from solcore import SOLCORE_ROOT

        paths = sorted((Path(SOLCORE_ROOT) / "material_data").glob("*.toml"))

        data: dict = {}
        for filename in paths:
            with filename.open("r") as f:
                data[filename.stem] = toml.load(f)

        return cls(data)

    def get_raw_parameter(self, param, name, db="all"):
        """Find the raw value of a parameter in the chosen database

        Args:
            param (str): Parameter name
            name (str): Material name
            db (str): Database to use to look for the parameter. Default: all.

        Returns:
            The raw value of the requested parameter.
        """
        return self.get_all_raw_parameters(name, db)[param]

    def get_all_raw_parameters(self, name, db="all"):
        """Find the raw value of all parameters in the chosen database.

        Args:
            name (str): Material name
            db (str): Database to use to look for the parameters. Default: all.

        Returns:
            A ChainMap object with all the raw parameters found for that material in
            the requested databases.
        """
        db_search = list(self._data.keys()) if db == "all" else (db,)
        return ChainMap(*(self._data[d].get(name, {}) for d in db_search)).new_child()


@dataclass(frozen=True)
class Material:
    name: str = "No name"
    composition: dict = field(default_factory=dict)
    T: float = 273.0
    Na: float = 0.0
    Nd: float = 0.0
    params: dict = field(default_factory=dict)
    nk: Optional[xr.DataArray] = None
    metadata: dict = field(default_factory=dict)

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
        metadata = {"param_db": param_db, "nk_db": nk_db}
        return cls(name, composition, T, Na, Nd, params, nk_data, metadata)

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

    # InGaAs = Material(name="InGaAs", composition={"In": 0.17})
    # print(InGaAs.material_str)
    # print(InGaAs.cat)

    # db = ParametersDB.load_databases()
    # GaAs = db.get_all_raw_parameters("GaAs")
    # pprint(GaAs)
    # pprint(db.get_raw_parameter("eg0_gamma", "GaAs"))
    # for k, v in GaAs.items():
    #     print(k, v)

    comp = Composition(("In", 0.52))
    pprint(comp)
