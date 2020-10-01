from __future__ import annotations

from typing import Optional, Any, Union, Hashable, Tuple
from dataclasses import dataclass, field

import xarray as xr
import pandas as pd


class MaterialParameterError(Exception):
    pass


class ReadOnlyDict(dict):
    def __readonly__(self, *args, **kwargs):
        raise RuntimeError("Cannot modify ReadOnlyDict")

    __setitem__ = __readonly__
    __delitem__ = __readonly__
    pop = __readonly__  # type: ignore
    popitem = __readonly__
    clear = __readonly__
    update = __readonly__  # type: ignore
    setdefault = __readonly__
    del __readonly__

    def __getattr__(self, item):
        return self[item]


@dataclass(frozen=True)
class Material:
    name: str = "No name"
    composition: ReadOnlyDict = field(default_factory=ReadOnlyDict)
    T: float = 273.0
    Na: float = 0.0
    Nd: float = 0.0
    params: ReadOnlyDict = field(default_factory=ReadOnlyDict)
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
            if "wavelength" not in nk.dims or "wavelength" not in nk.coords:
                msg = "'wavelength' is not a DataArray dimension and coordinate."
                raise ValueError(msg)
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

    @classmethod
    def from_dict(cls, data: dict) -> Material:
        """Construct a material object from a plain dictionary.

        This method of creating a Material does not rely on any information contained in
        Solcore's databases and it is up to the user to make sure the data is valid
        (units, format, etc.).

        Any entry that is not "name", "T", "Na", "Nd", "composition" or "nk" is bundled
        together as "params".

        Args:
            data (dict): A dictionary with all the material information. Composition
                should be a dictionary itself, and the nk data should be either None
                or a DataArray with a coordinate being the "wavelength".

        Returns:
            A new Material object.
        """
        d = data.copy()
        result = {k: d.pop(k) for k in ("name", "T", "Na", "Nd") if k in data}
        composition = d.pop("composition") if "composition" in d else {}
        result["composition"] = ReadOnlyDict(**composition)
        result["nk"] = d.pop("nk") if "nk" in d else None

        if result["nk"] is not None and (
            "wavelength" not in result["nk"].dims
            or "wavelength" not in result["nk"].coords
        ):
            msg = (
                "'nk' is not a DataArray or 'wavelength' is not a dimension and "
                "coordinate."
            )
            raise ValueError(msg)

        result["params"] = ReadOnlyDict(**d)
        return cls(**result)

    @classmethod
    def from_dataframe(
        cls,
        data: pd.DataFrame,
        index: Hashable = 0,
        nk_cols: Union[Tuple[str, str], bool] = False,
    ) -> Material:
        """Construct a material object from a pandas DataFrame.

        This method of creating a Material does not rely on any information contained in
        Solcore's databases and it is up to the user to make sure the data is valid
        (units, format, etc.).

        Args:
            data (pd.DataFrame): A DataFrame with all the material information.
                Composition entry should be a dictionary of key: value pairs.
            index (Hashable): Index label from where to retrieve the data. By default,
                index = 0 is used.
            nk_cols (Union[Tuple[str, str], bool]): It should be either a tuple with
                the name of the columns containing the wavelength and the nk data,
                respectively; or, if True, the column names are guessed
                from the material string name and data attempted to be retrieved; or if
                False, no nk data is retrieved from the dataframe.

        Returns:
            A new Material object.
        """
        result = {k: data.loc[index, k] for k in ("name", "T", "Na", "Nd") if k in data}
        composition = data.loc[index, "composition"] if "composition" in data else {}
        result["composition"] = ReadOnlyDict(**composition)

        _nk_cols: Tuple[str, str]
        if nk_cols == True:
            name = result.get("name", "No name")
            for k, v in result["composition"].items():
                name = name.replace(k, f"{k}{v:.2}")
            wl_label = f"wavelength {name}"
            nk_label = f"nk {name}"
            _nk_cols = (wl_label, nk_label)
        elif nk_cols == False:
            result["nk"] = None
            _nk_cols = ("", "")
        elif isinstance(nk_cols, tuple):
            _nk_cols = nk_cols
        else:
            raise TypeError("'nk_cols' must be bool or a tuple with two elements.")

        if _nk_cols[0] in data.columns and _nk_cols[1] in data.columns:
            result["nk"] = xr.DataArray(
                data.loc[:, _nk_cols[1]],
                dims=["wavelength"],
                coords={"wavelength": data.loc[:, _nk_cols[0]]},
            )
        else:
            msg = f"NK data columns {_nk_cols[0]} and {_nk_cols[1]} do not exist."
            raise ValueError(msg)

        param_cols = [
            k
            for k in data.columns
            if k not in ("name", "T", "Na", "Nd", "composition") + _nk_cols
            and not k.startswith("wavelength")
            and not k.startswith("nk")
        ]
        result["params"] = ReadOnlyDict(**data.loc[index, param_cols].to_dict())
        return cls(**result)

    def __getattr__(self, item: str) -> Any:
        """Retrieve attributes stored as parameters.

        Raises:
            MaterialParameterError: If the requested attribute does not exists in the
            parameters dictionary.
        """
        if item not in self.params:
            raise MaterialParameterError(
                f"Parameter '{item}' does not exist in "
                f"material '{self.material_str}'."
            )
        return self.params[item]

    @property
    def material_str(self) -> str:
        """Return the material name embedding the composition information."""
        result = self.name
        for k, v in self.composition.items():
            result = result.replace(k, f"{k}{v:.2}")
        return result

    def to_dict(self) -> dict:
        """ Provide all the Material information as a plain dictionary.

        Returns:
            A dictionary with all the material information.
        """
        result = dict(name=self.name, T=self.T, Na=self.Na, Nd=self.Nd)
        result["composition"] = self.composition
        result.update(self.params)
        result["nk"] = self.nk
        return result

    def to_dataframe(self, include_nk=False) -> pd.DataFrame:
        """Provide all the Material information as a pandas DataFrame.

        Args:
            include_nk (bool): If the nk data should be output. If True, there will be
                two columns named 'wavelength {material_str}' and 'nk {material_str}',
                where 'material_str' is the name of the material including composition.

        Returns:
            A DataFrame with the material information.
        """
        asdict = self.to_dict()
        asdict["composition"] = [asdict["composition"]]
        asdict.pop("nk")
        data = pd.DataFrame.from_dict(asdict)

        if include_nk and self.nk is not None:
            nk = pd.DataFrame(
                {
                    f"wavelength {self.material_str}": self.nk.wavelength.data,
                    f"nk {self.material_str}": self.nk.data,
                }
            )
            return pd.concat((data, nk), axis=1)

        return data
