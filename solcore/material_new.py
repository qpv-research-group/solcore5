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

    @classmethod
    def from_dict(cls, data: dict) -> Material:
        """Construct a material object from a plain dictionary.

        This method of creating a Material does not rely on any information contained in
        Solcore's databases and it is up to the user to make sure the data is valid
        (units, format, etc.).

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
        result["params"] = ReadOnlyDict(**d)
        return cls(**result)

    @classmethod
    def from_dataframe(
        cls,
        data: pd.DataFrame,
        index: Hashable = 0,
        nk_cols: Optional[Tuple[str, str]] = None,
    ) -> Material:
        """Construct a material object from a pandas DataFrame.

        This method of creating a Material does not rely on any information contained in
        Solcore's databases and it is up to the user to make sure the data is valid
        (units, format, etc.).

        Args:
            data (pd.DataFrame): A DataFrame with all the material information.
                Composition should be a dictionary itself.
            index (Hashable): Index label from where to retrieve the data. By default,
                index = 0 is used.
            nk_cols (Optional[Tuple[str, str]]): If provided, it should be a tuple with
                the name of the column containing the wavelength and the one containing
                the nk data, respectively. If not given, the column names are guessed
                from the material string name.

        Returns:
            A new Material object.
        """
        result = {k: data.loc[index, k] for k in ("name", "T", "Na", "Nd") if k in data}
        composition = data.loc[index, "composition"] if "composition" in data else {}
        result["composition"] = ReadOnlyDict(**composition)

        if nk_cols is None and "name" in result:
            name = result["name"]
            for k, v in result["composition"].items():
                name = name.replace(k, f"{k}{v:.2}")
            wl_label = f"wavelength{name}"
            nk_label = f"nk{name}"
            nk_cols = (wl_label, nk_label)
        else:
            nk_cols = ("", "")

        if nk_cols[0] in data.columns and nk_cols[1] in data.columns:
            result["nk"] = xr.DataArray(
                data.loc[:, nk_cols[1]],
                dims=["wavelength"],
                coords={"coords": data.loc[:, nk_cols[0]]},
            )
        else:
            result["nk"] = None

        param_cols = [
            k
            for k in data.cols
            if k not in ("name", "T", "Na", "Nd", "composition") + nk_cols
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

    def as_dict(self) -> dict:
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
        asdict = self.as_dict()
        asdict.pop("nk")
        if include_nk and self.nk is not None:
            asdict[f"wavelength{self.material_str}"] = self.nk.wavelength.data
            asdict[f"nk{self.material_str}"] = self.nk.data

        return pd.DataFrame.from_dict(asdict)
