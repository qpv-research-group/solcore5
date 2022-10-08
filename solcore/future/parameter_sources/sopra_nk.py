from typing import Tuple
from pathlib import Path
from configparser import ConfigParser
import pint_xarray  # noqa: F401

import xarray as xr
import pandas as pd

from ..parameter import (
    Parameter,
    ParameterSourceBase,
    InputArgumentMissing,
    ParameterMissing,
    MaterialMissing,
)


class SopraNKSource(ParameterSourceBase):

    name: str = "SopraNK"
    _root = Path(__file__).parent.parent / "material_data" / "sopra_nk"
    _instance = None

    def __new__(cls, contents: pd.DataFrame, compounds: ConfigParser, *args, **kwargs):
        if cls._instance is None:
            cls._instance = ParameterSourceBase.__new__(cls)
            cls._instance._contents = contents
            cls._instance._compounds = compounds

        return cls._instance

    def __init__(
        self, contents: pd.DataFrame, compounds: ConfigParser, *args, **kwargs
    ):
        self._contents: pd.DataFrame
        self._compounds: ConfigParser

    @classmethod
    def load_source(cls, source_name: str = "") -> ParameterSourceBase:
        """Factory method to initialise the source.

        Args:
            source_name: The name of the source, needed when a general base source
                might have several concrete sources.

        Returns:
            An instance of the source class
        """
        contents = pd.read_csv(cls._root / "SOPRA_DB_Updated.csv")
        compounds = ConfigParser()
        compounds.read(cls._root / "compounds.txt")
        return cls(contents, compounds)

    @property
    def materials(self) -> Tuple[str, ...]:
        """Materials this source provides parameters for.

        Returns:
            A tuple with the list of materials.
        """
        return tuple(
            (
                str(m)
                for m in self._contents.Symbol.unique()
                if m not in self._compounds.sections()
            )
        )

    def parameters(self, material: str) -> Tuple[str, ...]:
        """Parameters available in this source for the requested material.

        Args:
            material (str): The material whose parameters are of interests.

        Returns:
            A tuple with the parameters for this material that this source provides.
        """
        if material not in self.materials:
            return ()

        return ("nk",)

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
        raise ParameterMissing(self.name, material, parameter)

    def get_nk(self, material: str, **kwargs) -> xr.DataArray:
        """Retrieve the nk data for the material.

        Any arguments that obtaining this parameter requires must be included as
        keyword arguments in the call.

        Args:
            material (str): Material the enquiry is about.
            **kwargs: Any other argument needed to calculate the requested parameter.

        Raises
            MaterialMissing: if the material does not exist in the source
            InputArgumentMissing: if there is a problem when retrieving the parameter.

        Returns:
            A DataArray with the required data.
        """
        if material not in self.materials:
            raise MaterialMissing(self.name, material)

        if material in self._compounds.sections():
            comp = kwargs.get("comp", {})
            if len(comp) == 0:
                raise InputArgumentMissing("comp")

            x = self._compounds[material]["x"]
            c = comp.get(x, None)
            if c is None:
                raise InputArgumentMissing(x)

            data = self._load_alloy(self._root / f"{material.upper()}", comp=c)
        else:
            base_path = self._root / f"{material.upper()}.csv"
            data = pd.read_csv(base_path)

        data = xr.DataArray(
            data.n + 1.0j * data.k,
            name="nk",
            dims=["wavelength"],
            coords={"wavelength": data["wavelength[m]"]},
            attrs={"reference": self.name},
        )
        return data.pint.quantify({data.name: "dimensionless", "wavelength": "m"})

    def _load_alloy(self, folder: Path, comp: float) -> pd.DataFrame:
        """Loads the data - and calculates - the n and k of an alloy.

        Args:
            folder: Location of the input files.
            comp: Composition of the bowing element.

        Raises:
            NotImplementedError as this is not implemented yet.

        Returns:
            Dataframe with the wavelength, n and k data.
        """
        raise NotImplementedError
