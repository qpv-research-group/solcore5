from __future__ import annotations

import os

from pathlib import Path
from typing import Optional, Tuple, Union
from logging import getLogger
from tempfile import TemporaryDirectory

import pint_xarray  # noqa: F401
import xarray as xr
import numpy as np

from ..parameter import (
    InputArgumentMissing,
    MaterialMissing,
    Parameter,
    ParameterMissing,
    ParameterSourceBase,
    ParameterSourceError,
)

from .tools.dboperations import Database as DB


_DATABASE_URL: str = (
    "https://refractiveindex.info/download/database/rii-database-2020-01-19.zip"
)
"""URL to the zip file containing the database."""


class RefractiveindexInfoNKSource(ParameterSourceBase):

    name: str = "RefractiveindexInfoNK"
    _instance = None

    def __new__(cls, *args, **kwargs):
        if cls._instance is None:
            cls._instance = ParameterSourceBase.__new__(cls)
            cls._instance._db = None

        return cls._instance

    def __init__(self, *args, **kwargs):
        self._db: Optional[DB]

    @classmethod
    def download_db(
        cls,
        path: Union[Path, str, None] = None,
        overwrite: bool = False,
        url: str = _DATABASE_URL,
    ) -> RefractiveindexInfoNKSource:
        """Download the database to the selected location and load it.

        If a path is not provided, it downloads it to the present working directory.

        Args:
            path: Where to download the database.
            overwrite: If it shoould be overwriten if the database exists.
            url: URL pointing to the zip file of the database to download.

        Returns:
            An instance of the source class
        """
        path = Path(path) if path is not None else Path.cwd() / "nk.db"
        if path.is_file() and not overwrite:
            getLogger().warn(
                f"NK database already exists at {str(path)}. "
                "Choose a different location or set overwrite=True."
            )
            cls.select_db(path)

        db = DB(str(path))
        with TemporaryDirectory() as output:
            db.create_database_from_url(riiurl=url, outputfolder=output)

        cls()._db = db
        return cls()

    @classmethod
    def select_db(cls, path: Union[Path, str]) -> RefractiveindexInfoNKSource:
        """Select the database file and load it.

        Args:
            path: Location of the database.

        Returns:
            An instance of the source class
        """
        if not Path(path).is_file():
            raise FileNotFoundError(f"Database file {str(path)} could not be found.")

        cls()._db = DB(str(path))
        return cls()

    @classmethod
    def load_source(cls, source_name: str = "") -> ParameterSourceBase:
        """Factory method to initialise the source.

        Args:
            source_name: The name of the source, needed when a general base source
                might have several concrete sources.

        Returns:
            An instance of the source class
        """
        cwd_db = Path.cwd() / "nk.db"
        path = (
            str(cwd_db)
            if cwd_db.is_file()
            else os.environ.get("SOLCORE_REFRACTIVE_INDEX_DB")
        )
        if path is not None:
            db: Optional[DB] = DB(path)
        else:
            msg = "No refractiveindex.info database could be found. Ignoring source."
            getLogger().warn(msg)
            db = None

        cls()._db = db
        return cls()

    @property
    def materials(self) -> Tuple[str, ...]:
        """Materials this source provides parameters for.

        Returns:
            A tuple with the list of materials.
        """
        if self._db is None:
            return ()
        else:
            return tuple(self._db.get_available().unique())

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
            ParameterSourceError: if the database it not available
            MaterialMissing: if the material does not exist in the source
            InputArgumentMissing: if there is a problem when retrieving the parameter.

        Returns:
            A DataArray with the required data.
        """
        if self._db is None:
            raise ParameterSourceError("Refractiveindex.info database not aavailable.")

        if material not in self.materials:
            raise MaterialMissing(self.name, material)

        pageid = kwargs.get("pageid", None)
        if pageid is None:
            options = self._db.get_available_for_material(material)
            msg = (
                "A 'pageid' is required to retrieve materials from the "
                "'refractiveindex.info' database. "
                f"Valid 'pageid's for material {material} are:\n{options}"
            )
            print(msg)
            raise InputArgumentMissing("pageid")

        n_data = self._db.get_material_n_numpy(int(pageid))
        k_data = self._db.get_material_k_numpy(int(pageid))
        if n_data is None:
            getLogger().warn(f"No n data for material {material}. Setting equal to 1.")
            wl = k_data[:, 0]
            n = np.ones_like(wl)
            k = k_data[:, 1]
        elif k_data is None:
            getLogger().warn(f"No k data for material {material}. Setting equal to 0.")
            wl = n_data[:, 0]
            n = n_data[:, 1]
            k = np.zeros_like(wl)
        else:
            wl = n_data[:, 0]
            n = n_data[:, 1]
            k = k_data[:, 1]

        data = xr.DataArray(
            n + 1.0j * k,
            name="nk",
            dims=["wavelength"],
            coords={"wavelength": wl},
            attrs={"reference": self.name},
        )
        return data.pint.quantify({data.name: "dimensionless", "wavelength": "m"})
