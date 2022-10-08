from typing import Tuple, Dict
from pathlib import Path
from functools import lru_cache
import pint_xarray  # noqa: F401

import numpy as np
import xarray as xr

from ..parameter import (
    Parameter,
    ParameterSourceBase,
    InputArgumentMissing,
    ParameterMissing,
    MaterialMissing,
)


class BuiltinNKSource(ParameterSourceBase):

    name: str = "BuiltinNK"
    _root = Path(__file__).parent.parent / "material_data" / "builtin_nk"
    _instance = None

    def __new__(cls, paths, *args, **kwargs):
        if cls._instance is None:
            cls._instance = ParameterSourceBase.__new__(cls)
            cls._instance._paths = paths

        return cls._instance

    def __init__(self, paths, *args, **kwargs):
        self._paths: Dict[str, Path]

    @classmethod
    def load_source(cls, source_name: str = "") -> ParameterSourceBase:
        """Factory method to initialise the source.

        Args:
            source_name: The name of the source, needed when a general base source
                might have several concrete sources.

        Returns:
            An instance of the source class
        """
        paths = {p.stem.split("_nk")[0]: p for p in cls._root.glob("*_nk")}
        return cls(paths)

    @property
    def materials(self) -> Tuple[str, ...]:
        """Materials this source provides parameters for.

        Returns:
            A tuple with the list of materials.
        """
        return tuple(self._paths.keys())

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

        base_path = self._paths[material]
        if len(material) > 4:
            comp = kwargs.get("comp", None)
            if comp is None or len(comp) == 0:
                raise InputArgumentMissing("comp")

            c = tuple(comp.values())[0]
            wl, n = self._load_alloy(path=base_path, param="n", comp=c)
            _, k = self._load_alloy(path=base_path, param="k", comp=c)
        else:
            wl, n = self._load(path=base_path, param="n")
            _, k = self._load(path=base_path, param="k")

        data = xr.DataArray(
            n + 1.0j * k,
            name="nk",
            dims=["wavelength"],
            coords={"wavelength": wl},
            attrs={"reference": self.name},
        )
        return data.pint.quantify({data.name: "dimensionless", "wavelength": "m"})

    @lru_cache(maxsize=128)
    def _load(self, path: Path, param: str) -> Tuple[np.ndarray, np.ndarray]:
        """Loads the refractive index data from a txt file.

        The location of the file to load is expected to be

            BUILTIN_ROOT/{material}/{param}.txt

        Args:
            path: The root path to look for the file.
            param: Either 'n' or 'k', depending on which index is required.

        Returns:
            Tuple with the wavelength (m) and the refractive index (dimensionless)
        """
        return np.loadtxt(str(path / f"{param}.txt"), unpack=True)

    @lru_cache(maxsize=128)
    def _load_alloy(
        self, path: Path, param: str, comp: float
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Loads the data and then interpolates it to the required composition.

        The location of the files to load is expected to be in the subfolder:

            BUILTIN_ROOT/{material}/{param}

        and then inside

            - A file for each composition X.XXX with name: X.XXX_{material}_{param}.txt
            - A file critical_points.txt indicating the location of each critical point
                for each composition

        Args:
            path: The root path to look for the file.
            param: Either 'n' or 'k', depending on which index is required.
            comp: The composition of interest.

        Returns:
            Tuple with the wavelength (m) and the refractive index (dimensionless)
        """
        from ..parameter_sources.tools.critical_point_interpolate import (
            load_data_from_directory,
            critical_point_interpolate,
        )

        data, critical_points = load_data_from_directory(str(path / param))
        wl = list(data.values())[0][0]
        _, value, critical_points = critical_point_interpolate(
            data, critical_points, comp, wl
        )
        return wl, value
