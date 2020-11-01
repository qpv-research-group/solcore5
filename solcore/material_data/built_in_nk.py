from pathlib import Path
from functools import lru_cache
from typing import Optional

import xarray as xr
import numpy as np

from .nk_database import NK, MaterialNKDatabaseError


@NK.register_database(name="built-in")
def get_built_in_nk_data(name: str, composition: Optional[dict] = None):
    """Gets the nk data from the built-in database.

    Args:
        name (str): Name of the material.
        composition (dict): Composition of the material, eg. {"In": 0.17}.

    Returns:
        DataArray with the complex refractive index as function of wavelength, in m.
    """
    from solcore import config

    if name not in config.materials():
        msg = (
            f"No NK data in built-in database for material '{name}'. "
            f"Available materials are {config.materials()}"
        )
        raise MaterialNKDatabaseError(msg)

    material_directory = config.materials(name)
    extension = ""
    if not composition:
        extension = ".txt"
        comp = ()
    else:
        comp = tuple(composition.values())

    n_path = Path(material_directory) / f"n{extension}"
    k_path = Path(material_directory) / f"k{extension}"

    try:
        wl, n = load_data(n_path, comp)
        _, k = load_data(k_path, comp)
    except MaterialNKDatabaseError:
        msg = (
            f"Material '{name}' does not have composition "
            f"information in the built-in database. Check material's list and "
            f"composition element, currently '{tuple(composition.keys())[0]}'."
        )
        raise MaterialNKDatabaseError(msg)

    return xr.DataArray(n + 1.0j * k, dims=["wavelength"], coords={"wavelength": wl})


@lru_cache(maxsize=128)
def load_data(path: Path, composition: tuple = ()):
    from solcore.material_system.critical_point_interpolate import (
        load_data_from_directory,
        critical_point_interpolate,
    )

    if not composition:
        return np.loadtxt(str(path), unpack=True)
    else:
        data, critical_points = load_data_from_directory(str(path))
        if not data:
            raise MaterialNKDatabaseError
        wl = list(data.values())[0][0]
        _, value, critical_points = critical_point_interpolate(
            data, critical_points, composition[0], wl
        )
        return wl, value
