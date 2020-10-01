from pathlib import Path
from functools import lru_cache
from typing import Optional

import xarray as xr
import numpy as np

from .nk_database import NK, MaterialNKDatabaseError


@NK.register_database(name="built-in")
def get_built_in_nk_data(
    name: str, composition: Optional[dict] = None, T: float = 273.0
):
    from solcore import config

    if name not in config.materials():
        msg = (
            f"No NK data in built-in database for material '{name}'. "
            f"Available materials are {config.materials()}"
        )
        raise MaterialNKDatabaseError(msg)

    material_directory = config.materials(name)
    extension = ""
    if len(composition) == 0:
        extension = ".txt"

    n_path = Path(material_directory) / f"n{extension}"
    k_path = Path(material_directory) / f"k{extension}"

    wl, n = load_data(n_path, tuple(composition.values()))
    _, k = load_data(k_path, tuple(composition.values()))

    return xr.DataArray(n + 1.0j * k, dims=["wavelength"], coords={"wavelength": wl})


@lru_cache(maxsize=128)
def load_data(path: Path, composition: tuple = ()):
    from solcore.material_system.critical_point_interpolate import (
        load_data_from_directory,
        critical_point_interpolate,
    )

    if len(composition) == 0:
        return np.loadtxt(str(path), unpack=True)
    else:
        data, critical_points = load_data_from_directory(str(path))
        wl = list(data.values())[0][0]
        _, value, critical_points = critical_point_interpolate(
            data, critical_points, composition[0], wl
        )
        return wl, value
