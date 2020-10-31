from typing import Optional, Callable
from functools import partial

import xarray as xr


class MaterialNKDatabaseError(Exception):
    pass


class NK:
    _known_db: dict = {}
    """Dictionary with the NK data databases."""

    @staticmethod
    def register_database(
        database: Optional[Callable] = None, name: Optional[str] = None
    ):
        """Registers a database function."""
        if database is None:
            return partial(NK.register_database, name=name)

        name = name if name is not None else database.__name__
        NK._known_db[name] = database

        return database

    @staticmethod
    def get_data(
        database: str = "",
        name: str = "",
        composition: Optional[dict] = None,
        **kwargs,
    ) -> xr.DataArray:
        """Gets the complex refractive index from the database.

        Args:
            database (str): the name of the database from where to retrieve the data.
            name (str): Name of the material.
            composition (dict): Composition of the material, eg. {"In": 0.17}.
            kwargs (dict): Other keyword arguments needed by the chosen database.

        Returns:
            DataArray with the complex refractive index as function of wavelength, in m.
        """
        composition = composition if composition else {}
        if database not in NK._known_db:
            msg = (
                f"Unknown nk database '{database}'. "
                f"Available databases are '{list(NK._known_db.keys())}'."
            )
            raise MaterialNKDatabaseError(msg)

        return NK._known_db[database](name, composition, **kwargs)
