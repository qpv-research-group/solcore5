import json
import math
import os
import sys
from functools import lru_cache
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple, Union

from solcore.parameter import (
    Parameter,
    ParameterError,
    ParameterSource,
    ParameterSourceBase,
    alloy_parameter,
)

SAFE_BUILTINS = {k: v for k, v in math.__dict__.items() if not k.startswith("__")}
"""Only common mathematical opperations are allowed when evaluating expressions."""


@lru_cache()
def app_dir() -> Path:
    """Finds the application data directory for the current platform."""
    if sys.platform == "win32":
        path = Path.home() / "AppData" / "Local" / "solcore"
    elif sys.platform == "darwin":
        path = Path.home() / "Library" / "ApplicationSupport" / "solcore"
    else:
        path = Path.home() / ".solcore"

    if not path.is_dir():
        raise OSError(f"Solcore's application directory '{path}' does not exist.")

    return path


def locate_source_files_builtin() -> Iterator[Path]:
    """Locate the builtin parameter sources and return their names."""
    return (Path(__file__).parent.parent / "material_data").glob("*_simple_param.json")


def locate_source_files_in_solcore_app_dir() -> Iterator[Path]:
    """Locate the source files in Solcore's application directory."""
    try:
        return app_dir().glob("*_simple_param.json")
    except OSError:
        return ()


def locate_source_files_in_pwd() -> Iterator[Path]:
    """Locate sources in the present working directory."""
    return Path(os.getcwd()).glob("*_simple_param.json")


def locate_source_files() -> Tuple[Iterator[Path], ...]:
    """Locate the locations of the parameter sources and return their names."""
    return (
        locate_source_files_builtin(),
        locate_source_files_in_solcore_app_dir(),
        locate_source_files_in_pwd(),
    )


@lru_cache
def populate_sources(
    locations: Tuple[Iterator[Path], ...] = locate_source_files(),
) -> (List, Dict[str, Path], Dict[str, int]):
    """Create a subclass of the SimpleSource for each available parameters file.

    Args:
        locations: Tuple of iterators providing the paths for the parameters files.

    Returns:
        tuple with the:
            - list of concrete sources
            - dictionary to the source files
            - dictionary with the respective priorities
    """
    name: List[str] = []
    path: Dict[str, Path] = {}
    priority: Dict[str, int] = {}
    for i, location in enumerate(locations):
        for s in location:
            n = s.stem.split("_simple_param")[0]
            name.append(n)
            path[n] = s
            priority[n] = 5 * i

    return name, path, priority


class SimpleSource(ParameterSourceBase):

    name: Union[str, List[str]] = populate_sources()[0]
    _path: Dict[str, Path] = populate_sources()[1]
    _priority: Union[int, Dict[str, int]] = populate_sources()[2]

    @classmethod
    def load_source(cls, source_name: str = "") -> ParameterSource:
        """Factory method to initialise the source.

        Args:
            source_name: The name of the source, needed when a general base source
                might have several concrete sources.

        Returns:
            An instance of the source class
        """
        path = cls._path.get(source_name, None)
        if path is None:
            raise ValueError(
                f"Impossible to load the source. 'path' for '{source_name}' is 'None'."
            )

        priority = cls._priority.get(source_name, None)
        if priority is None:
            raise ValueError(
                f"Impossible to load the source. 'priority' for '{source_name}' "
                f"is 'None'."
            )

        with path.open("r") as f:
            data = json.load(f)

        reference = data.pop("reference", "")
        descriptions = data.pop("descriptions", {})

        return cls(source_name, priority, data, reference, descriptions)

    def __init__(
        self,
        name: str,
        priority: int,
        data: Dict[str, Dict],
        reference: str = "",
        descriptions: Optional[Dict[str, str]] = None,
    ):
        """Base class for all sources directly derived from data in a file.

        Args:
            name: The concrete source name
            priority: Integer indicating the priority of the source.
            data: Dictionary of materials and their properties.
            reference: Reference indicating the origin of the data.
            descriptions: Dictionary linking each property
                (a short name) with a description indicating what they are.
        """
        self.name = name
        self._priority = priority
        self.reference = reference
        self._data = data
        self._descriptions = descriptions if descriptions is not None else {}

    @property
    def materials(self) -> Tuple[str, ...]:
        """

        Returns:

        """
        return tuple(self._data.keys())

    def parameters(self, material: str) -> Tuple[str, ...]:
        """

        Args:
            material:

        Returns:

        """
        if material not in self.materials:
            raise ParameterError(f"Material '{material}' not in '{self.name}' source.")

        return tuple(self._data[material].keys())

    def get_parameter(self, material: str, parameter: str, **kwargs) -> Parameter:
        """

        Args:
            material:
            parameter:
            **kwargs:

        Returns:

        """
        if material not in self.materials:
            raise ParameterError(f"Material '{material}' not in '{self.name}' source.")

        if "x" in self._data[material]:
            return self._get_parameter_alloy(material, parameter, **kwargs)

        else:
            if parameter not in self.parameters(material):
                raise ParameterError(
                    f"Parameter '{parameter}' not available for '{material}' in source "
                    f"'{self.name}'."
                )

            raw = self._data[material][parameter]
            return self.to_param(raw, parameter, **kwargs)

    def _get_parameter_alloy(
        self, material: str, parameter: str, **kwargs
    ) -> Parameter:
        """

        Args:
            material:
            parameter:
            **kwargs:

        Returns:

        """
        p0 = self.get_parameter(self._data[material]["parent0"], parameter, **kwargs)
        p1 = self.get_parameter(self._data[material]["parent1"], parameter, **kwargs)
        b = self.to_param(self._data[material].get(parameter, 0), parameter, **kwargs)

        try:
            x = kwargs[self._data[material]["x"]]
        except KeyError:
            raise KeyError(
                f"Composition for element {self._data[material]['x']} is "
                f"required to get parameter '{parameter}' for material "
                f"'{material}'."
            )

        raw = alloy_parameter(p0, p1, x, b)
        return self.to_param(raw, parameter, **kwargs)

    def to_param(
        self, raw: Union[float, Parameter], parameter: str, **kwargs
    ) -> Parameter:
        """

        Args:
            raw:
            parameter:
            kwargs:

        Returns:

        """
        try:
            return Parameter(
                raw,
                description=self._descriptions.get(parameter, None),
                reference=self.name,
            )
        except Exception:
            value, units = raw.split(" ", 1)
            value = eval(value, {"__builtins__": SAFE_BUILTINS}, kwargs)
            return Parameter(
                value,
                units=units,
                description=self._descriptions.get(parameter, None),
                reference=self.name,
            )


if __name__ == "__main__":
    from solcore.parameter import ParameterManager

    v = ParameterManager()._load_source("vurgaftmanJAP2001")
    print(v.get_parameter("GaAs", "gamma1"))
    print(v.get_parameter("GaAs", "alpha_gamma"))
    print(v.get_parameter("GaAs", "lattice_constant", T=300))
    print(v.get_parameter("GaP", "eg_gamma", T=300))
    print(v.get_parameter("GaAsSb", "valence_band_offset", T=300))
