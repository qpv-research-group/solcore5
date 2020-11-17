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
        """Materials this source provides parameters for.

        Returns:
            A tuple with the list of materials.
        """
        return tuple(self._data.keys())

    def parameters(self, material: str) -> Tuple[str, ...]:
        """Parameters available in this source for the requested material.

        Args:
            material (str): The material whose parameters are of interests.

        Returns:
            A tuple with the parameters for this material that this source provides.
        """
        if material not in self.materials:
            raise ParameterError(f"Material '{material}' not in '{self.name}' source.")

        return tuple(self._data[material].keys())

    def get_parameter(self, material: str, parameter: str, **kwargs) -> Parameter:
        """Retrieve the parameter for the material.

        Any arguments that obtaining this parameter requires must be included as
        keyword arguments in the call.

        Args:
            material (str): Material the enquiry is about.
            parameter (str): The parameter of interest.
            **kwargs: Any other argument needed to calculate the requested parameter.

        Raises
            ParameterError if the material does not exist in this source or if the
                parameter does not exist for this material.

        Returns:
            A Parameter object with the requested parameter.
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
        """Retrieve the parameter for the material in the case of a ternary alloy.

        Any arguments that obtaining this parameter requires must be included as
        keyword arguments in the call.

        Args:
            material (str): Material the enquiry is about.
            parameter (str): The parameter of interest.
            **kwargs: Any other argument needed to calculate the requested parameter.

        Raises
            ParameterError if the information about any of the parents is missing.
            KeyError if the composition information is missing from the input arguments.

        Returns:
            A Parameter object with the requested parameter.
        """
        dmat = self._data[material]

        b = self.to_param(dmat.get(parameter, 0), parameter, **kwargs)
        try:
            p0 = self.get_parameter(dmat["parent0"], parameter, **kwargs)
            p1 = self.get_parameter(dmat["parent1"], parameter, **kwargs)
        except KeyError:
            raise ParameterError(
                "Ternary alloys must have 'parent0' and 'parent1' parameters defined "
                "as parent materials."
            )

        try:
            x = kwargs.get("comp", {})[dmat["x"]]
        except KeyError:
            raise KeyError(
                f"Composition for element {dmat['x']} is required to get parameter "
                f"'{parameter}' for material '{material}'."
            )

        raw = alloy_parameter(p0, p1, x, b)
        return self.to_param(raw, parameter, **kwargs)

    def to_param(
        self, raw: Union[float, str, Parameter], parameter: str, **kwargs
    ) -> Parameter:
        """Transform a raw input read from file into a Parameter object.

        If it cannot be transofrmed drectly because the parameter value is written as
        some sort of small mathematical expression (eg. 2*cos(4*T) ) then 'eval' is used
        to evaluate the expression. To be safe, the '__builtins__' passed to 'eval' are
        limited to the contents of the 'math' builtin library.

        Args:
            raw: The raw value of the parameter.
            parameter: The name of the parameter.
            **kwargs: Any other argument needed to calculate the requested parameter.

        Raises:
            NameError if 'eval' does not have all the information required to evaluate
                the expression.

        Returns:
            A Parameter object with the requested parameter.
        """
        if isinstance(raw, str):
            value, units = raw.split(" ", 1)
            try:
                value = float(value)
            except ValueError:
                value = eval(value, {"__builtins__": SAFE_BUILTINS}, kwargs)
        else:
            value = raw
            units = None

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
