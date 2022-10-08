import json
import math
import os
import re
import sys
from functools import lru_cache
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple, Union

from pint import Quantity
from ..parameter import (
    InputArgumentMissing,
    Parameter,
    ParameterMissing,
    ParameterSourceBase,
    alloy_parameter,
)

SAFE_BUILTINS = {k: v for k, v in math.__dict__.items() if not k.startswith("__")}
"""Only common mathematical opperations are allowed when evaluating expressions."""


def locate_source_files_builtin() -> Iterator[Path]:
    """Locate the builtin parameter sources and return their names."""
    return (Path(__file__).parent.parent / "material_data").glob("*_simple_param.json")


def locate_source_files_in_path() -> Tuple[Iterator[Path], ...]:
    """Locate the source files in SOLCORE_PARAMETERS locations."""
    sep = ";" if sys.platform == "win32" else ":"
    paths = os.environ.get("SOLCORE_PARAMETERS", "").split(sep)[::-1]
    return tuple((Path(p).glob("*_simple_param.json") for p in paths if p != ""))


def locate_source_files_in_pwd() -> Iterator[Path]:
    """Locate sources in the present working directory."""
    return Path(os.getcwd()).glob("*_simple_param.json")


def locate_source_files() -> Tuple[Iterator[Path], ...]:
    """Locate the locations of the parameter sources and return their names."""
    return (
        locate_source_files_builtin(),
        *locate_source_files_in_path(),
        locate_source_files_in_pwd(),
    )


@lru_cache(maxsize=128)
def populate_sources(
    locations: Tuple[Iterator[Path], ...] = locate_source_files(),
) -> Tuple[List, Dict[str, Path], Dict[str, int]]:
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


def register_simple_sources(cls):
    """Register multiple simple sources found in standard locations."""
    name, path, priority = populate_sources()
    for n in name:
        type(
            f"{n}{cls.__name__}",
            (cls,),
            {"name": n, "_path": path[n], "_priority": priority[n]},
        )
    return cls


@register_simple_sources
class SimpleSource(ParameterSourceBase):

    name: str = "_"
    _path: Path = Path()
    _instance = None

    def __new__(cls, reference, descriptions, data, *args, **kwargs):
        if cls._instance is None:
            cls._instance = ParameterSourceBase.__new__(cls)

            cls._instance.reference = reference
            cls._instance._descriptions = descriptions
            cls._instance._data = data

        return cls._instance

    @classmethod
    def load_source(cls, source_name: str = "") -> ParameterSourceBase:
        """Factory method to initialise the source.

        Args:
            source_name: The name of the source, needed when a general base source
                might have several concrete sources.

        Returns:
            An instance of the source class
        """
        with cls._path.open("r") as f:
            data = json.load(f)

        reference = data.pop("reference", "")
        descriptions = data.pop("descriptions", {})

        return cls(reference, descriptions, data)

    def __init__(
        self,
        reference: str,
        descriptions: Optional[Dict[str, str]],
        data: Dict[str, Dict],
    ):
        """Base class for all sources directly derived from data in a file.

        Args:
            data: Dictionary of materials and their properties.
            reference: Reference indicating the origin of the data.
            descriptions: Dictionary linking each property
                (a short name) with a description indicating what they are.
        """
        self.reference: str
        self._data: Dict[str, Dict]
        self._descriptions: Dict[str, str]

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
            return ()

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
            MaterialMissing: if the material does not exist in the source
            ParameterMissing: if the parameter does not exist for that material
            InputArgumentMissing: if there is a problem when retrieving the parameter.

        Returns:
            A Parameter object with the requested parameter.
        """
        if parameter not in self.parameters(material):
            raise ParameterMissing(self.name, material, parameter)

        if "x" in self._data[material]:
            return self._get_parameter_alloy(material, parameter, **kwargs)
        else:
            raw = self._data[material][parameter]
            return self.to_param(raw, parameter, **kwargs)

    def get_nk(self, material: str, **kwargs):
        raise ParameterMissing(self.name, material, "nk")

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
            InputArgumentMissing: if the information about any of the parents is
                missing. or if the composition information is missing from the input
                arguments.

        Returns:
            A Parameter object with the requested parameter.
        """
        dmat = self._data[material]

        p0name = dmat.get("parent0", None)
        p1name = dmat.get("parent1", None)
        x = kwargs.get("comp", {}).get(dmat["x"], None)

        if p0name is None:
            raise ParameterMissing(self.name, material, "parent0")
        elif p1name is None:
            raise ParameterMissing(self.name, material, "parent1")

        if x is None:
            raise InputArgumentMissing(dmat["x"])

        kwargs[dmat["x"]] = x

        p0 = self.parman.get_parameter(p0name, parameter, **kwargs)
        p1 = self.parman.get_parameter(p1name, parameter, **kwargs)
        b = self.to_param(dmat.get(parameter, 0), parameter, **kwargs)

        raw = alloy_parameter(p0, p1, x, b)
        return self.to_param(raw, parameter, **kwargs)

    def to_param(
        self, raw: Union[float, str, Parameter], parameter: str, **kwargs
    ) -> Parameter:
        """Transform a raw input read from file into a Parameter object.

        If it cannot be transformed directly because the parameter value is written as
        some sort of small mathematical expression (eg. 2*cos(4*T) ) then 'eval' is used
        to evaluate the expression. To be safe, the '__builtins__' passed to 'eval' are
        limited to the contents of the 'math' builtin library.

        Args:
            raw: The raw value of the parameter.
            parameter: The name of the parameter.
            **kwargs: Any other argument needed to calculate the requested parameter.

        Raises:
            InputArgumentMissing: if 'eval' does not have all the information required
                to evaluate the expression.

        Returns:
            A Parameter object with the requested parameter.
        """
        if isinstance(raw, str):
            value, units = raw.split(" ", 1)
            try:
                value = float(value)
            except ValueError:
                try:
                    kw = {
                        k: (v.m if isinstance(v, Quantity) else v)
                        for k, v in kwargs.items()
                    }
                    value = eval(value, {"__builtins__": SAFE_BUILTINS}, kw)
                except NameError as err:
                    var = re.search("(?<=')(.*?)(?=')", err.__str__())
                    if var is not None:
                        raise InputArgumentMissing(var[0])
                    else:
                        raise err
        else:
            value = raw
            units = None

        return Parameter(
            value,
            units=units,
            description=self._descriptions.get(parameter, ""),
            reference=self.name,
        )


if __name__ == "__main__":
    from ..parameter import ParameterManager

    v = ParameterManager()._load_source("vurgaftmanJAP2001")
    print(v.get_parameter("GaAs", "gamma1"))
    print(v.get_parameter("GaAs", "alpha_gamma"))
    print(v.get_parameter("GaAs", "lattice_constant", T=300))
    print(v.get_parameter("GaP", "eg_gamma", T=300))
    print(v.get_parameter("GaAsSb", "valence_band_offset", T=300))
