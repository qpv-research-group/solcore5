from typing import Tuple, Optional, Iterator, Dict, Union
from pathlib import Path
import json
import math
import logging

from solcore.parameter import (
    Parameter,
    ParameterSource,
    ParameterSourceBase,
    ParameterError,
    alloy_parameter,
)


SAFE_BUILTINS = {k: v for k, v in math.__dict__.items() if not k.startswith("__")}
"""Only common mathematical opperations are allowed when evaluating expressions."""


def locate_builtin() -> Iterator[Path]:
    """Locate the builtin parameter sources and return their names."""
    names = (Path(__file__).parent.parent / "material_data").glob("*_parameters.*")
    return (filename for filename in names if "mobility" not in filename.stem)


def populate_sources(cls: ParameterSource, builtins: Iterator[Path] = locate_builtin()):
    """Create a subclass of the BuiltInBaseSource for each available parameters file.

    Args:
        cls (cls): Class to subclass to create the builtin sources
        builtins (Iterator[Path]): Iterator providing the paths for the builtin
            parameters files.

    Returns:
        The same input class
    """
    for s in builtins:
        name = s.stem.split("_parameters")[0]
        type(f"{name.capitalize()}Source", (cls,), {"name": name, "path": s})

    return cls


@populate_sources
class BuiltInBaseSource(ParameterSourceBase):

    name = "_builtins"
    path: Optional[Path] = None

    @classmethod
    def load_source(cls) -> ParameterSource:
        """

        Returns:

        """
        if cls.path is None:
            raise ValueError(
                "Impossible to load the source. Class attribute 'path' is 'None'."
            )

        with cls.path.open("r") as f:
            data = json.load(f)

        reference = data.pop("reference", "")
        descriptions = data.pop("descriptions", {})

        return cls(data, reference, descriptions)

    def __init__(
        self,
        data: Dict[str, Dict],
        reference: str = "",
        descriptions: Optional[Dict[str, str]] = None,
    ):
        """Base class for all sources directly derived from data in a file.

        Args:
            data (Dict[str, Dict]): Dictionary of materials and their properties.
            reference (str): Reference indicating the origin of the data.
            descriptions (Dict[str, str]): Dictionary linking each property
                (a short name) with a description indicating what they are.
        """
        self._data = data
        self.reference = reference
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
            raise KeyError(f"Composition for element {self._data[material]['x']} is "
                           f"required to get parameter '{parameter}' for material "
                           f"'{material}'.")

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
    from solcore.parameter import ParameterSystem

    v = ParameterSystem()._load_source("vurgaftman_JAP_2001")
    print(v.get_parameter("GaAs", "gamma1"))
    print(v.get_parameter("GaAs", "alpha_gamma"))
    print(v.get_parameter("GaAs", "lattice_constant", T=300))
    print(v.get_parameter("GaP", "eg_gamma", T=300))
    print(v.get_parameter("GaAsSb", "valence_band_offset", T=300))
