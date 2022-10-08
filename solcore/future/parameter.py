from __future__ import annotations

from abc import ABC, abstractmethod
from functools import lru_cache
from typing import Dict, Optional, Set, Tuple, Type, Union
from warnings import warn

import xarray as xr
from pint import Quantity as Q_

from .constants import pi


class MaterialMissing(Exception):

    """Raised if a material does not exist in a source."""

    def __init__(self, source, material):
        self.source = source
        self.material = material

    def __str__(self):
        return f"Material '{self.material}' not in '{self.source}' source."


class ParameterMissing(Exception):

    """Raised if a parameter does not exist for a given material in that source."""

    def __init__(self, source, material, parameter):
        self.source = source
        self.material = material
        self.parameter = parameter

    def __str__(self):
        return (
            f"Parameter '{self.parameter}' not available for '{self.material}' in "
            f"source '{self.source}'."
        )


class InputArgumentMissing(Exception):

    """Raised if there is an error of missing input argument."""

    def __init__(self, argument):
        self.argument = argument

    def __str__(self):
        return f"Required argument '{self.argument}' could not be found in inputs."


class ParameterSourceError(Exception):

    """Raised if there is an error setting upt the source."""

    pass


class Parameter(Q_):
    def __new__(
        cls,
        value: Union[str, float, int, Q_],
        units: Optional[str] = None,
        description: str = "",
        reference: Union[str, Tuple[str, ...]] = (),
    ):
        """Wrapper of the pint.Quantity class adding 'description' and 'reference'.

        Args:
            value (str, number): Magnitude value or a string with the magnitude and its
                units.
            units (str): Units of the magnitude.
            description (str): Short description of the meaning of the parameter.
            reference (str): Keyword reference for the parameter, ie. where the
                parameter is coming from. Eg. 'Vurgaftman JAP 2001'
        """
        v = value
        u = units
        if isinstance(value, str):
            parsed = Q_(value, units)
            v = parsed.magnitude
            u = parsed.units
        elif isinstance(value, Q_):
            v = value.magnitude
            u = value.units
        out = Q_.__new__(cls, v, u)
        out._description = str(description)
        if isinstance(reference, str):
            out._reference = (reference,)
        elif isinstance(reference, tuple):
            out._reference = reference
        return out

    @property
    def description(self) -> str:
        """Parameter's description. Long form for `d`"""
        return self._description

    @property
    def d(self) -> str:
        """Parameter's description. Short form for `description`"""
        return self._description

    @property
    def reference(self) -> Tuple[str, ...]:
        """Parameter's reference. Long form for `r`"""
        return self._reference

    @property
    def r(self) -> Tuple[str, ...]:
        """Parameter's reference. Short form for `reference`"""
        return self._reference

    def __str__(self) -> str:
        out = super(Parameter, self).__str__()
        out = out.replace("Quantity", "Parameter")
        if self.d != "" and self.r != "":
            return f"{self.d}: {out} {self.r}"
        elif self.d != "":
            return f"{self.d}: {out}"
        elif self.r != "":
            return f"{out} {self.r}"
        else:
            return out

    def __repr__(self) -> str:
        out = super(Parameter, self).__repr__()
        out = out.replace("Quantity", "Parameter")
        return out.replace(")>", f", '{self.d}', '{self.r}')>")


class ParameterManager:

    _instance = None

    def __new__(cls):
        if cls._instance is None:
            inst = object.__new__(cls)
            inst._known_sources = {}
            inst.sources = {}
            cls._instance = inst
        return cls._instance

    def __init__(self):
        self._known_sources: Dict[str, Type[ParameterSourceBase]]
        self.sources: Dict[str, ParameterSourceBase]

    def gather_sources(self) -> None:
        """Scan several standard locations to register known sources.

        Returns:
            None
        """
        from . import parameter_sources  # noqa: F401

    def initialize(self) -> None:
        """Imports all known sources

        Returns:
            None
        """
        self.gather_sources()
        for source in self.known_sources:
            self.sources[source] = self._known_sources[source].load_source(source)

    def add_source(self, name: str, source_class: Type[ParameterSourceBase]) -> None:
        """Adds a parameters source class to the registry.

        Args:
            name (str): Name of the source
            source_class (ParameterSource): Class to register.

        Returns:
            None
        """
        if name in self.known_sources:
            warn(f"ParameterSource name '{name}' already exists.")

        self._known_sources[name] = source_class
        self._normalise_source.cache_clear()
        self._validate_source.cache_clear()

    @property
    def known_sources(self) -> Tuple[str, ...]:
        """Sources that have been registered and therefore are available."""
        return tuple(self._known_sources.keys())

    def get_parameter(
        self,
        material: str,
        parameter: str,
        source: Union[str, Tuple[str, ...]] = (),
        **kwargs,
    ) -> Parameter:
        """Retrieve the parameter for the material.

        Any arguments that obtaining this parameter requires must be included as
        keyword arguments in the call.

        Args:
            material (str): Material the enquiry is about.
            parameter (str): The parameter of interest.
            source (Union[str, Tuple[str], None]): Source name or list of source
                names in which to look for the information. By default, all available
                sources are used.
            **kwargs: Any other argument needed to calculate the requested parameter.

        Raises
            MaterialMissing: if the material does not exist in the selected sources
            ParameterMissing: if the parameter does not exist for that material in any
                of the selected sources
            InputArgumentMissing: if there is a problem when retrieving the parameter.

        Returns:
            A Parameter object with the requested parameter.
        """
        nsource = self._normalise_source(source)
        for s in nsource:
            if parameter in self.sources[s].parameters(material):
                return self.sources[s].get_parameter(material, parameter, **kwargs)
        else:
            raise ParameterMissing(nsource, material, parameter)

    def get_nk(
        self, material: str, source: Union[str, Tuple[str, ...]] = (), **kwargs,
    ) -> xr.DataArray:
        """Retrieve the refractive index data for the material.

        Any arguments that obtaining this parameter requires must be included as
        keyword arguments in the call.

        Args:
            material: Material the enquiry is about.
            source: Source name or list of source names in which to look for the
                information. By default, all available sources are used.
            **kwargs: Any other argument needed to calculate the requested parameter.

        Raises
            MaterialMissing: if the material does not exist in the selected sources
            ParameterMissing: if the parameter does not exist for that material in any
                of the selected sources
            InputArgumentMissing: if there is a problem when retrieving the parameter.

        Returns:
            A DataArray with the refractive index for the chosen material.
        """
        nsource = self._normalise_source(source)
        for s in nsource:
            if "nk" in self.sources[s].parameters(material):
                return self.sources[s].get_nk(material, **kwargs)
        else:
            raise ParameterMissing(nsource, material, "nk")

    def get_multiple_parameters(
        self,
        material: str,
        include: Union[Tuple[str, ...], None] = None,
        exclude: Union[str, Tuple[str, ...], None] = None,
        source: Union[str, Tuple[str, ...]] = (),
        **kwargs,
    ) -> Dict[str, Parameter]:
        """Retrieve multiple parameters for the material (defaults to all available).

        Any arguments that obtaining any of these parameter requires must be included as
        keyword arguments in the call.

        Args:
            material (str): Material the enquiry is about.
            include (Union[str, Tuple[str, ...], None] = None): The parameter or
                parameters of interest. If None, all available are retrieved.
            exclude (Union[str, Tuple[str, ...], None] = None): The parameter or
                parameters to be excluded. If None, none are excluded.
            source (Union[str, Tuple[str], None]): Source name or list of source
                names in which to look for the information. By default, all available
                sources are used.
            **kwargs: Any other argument needed to calculate the requested parameter.

        Raises
            ValueError: If both included and excluded arguments are set.
            MaterialMissing: if the material does not exist in the selected sources
            ParameterMissing: if the parameter does not exist for that material in any
                of the selected sources

        Returns:
            A dictionary of Parameter objects.
        """
        if include is not None and exclude is not None:
            raise ValueError(
                "Only 'include' or 'exclude' can be provided when asking "
                "for multiple parameters."
            )

        nsource = self._normalise_source(source)
        exclude = exclude if exclude is not None else ()
        parameters = include
        if parameters is None:
            param: Set[str] = set()
            for s in nsource:
                param = param | set(self._load_source(s).parameters(material))
            parameters = tuple(param)

        return {
            p: self.get_parameter(material, p, nsource, **kwargs)
            for p in parameters
            if p not in exclude
        }

    @lru_cache(maxsize=128)
    def _validate_source(self, source: str) -> None:
        """Checks if a source is a known source

        Args:
            source (str): The source to validate

        Raises:
            ParameterSourceError if any of the sources are not registered.

        Returns:
            None
        """
        if source not in self.known_sources:
            msg = (
                f"Unknown source '{source}'. "
                f"Known sources are: {self.known_sources}."
            )
            raise ParameterSourceError(msg)

    def _load_source(self, source: str) -> ParameterSourceBase:
        """Loads a known parameter source.

        Args:
            source (str): The name of the source to load.

        Returns:
            The instance of the loaded source.
        """
        if len(self.sources) == 0:
            self.initialize()

        self._validate_source(source)
        return self.sources[source]

    @lru_cache(maxsize=128)
    def _normalise_source(self, source: Union[str, Tuple[str]]) -> Tuple[str, ...]:
        """Normalise the sources to a standard sequence and prioritize them.

        Args:
            source: The source as a string or a sequence of strings. If it is an empty
                sequence, then all sources will be used

        Raises:
            ParameterSourceError if any of the sources are not registered.

        Returns:
            A tuple with the sources to check, even if it is just one.
        """
        if len(self.sources) == 0:
            self.initialize()

        if source == ():
            out = self.known_sources
            return tuple(
                sorted(out, reverse=True, key=lambda s: self.sources[s].priority)
            )
        elif isinstance(source, str):
            self._validate_source(source)
            return (source,)
        elif isinstance(source, tuple):
            out = ()
            for s in source:
                out = out + self._normalise_source(s)
            return out
        else:
            raise ParameterSourceError(
                f"Invalid type for source: {type(source)}. It "
                "must be a string, a tuple of strings or None."
            )


class ParameterSourceBase(ABC):

    name: str = ""
    _priority: int = 0

    def __init_subclass__(cls: Type[ParameterSourceBase]):
        if len(cls.name) == 0:
            raise ValueError(
                "A ParameterSource subclass cannot have an empty attribute 'name'."
            )
        elif not cls.name.startswith("_"):
            ParameterManager().add_source(cls.name, cls)

    @classmethod
    @abstractmethod
    def load_source(cls, source_name: str = "") -> ParameterSourceBase:
        """Factory method to initialise the source.

        Args:
            source_name: The name of the source, needed when a general base source
                might have several concrete sources.

        Returns:
            An instance of the source class
        """

    @property
    def parman(self) -> ParameterManager:
        """Convenience method to access the ParameterManager from within a source."""
        return ParameterManager()

    @property
    def priority(self) -> int:
        """Priority of the source. The higher the number, the higher the priority.

        Returns:
            The priority as a integer number
        """
        return self._priority

    @property
    @abstractmethod
    def materials(self) -> Tuple[str, ...]:
        """Materials this source provides parameters for.

        Returns:
            A tuple with the list of materials.
        """

    @abstractmethod
    def parameters(self, material: str) -> Tuple[str, ...]:
        """Parameters available in this source for the requested material.

        Args:
            material (str): The material whose parameters are of interests.

        Raises
            MaterialMissing: if the material does not exist in the source

        Returns:
            A tuple with the parameters for this material that this source provides.
        """

    @abstractmethod
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

    @abstractmethod
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
        raise ParameterMissing(self.name, material, "nk")


def alloy_parameter(
    p0: Union[Parameter, float],
    p1: Union[Parameter, float],
    x: float,
    b: Union[Parameter, float] = 0.0,
):
    """Calculate the parameter of an alloy including bowing, if present.

    If there is no bowing (b = 0), this formula becomes the standard Vegard's Law.

        result = p0 * (1 - x) + p1 * x - b * (1 - x) * x

    Args:
        p0: Value of the parameter for the first parent.
        p1: Value of the parameter for the second parent.
        x: Fractional composition of the second parent.
        b: Bowing parameter (default = 0)

    Returns:
        The parameter calculated for the alloy.
    """
    return p0 * (1 - x) + p1 * x - b * (1 - x) * x


def validate_nk(nk: xr.DataArray) -> None:
    """Validates if an nk entry is actually an xarray with the correct features"""
    if not isinstance(nk, xr.DataArray):
        raise TypeError("The value for 'nk' must be of type 'xr.DataArray'.")
    if "wavelength" not in nk.dims or "wavelength" not in nk.coords:
        msg = "'wavelength' is not a DataArray dimension and coordinate."
        raise ValueError(msg)


@xr.register_dataarray_accessor("alpha")
class Alpha:
    def __init__(self, nk: xr.DataArray):
        validate_nk(nk)
        self._nk = nk

    def __call__(self) -> xr.DataArray:
        return 4 * pi * self._nk.imag / self._nk.wavelength
