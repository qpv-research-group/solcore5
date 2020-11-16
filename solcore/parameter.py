from __future__ import annotations

from abc import ABC, abstractmethod
from functools import lru_cache
from typing import Dict, List, Optional, Tuple, TypeVar, Union

import pint


class ParameterError(Exception):
    pass


class ParameterSourceError(Exception):
    pass


class Parameter(pint.Quantity):
    def __new__(
        cls,
        value: Union[str, float, int],
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
            parsed = pint.Quantity(value, units)
            v = parsed.magnitude
            u = parsed.units
        out = pint.Quantity.__new__(cls, v, u)
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

    @staticmethod
    def gather_sources() -> None:
        """Scan several standard locations to register known sources.

        Returns:
            None
        """
        from . import parameter_system  # noqa: F401

    def initialize(self) -> None:
        """Imports all known sources

        Returns:
            None
        """
        self.gather_sources()
        for source in self.known_sources:
            self.sources[source] = self._known_sources[source].load_source(source)

    def add_source(self, name: str, source_class: ParameterSource) -> None:
        """Adds a parameters source class to the registry.

        Args:
            name (str): Name of the source
            source_class (ParameterSource): Class to register.

        Raises:
            ValueError if the source name already exists in the registry

        Returns:
            None
        """
        if name in self.known_sources:
            raise ValueError(f"ParameterSource name '{name}' already exists.")

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
        source: Union[str, Tuple[str, ...], None] = None,
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
            ParameterError if the material does not exist in the selected sources or if
            the parameter does not exist for this material.

        Returns:
            A Parameter object with the requested parameter.
        """
        nsource = self._normalise_source(source)
        for s in nsource:
            try:
                return self._load_source(s).get_parameter(material, parameter, **kwargs)
            except ParameterError:
                pass
        else:
            raise ParameterError(
                f"Parameter '{parameter}' not found for material "
                f"'{material}' in any of the following sources: {nsource}."
            )

    def get_multiple_parameters(
        self,
        material: str,
        include: Union[Tuple[str, ...], None] = None,
        exclude: Union[str, Tuple[str, ...], None] = None,
        source: Union[str, Tuple[str, ...], None] = None,
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
            ParameterError if the material does not exist in the selected sources or if
            the parameter does not exist for this material.

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
            parameters = set()
            for s in nsource:
                parameters = parameters | set(self._load_source(s).parameters(material))

        return {
            p: self.get_parameter(material, p, nsource, **kwargs)
            for p in parameters
            if p not in exclude
        }

    @lru_cache
    def _validate_source(self, source: str) -> None:
        """Checks if a source is a known source

        Args:
            source (str): The source to validate

        Raises:
            ParameterSourceError if any of the sources are nor registered.

        Returns:
            None
        """
        if source not in self.known_sources:
            msg = (
                f"Unknown source '{source}'. "
                f"Known sources are: {self.known_sources}."
            )
            raise ParameterSourceError(msg)

    def _load_source(self, source: str) -> ParameterSource:
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

    @lru_cache
    def _normalise_source(
        self, source: Union[str, Tuple[str], None]
    ) -> Tuple[str, ...]:
        """Normalise the sources to a standard sequence and prioritize them.

        Args:
            source (Union[str, Tuple[str], None]): The source as a string, a
                sequence of strings or None, in which case all sources should be used

        Raises:
            ParameterSourceError if any of the sources are nor registered.

        Returns:
            A tuple with the sources to check, even if it is just one.
        """
        if len(self.sources) == 0:
            self.initialize()

        if source is None:
            out = self.known_sources
            return tuple(
                sorted(out, reverse=True, key=lambda s: self.sources[s].priority)
            )
        elif isinstance(source, str):
            self._validate_source(source)
            return (source,)
        elif isinstance(source, Tuple):
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

    name: Union[str, List[str]] = ""
    _priority: Union[int, Dict[str, int]] = 0

    def __init_subclass__(cls: ParameterSource, **kwargs):
        if len(cls.name) == 0:
            raise ValueError(
                "A ParameterSource subclass cannot have an empty attribute 'name'."
            )
        elif isinstance(cls.name, str):
            ParameterManager().add_source(cls.name, cls)
        elif isinstance(cls.name, List):
            for n in cls.name:
                ParameterManager().add_source(n, cls)

    @classmethod
    @abstractmethod
    def load_source(cls, source_name: str = "") -> ParameterSource:
        """Factory method to initialise the source.

        Args:
            source_name: The name of the source, needed when a general base source
                might have several concrete sources.

        Returns:
            An instance of the source class
        """

    @property
    def parsys(self) -> ParameterManager:
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
            ParameterError if the material does not exist in this source or if the
                parameter does not exist for this material.

        Returns:
            A Parameter object with the requested parameter.
        """


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
    return p0 * (1.0 - x) + p1 * x - b * (1 - x) * x


ParameterSource = TypeVar("ParameterSource", bound=ParameterSourceBase)


if __name__ == "__main__":

    class NewSource(ParameterSourceBase):

        name = "Fancy"

        @classmethod
        def load_source(cls):
            return cls()

        def get_parameter(self, material: str, parameter: str, **kwargs) -> Parameter:
            pass

        @property
        def materials(self):
            return ()

        def parameters(self, material: str) -> Tuple[str, ...]:
            return ()

    class AnotherNewClass(NewSource):

        name = "NotFancy"

    print(ParameterManager()._known_sources)
    ParameterManager()._load_source("Fancy")
    print(ParameterManager().sources)
