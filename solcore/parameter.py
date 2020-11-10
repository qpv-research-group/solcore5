from __future__ import annotations

from numbers import Number
from typing import Optional, Union, Tuple, TypeVar, Dict
from abc import ABC, abstractmethod

import pint


class ParameterError(Exception):
    pass


class ParameterSourceError(Exception):
    pass


class Parameter(pint.Quantity):
    def __new__(
        cls,
        value: Union[str, Number],
        units: Optional[str] = None,
        description: Optional[str] = None,
        reference: Optional[str] = None,
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
        out._description = str(description) if description is not None else ""
        out._reference = str(reference) if reference is not None else ""
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
    def reference(self) -> str:
        """Parameter's reference. Long form for `r`"""
        return self._reference

    @property
    def r(self) -> str:
        """Parameter's reference. Short form for `reference`"""
        return self._reference

    def __str__(self) -> str:
        out = super(Parameter, self).__str__()
        if self.d != "" and self.r != "":
            return f"{self.d}: {out} ({self.r})"
        elif self.d != "":
            return f"{self.d}: {out}"
        elif self.r != "":
            return f"{out} ({self.r})"
        else:
            return out

    def __repr__(self) -> str:
        out = super(Parameter, self).__repr__()
        return out.replace(")>", f", '{self.d}', '{self.r}')>")


class ParameterSystem:

    __instance: Optional[ParameterSystem] = None

    def __new__(cls):
        if cls.__instance is None:
            inst = object.__new__(cls)
            inst._known_sources = {}
            inst.sources = {}
            ParameterSystem.__instance = inst
        return ParameterSystem.__instance

    @staticmethod
    def initialize() -> None:
        """Imports the materials data module to register the sources defined there

        Returns:
            None
        """

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
            ValueError(f"ParameterSource name '{name}' already exists.")

        self._known_sources[name] = source_class

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
            finally:
                raise ParameterError(
                    f"Parameter '{parameter}' nor found for material "
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
        self._validate_source(source)
        if source not in self.sources:
            self.sources[source] = self._known_sources[source].load_source()
        return self.sources[source]

    def _normalise_source(
        self, source: Union[str, Tuple[str], None]
    ) -> Tuple[str, ...]:
        """Normalise the source to a standard sequence of sources.

        Args:
            source (Union[str, Tuple[str], None]): The source as a string, a
                sequence of strings or None, in which case all sources should be used

        Raises:
            ParameterSourceError if any of the sources are nor registered.

        Returns:
            A tuple with the sources to check, even if it is just one.
        """
        if source is None:
            return self.known_sources
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

    name: str = ""

    def __init_subclass__(cls: ParameterSource, **kwargs):
        if cls.name == "":
            raise ValueError(
                "A ParameterSource subclass cannot have an empty string as the class "
                "attribute 'name'."
            )

        ParameterSystem().add_source(cls.name, cls)

    @property
    def parsys(self) -> ParameterSystem:
        """Convenience method to access the ParameterSystem from within a source."""
        return ParameterSystem()

    @property
    @abstractmethod
    def materials(self) -> Tuple[str, ...]:
        """Materials this source provides parameters for.

        Returns:
            A tuple with the list of materials.
        """

    @classmethod
    @abstractmethod
    def load_source(cls) -> ParameterSource:
        """Factory method to initialise the source.

        Returns:
            An instance of the source class
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

    print(ParameterSystem()._known_sources)
    ParameterSystem()._load_source("Fancy")
    print(ParameterSystem().sources)
