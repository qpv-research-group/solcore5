from __future__ import annotations

from numbers import Number
from typing import Optional, Union, Tuple
from abc import ABC, abstractmethod

import pint


class ParameterError(Exception):
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
                parameter is coming from. Eg. 'VurgaftmanJAP2001'
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
            inst.raw_sources = {}
            inst.sources = {}
            ParameterSystem.__instance = inst
        return ParameterSystem.__instance

    @staticmethod
    def initialize() -> None:
        """Imports the materials data module to register the sources defined there

        Returns:
            None
        """


class ParameterSourceBase(ABC):

    name: str = ""

    def __init_subclass__(cls, **kwargs):
        if cls.name == "":
            raise ValueError(
                "A ParameterSource subclass cannot have an empty string as the class "
                "attribute 'name'."
            )

        if cls.name in ParameterSystem().raw_sources:
            raise ValueError(f"ParameterSource name '{cls.name}' already exists.")

        ParameterSystem().raw_sources[cls.name] = cls

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


if __name__ == "__main__":

    class NewSource(ParameterSourceBase):

        name = "Fancy"

        def get_parameter(self, material: str, parameter: str, **kwargs) -> Parameter:
            pass

        @property
        def materials(self):
            return ()

        def parameters(self, material: str) -> Tuple[str, ...]:
            return ()


    class AnotherNewClass(NewSource):

        name = "NotFancy"

    print(ParameterSystem().raw_sources)
