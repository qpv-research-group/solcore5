from numbers import Number
from typing import Optional, Union

import pint


class Property(pint.Quantity):
    """Thin wrapper of the pint class Quantity adding a 'description' attribute."""

    def __new__(
        cls,
        value: Union[str, Number],
        units: Optional[str] = None,
        description: str = "",
    ):
        """

        Args:
            value:
            units:
            description:
        """
        v = value
        u = units
        if isinstance(value, str):
            parsed = pint.Quantity(value, units)
            v = parsed.magnitude
            u = parsed.units
        out = pint.Quantity.__new__(cls, v, u)
        out.description = description
        return out

    def __str__(self) -> str:
        out = super(Property, self).__str__()
        if self.description == "":
            return out
        return f"{self.description}: {out}"

    def __repr__(self) -> str:
        out = super(Property, self).__repr__()
        if self.description == "":
            return out
        return out.replace(")>", f", '{self.description}')>")


if __name__ == "__main__":
    a = Property("17 ms", description="Estimated time of arrival")
    b = pint.Quantity(30, "m")
    print(b / a)
