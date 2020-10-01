from __future__ import annotations

from typing import Optional, Sequence, Tuple
from dataclasses import dataclass

from .material_new import Material


@dataclass(frozen=True)
class Layer:
    width: float
    name: str = "No name"
    material: Optional[Material] = None
    geometry: Optional[dict] = None


class Structure(Tuple[Layer, ...]):
    def __new__(cls, layers: Sequence[Layer], name: str = "No name"):
        if not all(isinstance(layer, Layer) for layer in layers):
            raise TypeError("Layers must be of type Layer.")
        return super(Structure, cls).__new__(cls, tuple(layers))

    def __init__(self, layers: Sequence[Layer], name: str = "No name"):
        self.name = name

    @property
    def widths(self) -> Tuple[float]:
        return tuple(layer.width for layer in self)

    @property
    def width(self) -> float:
        return sum(self.widths)

    @property
    def relative_widths(self) -> Tuple[float]:
        return tuple(layer.width / self.width for layer in self)

    def insert_layer(
        self, idx: int, layer: Layer, name: Optional[str] = None
    ) -> Structure:
        if not isinstance(layer, Layer):
            raise TypeError("Layers must be of type Layer.")
        name = name if name is not None else self.name
        return Structure(self[:idx] + (layer,) + self[idx + 1 :], name)


if __name__ == "__main__":
    foo = Structure([Layer(10), Layer(20)], "Cat")
    print(foo)
    print(foo.widths)
    print(foo.width)
    print(foo.relative_widths)
