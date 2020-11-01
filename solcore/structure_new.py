from __future__ import annotations

from typing import Optional, Sequence, Tuple, Hashable, Union
from dataclasses import dataclass
from itertools import chain

import pandas as pd

from .material_new import Material


@dataclass(frozen=True)
class Layer:
    width: float
    material: Optional[Material] = None
    geometry: Optional[dict] = None

    @classmethod
    def from_dict(cls, data: dict) -> Layer:
        """Construct a Layer object from a plain dictionary.

        Internally, it calls Material.from_dict method with anything that is neither
        the layer width nor the geometry.

        Args:
            data (dict): A dictionary with all the layer information.

        Returns:
            A new Layer object.
        """
        d = data.copy()
        width = d.pop("width")
        geometry = d.pop("geometry", None)
        material = Material.from_dict(d) if d != {} else None
        return cls(width, material=material, geometry=geometry)

    @classmethod
    def from_dataframe(
        cls,
        data: pd.DataFrame,
        index: Hashable = 0,
        nk_cols: Union[Tuple[str, str], bool] = False,
    ) -> Layer:
        """Construct a layer object from a pandas DataFrame.

        Internally, it calls Material.from_dataframe method with anything that is
        neither the layer width nor the geometry.

        Args:
            data (pd.DataFrame): A DataFrame with all the layer information.
            index (Hashable): Index label from where to retrieve the data. By default,
                index = 0 is used.
            nk_cols (Union[Tuple[str, str], bool]): It should be either a tuple with
                the name of the columns containing the wavelength and the nk data,
                respectively; or, if True, the column names are guessed
                from the material string name and data attempted to be retrieved; or if
                False, no nk data is retrieved from the dataframe.

        Returns:
            A new Layer object.
        """
        width = data.loc[index, ["width"]]
        if "geometry" in data:
            geometry = data.loc[index, "geometry"]
            mat_data = data.drop(["wdith", "geometry"])
        else:
            geometry = None
            mat_data = data.drop("wdith")

        return cls(
            width,
            material=Material.from_dataframe(mat_data, index, nk_cols),
            geometry=geometry,
        )

    def to_dict(self) -> dict:
        """ Provide all the Layer information as a plain dictionary.

        Internally, it calls Material.to_dict.

        Returns:
            A dictionary with all the layer information.
        """
        result = dict(width=self.width, geometry=self.geometry)
        result.update(self.material.to_dict() if self.material is not None else {})
        return result

    def to_dataframe(self, include_nk=False) -> pd.DataFrame:
        """Provide all the Layer information as a pandas DataFrame.

        Internally, it calls Material.to_dataframe.

        Args:
            include_nk (bool): If the nk data should be output. If True, there will be
                two columns named 'wavelength {material_str}' and 'nk {material_str}',
                where 'material_str' is the name of the material including composition.

        Returns:
            A DataFrame with the layer information.
        """
        df = pd.DataFrame.from_dict(dict(width=self.width, geometry=self.geometry))

        if self.material is not None:
            mat_df = self.material.to_dataframe(include_nk)
            return pd.concat((df, mat_df), axis=1)

        return df


class Structure(Tuple[Layer, ...]):
    @classmethod
    def from_records(cls, data: Sequence[dict]) -> Structure:
        """Construct a Structure object from a sequence of dictionaries.

        Args:
            data (Sequence[dict]): A sequence of dictionaries. Each dictionary must be
            valid input for Layer.from_dict, which is called internally

        Returns:
            A new Structure object.
        """
        return cls([Layer.from_dict(d) for d in data])

    @classmethod
    def from_dataframe(
        cls,
        data: pd.DataFrame,
        nk_cols: Union[Sequence[Tuple[str, str]], bool] = False,
    ) -> Structure:
        """Construct a Structure object from a pandas DataFrame.

        Internally, it calls Layer.from_dataframe method.

        Args:
            data (pd.DataFrame): A DataFrame with all the structure information.
            nk_cols (Union[Sequence[Tuple[str, str]], bool]):
                - if False, no nk data is retrieved from the dataframe.
                - if True, the column names are guessed from the material string name
                    of each layer and data attempted to be retrieved.
                - if a Sequence of tuples, each tuple must have the name of the columns
                    containing the wavelength and the nk data, respectively. In this
                    case, there must be as many tuples as number of layers in the
                    Structure.

        Returns:
            A new Structure object.
        """
        if not nk_cols:
            return cls([Layer.from_dataframe(data, index=i) for i in data.index])
        elif isinstance(nk_cols, tuple):
            cols: list = list(nk_cols)
        else:
            cols = []
            for i in data.index:
                name = data.loc[i, "name"]
                composition = data.loc[i, "composition"]
                for k, v in composition.items():
                    name = name.replace(k, f"{k}{v:.2}")
                cols.append((f"wavelength {name}", f"nk {name}"))

        all_nk_cols = tuple(chain.from_iterable(cols))
        d = data.drop(all_nk_cols)
        return cls(
            [
                Layer.from_dataframe(
                    pd.concat([d, data[cols[i]]], axis=1), index=i, nk_cols=cols[i],
                )
                for i in d.dropna(how="all").index
            ]
        )

    def __new__(cls, layers: Sequence[Layer], name: str = "No name"):
        if not all(isinstance(layer, Layer) for layer in layers):
            raise TypeError("Layers must be of type Layer.")
        return super(Structure, cls).__new__(cls, tuple(layers))  # type: ignore

    def __init__(self, layers: Sequence[Layer], name: str = "No name"):
        self.name = name

    @property
    def widths(self) -> Tuple[float, ...]:
        return tuple(layer.width for layer in self)

    @property
    def width(self) -> float:
        return sum(self.widths)

    @property
    def relative_widths(self) -> Tuple[float, ...]:
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
