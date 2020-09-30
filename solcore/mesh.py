from __future__ import annotations

from typing import Sequence, Optional
from functools import reduce

import numpy as np
import xarray as xr


class NoMeshInformation(Exception):
    pass


class InvalidMesh(Exception):
    pass


def uniform(npoints: int, nodes: Sequence[float]):
    """Uniform mesh between the smallest and largest element of the nodes.

    Either the number of points or the spacing can be provided, but not both.

    Args:
        npoints (int): Number of total points in the mesh, including the endpoints.
        nodes (Sequence[int]): The end points of the mesh.

    Returns:
        A DataArray object with a uniform mesh.
    """
    p = np.linspace(np.min(nodes), np.max(nodes), npoints, endpoint=True)
    z = xr.DataArray(
        p,
        dims=["points"],
        coords={
            "interval": ("points", np.ones_like(p).astype(int)),
            "group": ("points", np.ones_like(p).astype(int)),
        },
    )
    return Mesh(z, nodes)


def piecewise_uniform(npoints: Sequence[int], nodes: Sequence[float]):
    """Uniform mesh within each interval between nodes.

    Args:
        npoints (Sequence[int]): Iterable providing the number of points per
            interval. Its length should be one less than the number of nodes.
        nodes (Sequence[int]): The nodes of the mesh.

    Returns:
        A new Mesh object with a uniform mesh in each interval.
    """
    if len(npoints) != len(nodes) - 1:
        msg = (
            "npoints should have more than 1 element and its length be one "
            f"less than the number of nodes. Found {len(npoints)} and "
            f"{len(nodes)}"
        )
        raise ValueError(msg)

    pieces = []
    intervals = []
    for i, n in enumerate(npoints):
        pieces.append(np.linspace(nodes[i], nodes[i + 1], n))
        intervals.append(np.full_like(pieces[-1], fill_value=i).astype(int))

    p = np.concatenate(pieces)
    z = xr.DataArray(
        p,
        dims=["points"],
        coords={
            "interval": ("points", np.concatenate(intervals)),
            "group": ("points", np.ones_like(p).astype(int)),
        },
    )
    return Mesh(z=z, nodes=nodes)


class Mesh:
    __slots__ = ("nodes", "z")

    @staticmethod
    def _validate(z: xr.DataArray) -> None:
        """Checks that the DataArray has the correct structure for being a Mesh1D.

        Args:
            z (DataArray): Array to check if it has the correct format.

        Raises:
            InvalidMesh: If the input DataArray has not the correct format.
        """
        if not all(["points" in z.dims, "interval" in z.coords, "group" in z.coords,]):
            raise InvalidMesh(
                "The format of the 'z' DataArray is not valid for a Mesh object."
            )

    def __init__(
        self, z: xr.DataArray, nodes: Optional[Sequence[float]] = None,
    ):
        if nodes is not None and (nodes[0], nodes[-1]) != (z[0].item(), z[-1].item()):
            msg = (
                "The first and last nodes must be identical to the first and "
                f"last points. Found {(nodes[0], nodes[-1])} and "
                f"{(z[0].item(), z[1].item())}, respectively."
            )
            raise ValueError(msg)

        self._validate(z)

        self.nodes = (
            np.array(nodes) if nodes is not None else np.array([min(z), max(z)])
        )
        self.z = z

    def __repr__(self) -> str:
        return self.z.__repr__() + "\nNodes:\n" + self.nodes.__repr__()

    @property
    def ninterval(self):
        """Number of available intervals."""
        return len(np.unique(self.z.interval))

    @property
    def ngroup(self):
        """Number of available groups."""
        return len(np.unique(self.z.group))

    def select_interval(self, i: int) -> Mesh:
        """Selects the mesh of the chosen interval.

        The group number of those coordinates is left unchanged. The nodes become the
        limits of the z values of the interval regardless of being exactly those before.

        Args:
            i (int): The chosen interval to select.

        Raises:
            ValueError: If the chosen interval does not exist.

        Returns:
            A Mesh with a view the z DataArray containing only the chosen interval.
        """
        if i not in self.z.interval:
            raise ValueError(
                f"The chosen interval, {i} does not exist. "
                f"Valid values are {set(self.z.interval)}"
            )

        return Mesh(z=self.z.where(self.z.interval == i, drop=True))

    def select_group(self, i: int) -> Mesh:
        """Selects the mesh of the chosen group.

        The group and interval numbers are left unchanged. As the limits of a group
        always match the limits of an interval and these are values of z by
        construction, it is guaranteed that the first and last nodes will match the
        first and last values of z.

        Args:
            i (int): The chosen group to select.

        Raises:
            ValueError: If the chosen group does not exist.

        Returns:
            A Mesh1D with a view the z DataArray containing only the chosen group.
        """
        if i not in self.z.group:
            raise ValueError(
                f"The chosen interval, {i} does not exist. "
                f"Valid values are {set(self.z.group)}"
            )

        z = self.z.where(self.z.group == i, drop=True)
        nodes = self.nodes[(self.nodes >= min(z)) * (self.nodes <= max(z))]
        return Mesh(z, nodes)

    def split_in_intervals(self) -> Sequence[Mesh]:
        """Splits the mesh in one Mesh object per interval.

        Returns:
            A tuple with one Mesh object for each interval.
        """
        return tuple(self.select_interval(i) for i in np.unique(self.z.interval))

    def split_in_groups(self) -> Sequence[Mesh]:
        """Splits the mesh in one Mesh object per group.

        Returns:
            A tuple with one Mesh object for each group.
        """
        return tuple(self.select_group(i) for i in np.unique(self.z.group))

    def broadcast(self, value, name=""):
        out = value * xr.ones_like(self.z,)
        out.name = name
        out.coords["points"] = self.z
        return out


def concat_mesh(meshes: Sequence[Mesh]) -> Mesh:
    """Concatenate 2 or more Mesh1D objects."""
    return reduce(concat2mesh, meshes[1:], meshes[0])


def concat2mesh(a: Mesh, b: Mesh) -> Mesh:
    """Concatenate the two Mesh1D objects."""
    values = np.concatenate([a.z.values, b.z.values - b.z[0].values + a.z[-1].values])

    intervals = np.concatenate(
        [a.z.interval.values, b.z.interval.values + a.z.interval.values.max() + 1,]
    )
    for i, val in enumerate(np.unique(intervals)):
        intervals[intervals == val] = i

    groups = np.concatenate(
        [a.z.group.values, b.z.group.values + a.z.group[:-1].values.max() + 1,]
    )
    for i, val in enumerate(np.unique(groups)):
        groups[groups == val] = i

    nodes = np.unique(np.concatenate([a.nodes, b.nodes - b.nodes[0] + a.nodes[-1]]))

    z = xr.DataArray(
        values,
        dims=["points"],
        name="z",
        coords={"interval": ("points", intervals), "group": ("points", groups),},
    )
    return Mesh(z, nodes)
