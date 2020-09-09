from __future__ import annotations

from typing import Sequence, Union, Optional

import numpy as np
import xarray as xr


class Mesh1D(xr.Dataset):
    __slots__ = ()

    @classmethod
    def uniform(cls, npoints: int, nodes: Sequence[float]):
        """Uniform mesh between the smallest and largest element of the nodes.

        Either the number of points or the spacing can be provided, but not both.

        Args:
            npoints (int): Number of total points in the mesh, including the endpoints.
            nodes (Sequence[int]): The nodes of the mesh. Only the minimum and
                maximum values are used in this method.

        Returns:
            A new Mesh1D object with a uniform mesh.
        """
        return cls(
            points=np.linspace(np.min(nodes), np.max(nodes), npoints, endpoint=True),
            nodes=nodes,
        )

    @classmethod
    def piecewise_uniform(cls, npoints: Sequence[int], nodes: Sequence[float]):
        """Uniform mesh within each interval between nodes.

        Args:
            npoints (Sequence[int]): Iterable providing the number of points per
                interval. Its length should be one less than the number of nodes.
                Except for the last interval, the number of points per interval will
                exclude the endpoint.
            nodes (Sequence[int]): The nodes of the mesh.

        Returns:
            A new Mesh1D object with a uniform mesh in each interval.
        """
        if len(npoints) > 1 and len(npoints) != len(nodes) - 1:
            msg = (
                "npoints should have more than 1 element and its length be one "
                f"less than the number of nodes. Found {len(npoints)} and "
                f"{len(nodes)}"
            )
            raise ValueError(msg)

        pieces = []
        for i, n in enumerate(npoints[:-1]):
            pieces.append(np.linspace(nodes[i], nodes[i + 1], n, endpoint=False))
        pieces.append(np.linspace(nodes[-2], nodes[-1], npoints[-1], endpoint=True))

        return cls(points=np.concatenate(pieces), nodes=nodes,)

    def __init__(
        self,
        points: Union[np.ndarray, Sequence],
        nodes: Optional[Sequence[float]] = None,
    ):
        super(Mesh1D, self).__init__()

        self["points"] = xr.DataArray(points, dims="z")

        nodes = nodes if nodes is not None else [self.points.min(), self.points.max()]
        self["nodes"] = xr.DataArray(nodes, dims="z_nodes")

    def offset(self, offset: float) -> Mesh1D:
        """New Mesh1D obtained by offsetting this one.

        Args:
            offset (number): The offset.

        Returns:
            A new Mesh1D object
        """
        return Mesh1D(self.points.values + offset, self.nodes.values + offset)
