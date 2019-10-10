"""
Calculates the k-points of the Brillouin zone in a given direction
"""
import numpy as np
from typing import List, Dict, Tuple

label_map: Dict[str, str] = {
    "Gamma": "$\Gamma$",
    "X": "X",
    "K": "K",
    "L": "L",
    "W": "W",
    "U": "U",
}


def brillouin_critical_points(a: float) -> Dict[str, np.ndarray]:
    return {
        "Gamma": np.array((0, 0, 0)),
        "X": np.pi / a * np.array((0, 2, 0)),
        "W": np.pi / a * np.array((1, 2, 0)),
        "L": np.pi / a * np.array((1, 1, 1)),
        "U": np.pi / a * np.array((1 / 2, 2, 1 / 2)),
        "K": np.pi / a * np.array((3 / 2, 3 / 2, 0))
    }


def traverse_brillouin(a: float,
                       traverse_order: tuple = ("L", "Gamma", "X",
                                                "W", "K", "Gamma"),
                       steps: int = 30) -> Tuple[np.ndarray, float,
                                                 List[Tuple[float, str]]]:

    critical_points = brillouin_critical_points(a)
    traverse_list = list(traverse_order)

    start_point = traverse_list.pop(0)

    k_coords = critical_points[start_point]

    coords = np.array([k_coords])

    xticks = [0]
    graph_cords = np.array([0])

    traverse_widths = []

    for i, target in enumerate(traverse_list):
        l_vec = (critical_points[target] - k_coords) * a
        traverse_widths.append(np.sqrt(sum(l_vec * l_vec)))
        k_coords = critical_points[target]

    k_coords = critical_points[start_point]

    for i, target in enumerate(traverse_list):
        coords = np.vstack(
            (coords, k_coords + np.matrix(np.linspace(0, 1, steps)[1:]).transpose() * (critical_points[target] - k_coords)))
        graph_cords = np.hstack((graph_cords, np.linspace(graph_cords[-1], graph_cords[-1] + traverse_widths[i], steps)[1:]))
        xticks.append(graph_cords[-1])
        k_coords = critical_points[target]

    scale = max(graph_cords)
    return np.array(coords), graph_cords / scale, list(zip(xticks / scale, [label_map[s] for s in traverse_order]))


def kvector(a: float, t: float = 0, p: float = np.pi,
            fraction: float = 0.2, points: int = 50,
            vin: np.ndarray = None) -> np.ndarray:
    """ Calculates the k points in a direction given by the spheric angles theta (t) and phi (p).

    The fraction of the Brilluin zone calculated is given by "fraction" and the number of points by "points".
    If "vin" is given, the direction of interest is taken from this vector."""

    s3 = np.sqrt(3)

    # Maximum k in the high symmetry directions
    Xkmax = 2 * np.pi / a
    Lkmax = s3 * np.pi / a

    # Unity vectors of the high symmetry directions (normalised)
    X = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]])
    L = 1 / s3 * np.array(
        [[1, 1, 1], [1, 1, -1], [1, -1, 1], [1, -1, -1], [-1, 1, 1], [-1, 1, -1], [-1, -1, -1], [-1, -1, 1]])

    # Unity vector in the direction of interest
    if vin is not None:
        u = vin / np.linalg.norm(vin)
    else:
        u = np.array([np.cos(t) * np.sin(p), np.sin(t) * np.sin(p), np.cos(p)])

    # We calculate the angle between u and the X and L vectors. The one producing the smallest norm will give the bounding plane.
    Xcos = Xkmax / max(abs(np.dot(X, u)))
    Lcos = Lkmax / max(abs(np.dot(L, u)))

    kmax = min(Xcos, Lcos)

    # Finally, we create the kx, ky and kz values in the u direction with the maximum norm kmax.
    magnitudes = kmax * np.linspace(0, fraction, points)
    output = []

    for i in range(len(magnitudes)):
        output.append(magnitudes[i] * u)

    return np.array(output)
