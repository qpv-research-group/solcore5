import types
from typing import Optional, Tuple

import numpy as np
from numpy.typing import NDArray

try:
    from ..absorption_calculator import (
        calculate_absorption_profile_rcwa,
        calculate_rat_rcwa,
    )

    reason_to_exclude = None
except ImportError:
    reason_to_exclude = "An installation of S4 has not been found."

from ..registries import register_optics
from ..solar_cell import SolarCell
from ..state import State
from ..structure import Junction, Layer, TunnelJunction

rcwa_options = State()
rcwa_options.size = ((500, 500), (500, 500))
rcwa_options.orders = 4
rcwa_options.theta = 0
rcwa_options.phi = 0
rcwa_options.pol = "u"
rcwa_options.parallel = False
rcwa_options.n_jobs = -1


@register_optics(name="RCWA", reason_to_exclude=reason_to_exclude)
def solve_rcwa(
    solar_cell: SolarCell,
    wavelength: NDArray,
    position: NDArray,
    parallel: bool = False,
    size: Tuple[int, int] = ((500, 500), (500, 500)),
    orders: int = 4,
    theta: float = 0,
    phi: float = 0,
    pol: str = "u",
    n_jobs: int = -1,
    rcwa_options: Optional[State] = None,
    **kwargs
) -> None:
    """Calculates the RAT of a solar cell object using the RCWA solver.

    The SolarCell object is updated with the wavelength, the calculated reflected,
    transmitted and absorved ligth. Additionally, each a method to calculate the
    absorved ligth per junction and the ligth absorved per later are also added to the
    individual junctions.

    Args:
        solar_cell (SolarCell):  A solar_cell object
        wavelength (NDArray): Array of wavelegth at which the optics are calculated.
        position (NDArray): Array of positions in the z direction to calculate the
            absorption vs depth.
        parallel (bool, optional): whether or not to execute calculation in parallel
            (over wavelengths). Defaults to False.
        size (Tuple[int, int], optional):  a tuple of 2-D vectors in the format ((ux,
            uy), (vx, vy)) giving the x and y components of the lattice unit vectors in
            nm. Defaults to ((500, 500), (500, 500)).
        orders (int, optional): number of orders to retain in the RCWA calculations.
            Defaults to 4.
        theta (float, optional): the polar incidence angle, in degrees, with 0 degrees
            being normal incidence. Defaults to 0.
        phi (float, optional): azimuthal incidence angle (in degrees). Defaults to 0.
        pol (str, optional): the polarization of the light ('s', 'p' or 'u'). Defaults
            to "u".
        n_jobs (int, optional): the 'n_jobs' argument passed to Parallel from the joblib
            package. If set to -1, all available CPUs are used, if set to 1 no parallel
            computing is executed. The number of CPUs used is given by n_cpus + 1 +
            n_jobs. Defaults to -1.
        rcwa_options (Optional[State], optional): dictionary of options for S4. Defaults
            to None. The list of possible entries and their values is:

            * LatticeTruncation: 'Circular' or 'Parallelogramic' (default 'Circular')
            * DiscretizedEpsilon: True or False (default False)
            * DiscretizationResolution: integer (default value 8)
            * PolarizationDecomposition: True or False (default False)
            * PolarizationBasis: 'Default' or 'Normal' or 'Jones' (default 'Default')
            * LanczosSmoothing: True or False (default False)
            * SubpixelSmoothing: True or False (default False)
            * ConserveMemory: True or False (default False)
            * WeismannFormulation: True or False (default False)

    Return:
        None
    """
    solar_cell.wavelength = wavelength

    # We include the shadowing losses
    initial = (1 - solar_cell.shading) if hasattr(solar_cell, "shading") else 1

    # Now we calculate the absorbed and transmitted light. We first get all the relevant
    # parameters from the objects
    all_layers = []
    n_layers_junction = []
    for j, layer_object in enumerate(solar_cell):

        # Attenuation due to absorption in the AR coatings or any layer in the front
        # that is not part of the junction
        if type(layer_object) is Layer:
            all_layers.append(layer_object)
            n_layers_junction.append(1)

        # For each junction, and layer within the junction, we get the absorption
        # coefficient and the layer width.
        elif type(layer_object) in [TunnelJunction, Junction]:
            n_layers_junction.append(len(layer_object))
            for i, layer in enumerate(layer_object):
                all_layers.append(layer)

    # With all the information, we create the optical stack
    stack = all_layers
    position = position * 1e9
    substrate = solar_cell.substrate
    incidence = solar_cell.incidence

    print("Calculating RAT...")
    RAT = calculate_rat_rcwa(
        stack,
        size,
        orders,
        wavelength * 1e9,
        incidence,
        substrate,
        theta=theta,
        phi=phi,
        pol=pol,
        parallel=parallel,
        user_options=rcwa_options,
    )

    print("Calculating absorption profile...")
    out = calculate_absorption_profile_rcwa(
        stack,
        size,
        orders,
        wavelength * 1e9,
        RAT["A_pol"],
        dist=position,
        theta=theta,
        phi=phi,
        pol=pol,
        incidence=incidence,
        substrate=substrate,
        parallel=parallel,
        n_jobs=n_jobs,
        user_options=rcwa_options,
    )
    # With all this information, we are ready to calculate the differential absorption
    # function
    diff_absorption, _ = calculate_absorption_rcwa(out, initial)

    layer = 0
    A_per_layer = np.array(RAT["A_per_layer"].T)

    # Each building block (layer or junction) needs to have access to the absorbed light
    # in its region. We update each object with that information.
    for j in range(len(solar_cell)):
        solar_cell[j].diff_absorption = diff_absorption
        solar_cell[j].absorbed = types.MethodType(absorbed, solar_cell[j])

        solar_cell[j].layer_absorption = initial * np.sum(
            A_per_layer[layer : (layer + n_layers_junction[j])], axis=0
        )
        layer = layer + n_layers_junction[j]

    solar_cell.reflected = RAT["R"] * initial
    solar_cell.absorbed = sum(
        [solar_cell[x].layer_absorption for x in np.arange(len(solar_cell))]
    )
    solar_cell.transmitted = initial - solar_cell.reflected - solar_cell.absorbed


def absorbed(self, z):
    out = self.diff_absorption(self.offset + z) * (z < self.width)
    return out.T


def calculate_absorption_rcwa(tmm_out, initial=1):
    all_z = tmm_out["position"] * 1e-9
    all_abs = initial * tmm_out["absorption"] / 1e-9

    def diff_absorption(z):
        idx = all_z.searchsorted(z)
        idx = np.where(idx <= len(all_z) - 2, idx, len(all_z) - 2)
        try:
            z1 = all_z[idx]
            z2 = all_z[idx + 1]

            f = (z - z1) / (z2 - z1)

            out = f * all_abs[:, idx] + (1 - f) * all_abs[:, idx + 1]

        except IndexError:
            out = all_abs[:, idx]

        return out

    all_absorbed = np.trapz(diff_absorption(all_z), all_z)

    return diff_absorption, all_absorbed
