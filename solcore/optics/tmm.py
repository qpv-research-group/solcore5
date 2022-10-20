import types
from typing import List, Optional

import numpy as np
from numpy.typing import NDArray

from ..absorption_calculator import (
    OptiStack,
    calculate_absorption_profile,
    calculate_rat,
)
from ..registries import register_optics
from ..solar_cell import SolarCell
from ..structure import Junction, Layer, TunnelJunction


@register_optics(name="TMM")
def solve_tmm(
    solar_cell: SolarCell,
    wavelength: NDArray,
    position: NDArray,
    BL_correction: bool = True,
    theta: float = 0.0,
    pol: str = "u",
    zero_threshold: float = 1e-5,
    no_back_reflection: bool = True,
    coherency_list: Optional[List[str]] = None,
    **kwargs
) -> None:
    """Calculates the RAT of a solar cell object using the transfer matrix method.

    Internally, it creates an OptiStack and then it calculates the optical properties of
    the whole structure. A substrate can be specified in the SolarCell object, which is
    treated as a semi-infinite transmission medium. Shading can also be specified (as a
    fraction).

    A coherency_list option can be provided:

    Args:
        solar_cell: A solar_cell object
        wavelength: Array of wavelegth at which the optics are calculated.
        position: Array of positions in the z direction to calculate the absorption vs
            depth.
        BL_correction: If is set to True, thick layers (thickness > 10*maximum
            wavelength) are treated incoherently using the Beer-Lambert law, to avoid
            the calculation of unphysical interference oscillations in the R/A/T
            spectra.
        theta: the polar incidence angle, in degrees, with 0 degrees being normal
            incidence.
        pol: the polarization of the light ('s', 'p' or 'u')
        zero_threshold: when the fraction of incident light absorbed in a layer is
            less than this value, the absorption profile is completely set to zero for
            both coherent and incoherent calculations. This is applied on a
            wavelength-by-wavelength basis and is intended to prevent errors where
            integrating a weak absorption profile in a layer over many points leads to
            calculated EQE > total absorption in that layer.
        no_back_reflection: Sets whether reflections from the back surface are
            suppressed (if set to True, the default), or taken into account (if set to
            False).
        coherency_list: If present, this should have the same number of elements than
            number of layers (if a Junction contains multiple Layers, each should have
            its own entry in the coherency list). Each element is either 'c' for
            coherent treatment of that layer or 'i' for incoherent treatment.

    Return:
        None
    """
    # We include the shadowing losses
    initial = (1 - solar_cell.shading) if hasattr(solar_cell, "shading") else 1

    # Now we calculate the absorbed and transmitted light. We first get all the relevant
    # parameters from the objects
    all_layers = []
    widths = []
    n_layers_junction = []

    for j, layer_object in enumerate(solar_cell):

        # Attenuation due to absorption in the AR coatings or any layer in the front
        # that is not part of the junction
        if type(layer_object) is Layer:
            all_layers.append(layer_object)
            widths.append(layer_object.width)
            n_layers_junction.append(1)

        # For each junction, and layer within the junction, we get the absorption
        # coefficient and the layer width.
        elif type(layer_object) in [TunnelJunction, Junction]:
            n_layers_junction.append(len(layer_object))
            for i, layer in enumerate(layer_object):
                all_layers.append(layer)
                widths.append(layer.width)

    # With all the information, we create the optical stack
    full_stack = OptiStack(
        all_layers,
        no_back_reflection=no_back_reflection,
        substrate=solar_cell.substrate,
        incidence=solar_cell.incidence,
    )

    if coherency_list is not None:
        coherent = False
        if len(coherency_list) != full_stack.num_layers:
            raise ValueError(
                "Error: The coherency list must have as many elements (now {}) as the "
                "number of layers (now {}).".format(
                    len(coherency_list), full_stack.num_layers
                )
            )
    else:
        coherent = True

    # assume it's safe to ignore interference effects
    if BL_correction and any(widths > 10 * np.max(wavelength)):
        make_incoherent = np.where(np.array(widths) > 10 * np.max(wavelength))[0]
        print("Treating layer(s) " + str(make_incoherent).strip("[]") + " incoherently")
        if not coherency_list:
            coherency_list_ = np.array(len(all_layers) * ["c"])
            coherent = False
        else:
            coherency_list_ = np.array(coherency_list)
        coherency_list_[make_incoherent] = "i"
        coherency_list = coherency_list_.tolist()

    position = position * 1e9
    profile_position = position[position < sum(full_stack.widths)]

    print("Calculating RAT...")
    RAT = calculate_rat(
        full_stack,
        wavelength * 1e9,
        angle=theta,
        coherent=coherent,
        coherency_list=coherency_list,
        no_back_reflection=no_back_reflection,
        pol=pol,
    )

    print("Calculating absorption profile...")
    out = calculate_absorption_profile(
        full_stack,
        wavelength * 1e9,
        dist=profile_position,
        angle=theta,
        no_back_reflection=no_back_reflection,
        pol=pol,
        coherent=coherent,
        coherency_list=coherency_list,
        zero_threshold=zero_threshold,
        RAT_out=RAT,
    )

    # With all this information, we are ready to calculate the differential absorption
    # function
    diff_absorption, all_absorbed = calculate_absorption_tmm(out, initial)

    # Each building block (layer or junction) needs to have access to the absorbed light
    # in its region. We update each object with that information.
    # first entry is R, last entry is T
    layer = 0
    A_per_layer = np.array(RAT["A_per_layer"][1:-1])

    for j in range(len(solar_cell)):

        solar_cell[j].layer_absorption = initial * np.sum(
            A_per_layer[layer : (layer + n_layers_junction[j])], axis=0
        )

        solar_cell[j].diff_absorption = diff_absorption
        solar_cell[j].absorbed = types.MethodType(absorbed, solar_cell[j])

        layer = layer + n_layers_junction[j]

    solar_cell.reflected = RAT["R"] * initial
    solar_cell.absorbed = sum(
        [solar_cell[x].layer_absorption for x in np.arange(len(solar_cell))]
    )
    solar_cell.transmitted = initial - solar_cell.reflected - solar_cell.absorbed


def absorbed(self, z):
    out = self.diff_absorption(self.offset + z) * (z < self.width)

    return out.T


def calculate_absorption_tmm(tmm_out, initial=1):
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
