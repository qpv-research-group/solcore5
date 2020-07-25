from __future__ import annotations

from typing import NamedTuple, Optional, Dict, Sequence, Callable, Union
from dataclasses import dataclass
from warnings import warn

import xarray as xr
import numpy as np

from ..junction_base import JunctionBase
from ..structure import Structure
from ..constants import kb, q, vacuum_permittivity as e0
from ..analytic_solar_cells.depletion_approximation import get_Jsrh


NK_DATA_SIGNATURE = Callable[[np.ndarray], xr.DataArray]


class AbruptHomojunctionData(NamedTuple):
    """Parameters of an abrupt homojunction for the depletion approximation.

    The subindex always refer to the layer, so 'xp' is the width of the P region,
    'ln' represents the diffusion length of minority carriers in the N region,
    i.e. holes, etc.

    polarity (str): Polarity of the junction, either 'pn', 'np', 'pin' or 'nip'.
    t (float): Junction temperature, in K.
    xn (float): Width of the N region, in m.
    xp (float): Width of the P region, in m.
    xi (float): Width of the I region, in m.
    sp (float): Surface recombination velocity at the P side, in m/s.
    sn (float): Surface recombination velocity at the N side, in m/s.
    lp (float): Diffusion length of electrons in the P region, in m.
    ln (float): Diffusion length of holes in the N region, in m.
    dp (float): Diffusion coefficient on the P region, in m^2/s.
    dn (float): Diffusion coefficient on the N region, in m^2/s.
    Nd (float): Doping of the N region, in m^-3.
    Na (float): Doping of the P region, in m^-3.
    ni (float): Intrinsic carrier concentration, in m^-3.
    es (float): Relative permittivity.
    bandgap (float): Bandgap of the material, in eV.
    rs (float): Series resistance, in Ohm.
    rsh (float): Shunt resistance, in Ohm.
    """

    polarity: str = "pn"
    t: float = 297.0
    xn: float = 2e-7
    xp: float = 3e-6
    xi: float = 0.0
    sn: float = 0.0
    sp: float = 0.0
    ln: float = 5e-6
    lp: float = 5e-7
    dn: float = 8.7e-5
    dp: float = 1.2e-3
    Nd: float = 3e24
    Na: float = 1e23
    ni: float = 2.25e12
    es: float = 9.0
    bandgap: float = 1.4
    rs: float = 0.0
    rsh: float = 1e14


@dataclass(frozen=True)
class AbruptHomojunction(JunctionBase):
    """Junction class for abrupt homojunctions using the depletion approximation.

    The junction is assumed to be abrupt, with clearly defined n, p and, optionally,
    intrinsic regions.

    Attributes:
        data (TwoDiodeData): Object storing the parameters of a 2 diode equation.
        solver (str): Name of the IV solver to use. Options are 'builtin' and 'spice'.
        params (Dict): Other parameters with physical information, possibly used
            during the construction of the TwoDiodeData object.
        options (Dict): Options to pass to the iv calculator.
    """

    data: AbruptHomojunctionData = AbruptHomojunctionData()
    nk_data: Union[xr.DataArray, NK_DATA_SIGNATURE, None] = None
    layer_widths: Optional[xr.DataArray] = None
    structure: Optional[Structure] = None

    @classmethod
    def from_structure(
        cls,
        structure: Structure,
        sn: float = 0.0,
        sp: float = 0.0,
        rs: float = 0.0,
        rsh: float = 1e14,
    ):
        """Creates a DepletionAbruptHomojunction from layers and materials.

        This method scans the structure, identifies the different regions, their
        polarity and extract the values of the relevant parameters. Only adjacent n-i-p
        regions will be considered: if there is, for example, an n++ window layer on
        top of the emitter, this will be ignored. The same applies to any back surface
        field layers. Any expected impact on the electrical properties of those layers
        should be included bundled in the surface recombination velocities.

        While not used in the electrical calculations, those extra layers are indeed
        considered in the calculation of the optical properties of the junction.

        Args:
            structure (Structure): A Structure object containing the layers of the
                junction with the different materials.
            sn (float): Surface recombination velocity at the N side, in m/s.
            sp (float): Surface recombination velocity at the P side, in m/s.
            rs (float): Series resistance, in Ohm.
            rsh (float): Shunt resistance, in Ohm.

        Raises:
            ValueError: If it cannot identify an n-(i)-p sequence of materials.
            ValueError: If the n-(i)-p materials are different.
            ValueError: If the n-(i)-p have different temperatures.
            AttributeError: If any of the materials lack of an essential attribute.

        Returns:
            New instance of AbruptHomojunction
        """
        idx, polarity = cls._find_polarity(structure)

        mats = tuple(structure[i].material.name for i in idx)
        comp = tuple(structure[i].material.composition for i in idx)
        if len(set(mats)) != 1 or len(set(comp)) != 1:
            raise ValueError(
                "The materials are different: Not an homojunction! "
                f"Found materials and compositions are: {list(zip(mats, comp))}"
            )

        temps = tuple(structure[i].material.T for i in idx)
        if len(set(temps)) != 1:
            raise ValueError(
                "The materials have different temperatures. "
                f"Found temperatures are: {temps}"
            )

        layers = structure[slice(idx[0], idx[-1] + 1, 1)]
        pregion = layers[0] if polarity[0] == "p" else layers[-1]
        nregion = layers[-1] if polarity[0] == "p" else layers[0]

        d: Dict = dict(polarity=polarity, t=temps[0], sn=sn, sp=sp, rs=rs, rsh=rsh)
        d["xp"] = pregion.width
        d["xi"] = layers[1].width if len(layers) == 3 else 0.0
        d["xn"] = nregion.width
        d["lp"] = pregion.material.electron_diffusion_length
        d["ln"] = nregion.material.hole_diffusion_length
        d["dp"] = pregion.material.electron_mobility * kb * d["t"] / q
        d["dn"] = nregion.material.hole_mobility * kb * d["t"] / q
        d["Nd"] = nregion.material.Nd
        d["Na"] = pregion.material.Na
        d["ni"] = pregion.material.ni
        d["es"] = pregion.material.relative_permitivity
        d["band_gap"] = pregion.material.band_gap

        data = AbruptHomojunctionData(**d)
        layer_widths = xr.DataArray([layer.width for layer in structure], dims="layer")

        def nk_data(wl: np.ndarray) -> xr.DataArray:
            return xr.DataArray(
                [
                    layer.material.n(wl) + 1.0j * layer.material.p(wl)
                    for layer in structure
                ],
                dims=("layer", "wavelength"),
                coords={"wavelength": wl},
            )

        return cls(
            data=data, nk_data=nk_data, layer_widths=layer_widths, structure=structure
        )

    @staticmethod
    def _find_polarity(structure: Structure) -> (Sequence, str):
        """Finds the polarity of a Structure of layers with materials.

        Args:
            structure (Structure): A Structure object containing the layers of the
                junction with the different materials.

        Raises:
            ValueError: If it cannot identify an n-(i)-p sequence of materials.

        Returns:
            (Sequence, str) Tuple with the indices of the n-i-p layers and the polarity.
        """
        pol = []
        idx = []
        for i, layer in enumerate(structure):
            acceptor = getattr(layer.material, "Na", 0)
            donor = getattr(layer.material, "Nd", 0)
            if acceptor > donor:
                doping = "p"
            elif acceptor < donor:
                doping = "n"
            else:
                doping = "i"

            if pol and pol[-1] == doping:
                idx[-1] = i
            else:
                pol.append(doping)
                idx.append(i)

            if "".join(pol) in ("pn", "np", "pin", "nip"):
                polarity = "".join(pol)
                break
        else:
            raise ValueError("A valid n-(i)-p sequence of layers could not be found.")

        return idx, polarity

    def widths(self) -> Optional[xr.DataArray]:
        """Provides the widths of all layers the junction, in m.

        Returns:
            An xr.DataArray with the widths or None. The only coordinate must be
            'layer'.
        """
        return self.layer_widths if self.layer_widths else None

    def nk(self, wavelength: np.ndarray) -> Optional[xr.DataArray]:
        """Provides the complex refractive index of all layers of the junction.

        Args:
            wavelength (np.ndarray): Array with the wavelengths in meters.

        Returns:
            A xr.DataArray with the complex refractive index as a function of two
            coordinates, 'wavelength' and 'layer'.
        """
        if isinstance(self.nk_data, xr.DataArray):
            return self.nk_data.interp({"wavelength": wavelength})
        elif isinstance(self.nk_data, Callable):
            return self.nk_data(wavelength)

    def solve_iv(
        self,
        voltage: np.ndarray,
        absorption: Optional[xr.DataArray] = None,
        source: Optional[xr.DataArray] = None,
    ) -> xr.Dataset:
        """Calculates the IV curve of the junction.

        If absorption is provided, then light_source must also be provided and the
        light IV curve should be calculated instead. In this case, parameters like
        Voc, Isc, fill factor, etc. are also calculated.

        None of the currents provided are affected by any series resistance.

        Args:
            voltage (np.ndarray): Array of voltages at which to calculate the IV curve.
            absorption (xr.DataArray, optional): Array with the fraction of absorbed
                light as a function of 'wavelength' and 'position'.
            source (xr.DataArray, optional): Light source to use providing the number
                of photons as a junction of 'wavelength'.

        Returns:
            A xr.Dataset with the output of the calculation. It contains the 'current'
            DataArray giving the total current in A as a function of the input 'voltage'
            as well as:
            - The dark currents from the quasi neutral regions and the depletion
                region.
            - The shunt current.

            If light IV is calculated:
            - The short circuit currents collected from the three regions.
            - The curve parameters (Voc, Isc, FF, Vmpp, Impp and Pmpp)
            """
        vbi = built_in_voltage(self.data.t, self.data.Nd, self.data.Na, self.data.ni)
        v = np.clip(voltage, a_min=np.min(voltage), a_max=vbi - 0.001)
        wn = depletion_widths(
            v, vbi, self.data.Nd, self.data.Na, self.data.es, self.data.xi
        )
        wp = depletion_widths(
            v, vbi, self.data.Na, self.data.Nd, self.data.es, self.data.xi
        )

        # Dark current from the N and P quasi neutral regions.
        jn_dark = xr.DataArray(
            current_quasi_neutral(
                v,
                self.data.xn,
                wn,
                self.data.ln,
                self.data.sn,
                self.data.dn,
                self.data.ni ** 2 / self.data.Nd,
                self.data.t,
            ),
            dims="voltage",
            coords={"voltage": voltage},
        )
        jp_dark = xr.DataArray(
            current_quasi_neutral(
                v,
                self.data.xp,
                wp,
                self.data.lp,
                self.data.sp,
                self.data.dp,
                self.data.ni ** 2 / self.data.Na,
                self.data.t,
            ),
            dims="voltage",
            coords={"voltage": voltage},
        )

        # Dark current from the depleted region
        j_depletion = xr.DataArray(
            get_Jsrh(
                self.data.ni,
                v,
                vbi,
                self.data.lp ** 2 / self.data.dp,
                self.data.ln ** 2 / self.data.dn,
                wn + wp + self.data.xi,
                kb * self.data.t,
            ),
            dims="voltage",
            coords={"voltage": voltage},
        )

        # Shunt current
        j_shunt = xr.DataArray(
            v / self.data.rsh, dims="voltage", coords={"voltage": voltage},
        )

        return xr.Dataset(
            {
                "current": jn_dark + jp_dark + j_depletion + j_shunt,
                "jn_dark": jn_dark,
                "jp_dark": jp_dark,
                "j_depletion": j_depletion,
                "j_shunt": j_shunt,
            }
        )

    def solve_qe(
        self, absorption: xr.DataArray, source: Optional[xr.DataArray] = None
    ) -> xr.Dataset:
        pass

    def solve_equilibrium(self):
        pass

    def solve_short_circuit(
        self, absorption: xr.DataArray, source: xr.DataArray
    ) -> xr.Dataset:
        raise NotImplementedError


def built_in_voltage(t: float, Nd: float, Na: float, ni: float) -> float:
    """Calculates the built-in voltage of the junction.

    Args:
        t (float): Junction temperature, in K.
        Nd (float): Doping of the N region, in m^-3.
        Na (float): Doping of the P region, in m^-3.
        ni (float): Intrinsic carrier concentration, in m^-3.

    Returns:
        The built-in voltage, in V.
    """
    return (kb * t / q) * np.log(Nd * Na / ni ** 2)


def depletion_widths(
    v: np.ndarray, vbi: float, Nthis: float, Nother: float, es: float, xi: float = 0.0,
) -> np.ndarray:
    """Calculates the depletion width of a region.

    Args:
        v (np.ndarray): Array of voltages.
        vbi (float): Built-in voltage.
        Nthis (float): Doping of this region, in m^-3.
        Nother (float): Doping of the opposite region, in m^-3.
        es (float): Relative permittivity.
        xi (float): Width of the I region, in m.

    Returns:
        The depletion width of the region
    """
    return -xi + np.sqrt(
        xi ** 2
        + 2.0 * es * e0 * (vbi - v) / q * Nother / Nthis * (1 / (Nother + Nthis))
    )


def current_quasi_neutral(
    v: np.ndarray,
    x: float,
    w: np.ndarray,
    ll: float,
    s: float,
    d: float,
    n_m: float,
    t: float,
):
    """Calculate the current of the quasi-neutral regions.

    Implements equations 6.34/6.39 of:
    J. Nelson, “The Physics of Solar Cells”, Imperial College Press (2003).

    Args:
        v (np.ndarray): Applied voltage, in v.
        x (float): Total width of the region, in m.
        w (np.ndarray): Depleted width of the region as a function of voltage, in m.
        ll (float): Minority carriers diffusion length, in m.
        s (float): Surface recombination velocity in m/2.
        d (float): Diffusion coefficient, in m^2/s
        n_m (float): Concentration of minority carriers.
        t (float): Junction temperature, in K.

    Raises:
        ValueError: If the quasi-neutral region width is reduced to less than zero.

    Returns:
        The current of the quasi-neutral region at the applied voltages.
    """
    harg = (x - w) / ll
    if (harg < 0).any():
        warn(
            "Negative value found for the quasi-neutral region width. "
            "Clipping to zero.",
            UserWarning,
        )
        harg[harg < 0] = 0

    sinh_harg = np.sinh(harg)
    cosh_harg = np.cosh(harg)
    lsod = (ll * s) / d

    return (
        (q * d * n_m / ll)
        * (np.exp(q * v / kb / t) - 1)
        * ((lsod * cosh_harg + sinh_harg) / (lsod * sinh_harg + cosh_harg))
    )
