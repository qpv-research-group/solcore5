from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Dict, NamedTuple, Callable
from copy import deepcopy
import numpy as np
import xarray as xr
from scipy.optimize import root

from ..junction_base import JunctionBase
from ...constants import q, kb, hbar, c


class TwoDiodeData(NamedTuple):
    """Object storing the parameters of a 2 diode equation.

    t (float): Junction temperature.
    j01 (float): Saturation current 1.
    j02 (float): Saturation current 2.
    n1 (float): Ideality factor 1.
    n2 (float): Ideality factor 2.
    rs (float): Series resistance.
    rsh (float): Shunt resistance.
    js (float, optional): Short circuit current. If provided, it will be used in the IV
        calculations instead of calculating it from the absorption and the spectrum.
    """

    t: float = 297
    j01: float = 1e-6
    j02: float = 0.0
    n1: float = 1.0
    n2: float = 2.0
    rs: float = 0.0
    rsh: float = 1e14
    jsc: Optional[float] = None


IV_SOLVERS: Dict[str, Callable[[np.ndarray, TwoDiodeData, Dict], np.ndarray]] = {}
"""Registry of IV solvers."""


def register_iv_solver(fun: Optional[Callable] = None, name: Optional[str] = None):

    if fun is None:
        return lambda x: register_iv_solver(x, name=name)

    name = name if name is not None else fun.__name__
    IV_SOLVERS[name] = fun

    return fun


@dataclass(frozen=True)
class TwoDiodeJunction(JunctionBase):
    """Junction class for the two diode model.

    Attributes:
        data (TwoDiodeData): Object storing the parameters of a 2 diode equation.
        solver (str): Name of the IV solver to use. Options are 'builtin' and 'spice'.
        params (Dict): Other parameters with physical information, possibly used
            during the construction of the TwoDiodeData object.
        options (Dict): Options to pass to the iv calculator.
    """

    data: TwoDiodeData = TwoDiodeData()
    solver: str = "builtin"
    params: Optional[Dict] = field(default_factory=dict)
    options: Optional[Dict] = field(default_factory=dict)

    def __post_init__(self):
        if self.solver not in IV_SOLVERS:
            raise KeyError(
                f"Invalid IV solver: '{self.solver}'. "
                f"Valid options are: {list(IV_SOLVERS.keys())}"
            )
        if self.params is None:
            object.__setattr__(self, "params", {})
        if self.options is None:
            object.__setattr__(self, "options", {})

    @classmethod
    def from_reference_temperature(
        cls,
        band_gap: float,
        t: float,
        data: TwoDiodeData = TwoDiodeData(),
        solver: str = "builtin",
        params: Optional[Dict] = None,
        options: Optional[Dict] = None,
    ):
        """Creates a TwoDiodeJunction from TwoDiodeData at different temperature.

        We want the junction to be defined at a certain temperature t, but we know
        the parameters at a different temperature, contained in the data object.
        Assuming the temperatures are not too different, we can estimate a new set of
        parameters using the bandgap of the junction.

        Args:
            band_gap: Bandgap associated to the junction, in J.
            t: Temperature of interest.
            data (TwoDiodeData): Object storing the parameters of a 2 diode equation at
                a reference temperature.
            solver (str): Name of the IV solver to use. Options are 'builtin' and
                'spice'.
            params (Dict): Other parameters with physical information.
            options (Dict): Options to pass to the iv calculator - used if series
                resistance is not None.

        Returns:
            New instance of a TwoDiodeJunction
        """
        j01 = (
            data.j01
            * (t / data.t) ** 3
            * np.exp(-band_gap / (data.n1 * kb) * (1 / t - 1 / data.t))
        )
        j02 = (
            data.j02
            * (t / data.t) ** (5.0 / 3.0)
            * np.exp(-band_gap / (data.n2 * kb) * (1 / t - 1 / data.t))
        )

        parameters = deepcopy(params) if params is not None else {}
        parameters.update({"band_gap": band_gap, "t_ref": data.t})
        options = options if options else {}
        return cls(data._replace(t=t, j01=j01, j02=j02), solver, parameters, options)

    @classmethod
    def from_bandgap_absorption(
        cls,
        band_gap: float,
        n: float,
        data: TwoDiodeData = TwoDiodeData(),
        solver: str = "builtin",
        params: Optional[Dict] = None,
        options: Optional[Dict] = None,
    ):
        """Creates a TwoDiodeJunction calculating j01 out of the bandgap absorption.

        Calculate the reverse saturation current j01, assumed radiative, considering an
        absorption equal to 1 above the bandgap. Light trapping is included by
        considering the refractive index of the material:

        .. math:: J_{01} = \\frac {q n^2 k_b T} {2 \\pi ^2 c^2 \\hbar ^3}
                    e^{\\frac{-E_g}{k_b T}} (E_g^2 + 2 k_b T E_g + 2 k_b^2 T^2)

        Args:
            band_gap (float): Bandgap associated to the junction, in J.
            n (float): Refractive index of the junction.
            data (TwoDiodeData): Object storing the parameters of a 2 diode equation.
            solver (str): Name of the IV solver to use. Options are 'builtin' and
                'spice'.
            params (Dict): Other parameters with physical information.
            options (Dict): Options to pass to the iv calculator.

        Returns:
            New instance of a TwoDiodeJunction
        """
        term1 = 2 * np.pi * n ** 2 * q / (4 * np.pi ** 3 * hbar ** 3 * c ** 2)
        term2 = kb * data.t * np.exp(-band_gap / (kb * data.t))
        term3 = (
            band_gap ** 2 + (2 * band_gap * kb * data.t) + (2 * kb ** 2 * data.t ** 2)
        )

        parameters = deepcopy(params) if params is not None else {}
        parameters.update({"band_gap": band_gap})
        return cls(
            data._replace(j01=term1 * term2 * term3), solver, parameters, options
        )

    @classmethod
    def from_voc(
        cls,
        voc: float,
        data: TwoDiodeData = TwoDiodeData(),
        solver: str = "builtin",
        params: Optional[Dict] = None,
        options: Optional[Dict] = None,
    ):
        """Creates a TwoDiodeJunction calculating j02 based on the j01, jsc and the Voc.

        It is just the result of solving the 2-diode equation for j02.

        Args:
            voc (float): Open circuit voltage, in V.
            data (TwoDiodeData): Object storing the parameters of a 2 diode equation.
            solver (str): Name of the IV solver to use. Options are 'builtin' and
                'spice'.
            params (Dict): Other parameters with physical information.
            options (Dict): Options to pass to the iv calculator.

        Returns:
            New instance of a TwoDiodeJunction
        """
        if data.jsc is None:
            msg = "'jsc' cannot be None when creating a TwoDiodeJunction from Voc."
            raise ValueError(msg)

        term1 = (
            data.jsc
            - data.j01 * (np.exp(q * voc / (data.n1 * kb * data.t)) - 1)
            - voc / data.rsh
        )
        term2 = np.exp(q * voc / (data.n2 * kb * data.t)) - 1

        parameters = deepcopy(params) if params is not None else {}
        parameters.update({"Voc": voc})
        return cls(data._replace(j02=term1 / term2), solver, parameters, options)

    @classmethod
    def from_radiative_efficiency(
        cls,
        v_ref: float,
        reff: float,
        data: TwoDiodeData = TwoDiodeData(),
        solver: str = "builtin",
        params: Optional[Dict] = None,
        options: Optional[Dict] = None,
    ):
        """Creates 2D junction calculating j02 based on j01 and radiative efficiency.

        Args:
            v_ref (float): Reference voltage at which the radiative efficiency is known.
            reff (float): Radiative efficiency.
            data (TwoDiodeData): Object storing the parameters of a 2 diode equation.
            solver (str): Name of the IV solver to use. Options are 'builtin' and
                'spice'.
            params (Dict): Other parameters with physical information.
            options (Dict): Options to pass to the iv calculator.

        Returns:
            New instance of a TwoDiodeJunction
        """

        term1 = data.j01 * (np.exp(q * v_ref / (data.n1 * kb * data.t)) - 1)
        term2 = 1 / reff - 1
        term3 = np.exp(q * v_ref / (data.n2 * kb * data.t)) - 1

        parameters = deepcopy(params) if params is not None else {}
        parameters.update({"Vref": v_ref, "Reff": reff})
        return cls(
            data._replace(j02=(term1 * term2 + v_ref / data.rsh) / term3),
            solver,
            parameters,
            options,
        )

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

        Args:
            voltage (np.ndarray): Array of voltages at which to calculate the IV curve.
            absorption (xr.DataArray, optional): Array with the fraction of absorbed
                light as a function of 'wavelength' and 'position'.
            source (xr.DataArray, optional): Light source to use providing the number
                of photons as a junction of 'wavelength'.

        Returns:
            A xr.Dataset with the output of the calculation. Contains a 'current'
            DataArray giving the current in amps as a function of the input 'voltage'.
            If light IV is calculated, the curve parameters (Voc, Isc, FF, Vmpp, Impp,
            Pmpp and eta) are provided as attributes of the Dataset.
        """
        jsc = (
            self.data.jsc
            if self.data.jsc is not None
            else self._get_jsc(absorption, source)
        )

        i = IV_SOLVERS[self.solver](
            voltage, self.data._replace(jsc=jsc), **self.options
        )

        current = xr.DataArray(i, dims=["voltage"], coords={"voltage": voltage})
        parameters = {} if jsc == 0.0 else self.iv_parameters(voltage, current.values)
        return xr.Dataset({"current": current}, attrs=parameters)

    def solve_qe(
        self, absorption: xr.DataArray, source: Optional[xr.DataArray] = None
    ) -> xr.Dataset:
        """Calculates the external and internal quantum efficiency of the junction.

        Args:
            absorption (xr.DataArray, optional): Array with the fraction of absorbed
                light as a function of 'wavelength' and 'position'.
            source (xr.DataArray, optional): Light source to use providing the number
                of photons as a junction of 'wavelength'.

        Returns:
            A xr.Dataset with the 'eqe' and 'iqe' DataArrays as a function of
            'wavelength'.
        """
        eqe = absorption.integrate("position")
        iqe = xr.ones_like(eqe)

        return xr.Dataset({"eqe": eqe, "iqe": iqe})

    def solve_equilibrium(self):
        raise NotImplementedError

    def solve_short_circuit(
        self, absorption: xr.DataArray, light_source: xr.DataArray
    ) -> xr.Dataset:
        raise NotImplementedError

    def _get_jsc(
        self,
        absorption: Optional[xr.DataArray] = None,
        source: Optional[xr.DataArray] = None,
    ) -> float:
        """Calculates the short circuit current out of the absorption.

        Args:
            absorption (xr.DataArray, optional): Array with the fraction of absorbed
                light as a function of 'wavelength' and 'position'.
            source (xr.DataArray, optional): Light source to use providing the number
                of photons as a junction of 'wavelength'.

        Returns:
            The short circuit current.
        """
        if absorption is None or source is None:
            return 0.0

        sp = source.interp({"wavelength": absorption.wavelength})
        return q * (self.solve_qe(absorption).eqe * sp).integrate("wavelength")


def iv_no_rs(v: np.ndarray, data: TwoDiodeData = TwoDiodeData()) -> np.ndarray:
    """Calculates the current at the chosen voltages using the given parameters.

    Args:
        v (np.ndarray): Voltages at which to calculate the currents.
        data (TwoDiodeData): Object storing the parameters of a 2 diode equation.

    Returns:
        Numpy array of the same length as the voltages with the currents.
    """
    return (
        data.j01 * (np.exp(q * v / (data.n1 * kb * data.t)) - 1)
        + data.j02 * (np.exp(q * v / (data.n2 * kb * data.t)) - 1)
        + v / data.rsh
        - data.jsc
    )


@register_iv_solver(name="builtin")
def iv2diode_default(
    v: np.ndarray, data: TwoDiodeData = TwoDiodeData(), **kwargs,
) -> np.ndarray:
    """Calculates the current using the 2-diodes equation.

    If series resistance is zero, it just return the result of replacing all the
    parameters in the 2 diode equation. Otherwise, scipy.optimize.root is used to
    numerically solve the transcendental equation.

    Args:
        v (np.ndarray): Voltages at which to calculate the currents.
        data (TwoDiodeData): Object storing the parameters of a 2 diode equation.
        kwargs: Options to be passed to the root finding algorithm.

    Returns:
        Numpy array of the same length as the voltages with the currents.
    """
    current = iv_no_rs(v, data)
    if data.rs == 0.0:
        return current

    def fun(j):
        return j - iv_no_rs(v - data.rs * j, data)

    jguess = np.clip(current, np.min(v) / data.rs, np.max(v) / data.rs)
    sol = root(fun, jguess, **kwargs)
    if not sol.success:
        msg = f"The IV calculation failed to converge. Solver info:\n{sol}"
        raise RuntimeError(msg)

    return sol.x


@register_iv_solver(name="spice")
def iv2diode_spice(
    v: np.ndarray, data: TwoDiodeData = TwoDiodeData(), **kwargs,
) -> np.ndarray:
    """Solves the 2-diodes equation in SPICE.

    This function uses a SPICE engine to solve the 2-diode equation. The SPICE engine
    needs to have been configured separately.

    Args:
        v (np.ndarray): Voltages at which to calculate the currents.
        data (TwoDiodeData): Object storing the parameters of a 2 diode equation.
        kwargs: Options to be passed to the SPICE solver.

    Returns:
        Numpy array of the same length as the voltages with the currents.
    """
    raise NotImplementedError("The SPICE solver is not implemented, yet.")
