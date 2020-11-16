"""

"""
from functools import partial
from inspect import Parameter as insParam
from inspect import signature
from itertools import chain
from logging import getLogger
from typing import Callable, Dict, List, Optional, Tuple, Union

from pint import Quantity
from solcore.constants import electron_mass, vacuum_permittivity
from solcore.parameter import (
    Parameter,
    ParameterError,
    ParameterSource,
    ParameterSourceBase
)


class CalculableParameters(ParameterSourceBase):

    name: Union[str, List[str]] = "Calculable"
    _priority: Union[int, Dict[str, int]] = -10
    __instance = None

    def __new__(cls, *args, **kwargs):
        if cls.__instance is None:
            cls.__instance = ParameterSourceBase.__new__(cls)
            cls.__instance._params = {}
            cls.__instance._descriptions = {}
            cls.__instance._warned = False
        return cls.__instance

    def __getitem__(self, parameter: str) -> Callable:
        """Dictionary-like access to the calculable source.

        Args:
            parameter (str): The parameter of interest.

        Raises:
            ParameterError: If the parameter does not exists in the source

        Returns:
            The callable for this parameter
        """
        if parameter not in self._params:
            raise ParameterError(
                f"Parameter '{parameter}' not in '{self.name}' parameters source"
            )
        return self._params[parameter]

    def __getattr__(self, parameter: str) -> Callable:
        """Attribute-like access to the calculable source.

        Args:
            parameter (str): The parameter of interest.

        Raises:
            ParameterError: If the parameter does not exists in the source

        Returns:
            The callable for this parameter
        """
        return self[parameter]

    @classmethod
    def load_source(cls, source_name: str = "") -> ParameterSource:
        return cls.__instance

    @classmethod
    def register_calculable(
        cls, function: Optional[Callable] = None, description: Optional[str] = ""
    ):
        """Register a calculable with the calculable parameter source.

        Args:
            function: Calculable to be registered
            description: Description of the calculable parameter

        Returns:
            The same input
        """
        if function is None:
            return partial(cls.register_calculable, description=description)

        cls()._params[function.__name__] = function
        cls()._descriptions[function.__name__] = description
        return function

    @property
    def materials(self) -> Tuple[str, ...]:
        """Materials this source provides parameters for.

        Returns:
            An empty tuple
        """
        warning_msg = (
            "The parameters in the the CalculableParameter source are available for "
            "any material that can provide the required input arguments to the "
            "calculable. Use  the method 'list_arguments' to get a list of arguments "
            "required for a given calculable parameter."
        )
        if not self._warned:
            getLogger().warning(warning_msg)
            self._warned = True

        return ()

    def parameters(self, material: str = "") -> Tuple[str, ...]:
        """Parameters available in this source for the requested material.

        Args:
            material (str): The material whose parameters are of interests. It is
                ignored in this source.

        Returns:
            A tuple with the parameters for this material that this source provides.
        """
        return tuple(self._params.keys())

    def list_arguments(self, parameter: str) -> Tuple[str, ...]:
        """Provides a list of arguments required to calculate the parameter.

        Args:
            parameter (str): The parameter of interest.

        Returns:
            A tuple with the arguments.
        """
        sig = signature(self[parameter]).parameters
        return tuple(sig.keys())

    def get_parameter(self, material: str, parameter: str, **kwargs) -> Parameter:
        """Retrieve the parameter for the material.

        Any arguments that obtaining this parameter requires must be included as
        keyword arguments in the call.

        Args:
            material (str): Material the enquiry is about.
            parameter (str): The parameter of interest.
            **kwargs: Any other argument needed to calculate the requested parameter.

        Raises
            ParameterError if the material does not exist in this source or if the
                parameter does not exist for this material.

        Returns:
            A Parameter object with the requested parameter.
        """
        calculable = self[parameter]

        # These are the input arguments required by the calculable
        sig = signature(calculable).parameters

        # Some might come from the input kwargs of this method
        params = {p: kwargs[p] for p in sig.keys() if p in kwargs}

        # Others are parameters that need to be retrieved for the material
        to_retrieve = [p for p in sig.keys() if p not in kwargs]
        for p in to_retrieve:
            # Is it another calculable?
            if p in self.parameters():
                params[p] = self.get_parameter(material, p, **kwargs)

            else:
                # If not, is it in another source?
                try:
                    params[p] = self.parsys.get_parameter(material, p, **kwargs)

                except ParameterError as err:
                    # If not... is there a default value??
                    if sig[p].default != insParam.empty:
                        params[p] = sig[p].default

                    # Well, bad luck then...
                    raise ParameterError(err)

        out = calculable(**params)
        references = chain.from_iterable(
            {p.r for p in params.values() if isinstance(p, Parameter)}
        )
        references = tuple(set(chain(("Calculable",), references)))
        return Parameter(
            out, description=self._descriptions[parameter], reference=references
        )


def eg(T: Quantity, eg0: Quantity, alpha: Quantity, beta: Quantity) -> Quantity:
    """Energy gap as a function of temperature

    Calculate the energy gap for temperature T using the Varshni relationship:

        Eg = Eg0 + alpha * T**2 / (T + beta)

    Args:
        T: Temperature (in K)
        eg0: Energy gap at T=0 (in eV)
        alpha: Proportionality constant (in eV/K)
        beta: Temperature offset (in eV)

    Returns:
        The gap at the chosen temperature (in eV)
    """
    Tk = T.to("K")
    out = eg0 - alpha * Tk ** 2 / (Tk + beta)
    return out.to("eV")


@CalculableParameters.register_calculable(description="Band gap at the Gamma point")
def eg_gamma(
    T: Quantity, eg0_gamma: Quantity, alpha_gamma: Quantity, beta_gamma: Quantity
) -> Quantity:
    """Energy gap at the Gamma point as a function of temperature.

    Calculate the energy gap for temperature T using the Varshni relationship:

        Eg = Eg0 + alpha * T**2 / (T + beta)

    Args:
        T: Temperature (in K)
        eg0_gamma: Energy gap at T=0 (in eV)
        alpha_gamma: Proportionality constant (in eV/K)
        beta_gamma: Temperature offset (in eV)

    Returns:
        The gap at the chosen temperature (in eV)
    """
    return eg(T, eg0_gamma, alpha_gamma, beta_gamma)


@CalculableParameters.register_calculable(description="Band gap at the X point")
def eg_x(T: Quantity, eg0_x: Quantity, alpha_x: Quantity, beta_x: Quantity) -> Quantity:
    """Energy gap at the X point as a function of temperature.

    Calculate the energy gap for temperature T using the Varshni relationship:

        Eg = Eg0 + alpha * T**2 / (T + beta)

    Args:
        T: Temperature (in K)
        eg0_x: Energy gap at T=0 (in eV)
        alpha_x: Proportionality constant (in eV/K)
        beta_x: Temperature offset (in eV)

    Returns:
        The gap at the chosen temperature (in eV)
    """
    return eg(T, eg0_x, alpha_x, beta_x)


@CalculableParameters.register_calculable(description="Band gap at the L point")
def eg_l(T: Quantity, eg0_l: Quantity, alpha_l: Quantity, beta_l: Quantity) -> Quantity:
    """Energy gap at the L point as a function of temperature.

    Calculate the energy gap for temperature T using the Varshni relationship:

        Eg = Eg0 + alpha * T**2 / (T + beta)

    Args:
        T: Temperature (in K)
        eg0_l: Energy gap at T=0 (in eV)
        alpha_l: Proportionality constant (in eV/K)
        beta_l: Temperature offset (in eV)

    Returns:
        The gap at the chosen temperature (in eV)
    """
    return eg(T, eg0_l, alpha_l, beta_l)


@CalculableParameters.register_calculable(description="Band gap energy")
def band_gap(
    eg_gamma: Optional[Quantity] = None,
    eg_x: Optional[Quantity] = None,
    eg_l: Optional[Quantity] = None,
) -> Quantity:
    """Band gap energy, taken as the minimum of the Gamma, X and L points.

    Raises
        ValueError: If the gap of at least one of the points is not provided

    Args:
        eg_gamma: Energy gap at the gamma point (in eV)
        eg_x: Energy gap at the x point (in eV)
        eg_l: Energy gap at the l point (in eV)

    Returns:
        The band gap (in eV)
    """
    gaps = [g.to("eV") for g in (eg_gamma, eg_x, eg_l) if g is not None]

    if len(gaps) == 0:
        raise ValueError("The gap for at least one of the points need to be provided.")

    return min(gaps)


@CalculableParameters.register_calculable(
    description="Lowest point in the conduction band"
)
def lowest_band(
    band_gap: Quantity, eg_gamma: Quantity, eg_x: Quantity, eg_l: Quantity
) -> str:
    """Label indicating what is the lowest band out of [Gamma, X, L].

    Args:
        band_gap: Band gap energy (in eV)
        eg_gamma: Energy gap at the gamma point (in eV)
        eg_x: Energy gap at the x point (in eV)
        eg_l: Energy gap at the l point (in eV)

    Returns:
        The band gap (in eV)
    """
    return ["Gamma", "X", "L"][[eg_gamma, eg_x, eg_l].index(band_gap)]


@CalculableParameters.register_calculable(description="Split-off hole effective mass")
def eff_mass_split_off(
    gamma1: Union[Quantity, float],
    interband_matrix_element: Quantity,
    spin_orbit_splitting: Quantity,
    band_gap: Quantity,
) -> Parameter:
    """Split-off hole effective mass, m_so.

    Provided by Eq. 2.18 of Vurgaftman et al. JAP, 2001:

        m0/m_so = g1 - Ep*Delta_so / (3*Eg*(Eg+Delta_so))

    Args:
        gamma1: Luttinger parameter gamma1 (dimensionless)
        interband_matrix_element: Interband matrix element (eV)
        spin_orbit_splitting: Splitting of the split-off band (eV)
        band_gap: Bandgap of the material (eV)

    Returns:
        The effective mass (in kg)
    """
    mr = gamma1 - interband_matrix_element * spin_orbit_splitting / (
        3 * band_gap * (band_gap + spin_orbit_splitting)
    )
    return electron_mass * mr


@CalculableParameters.register_calculable(
    description="Heavy hole effective mass along the z [110] direction"
)
def eff_mass_hh_z(
    gamma1: Union[Quantity, float], gamma2: Union[Quantity, float]
) -> Parameter:
    """Heavy hole effective mass along the z direction [100].

    Provided by Eq. 2.16 of Vurgaftman et al. JAP, 2001:

        m0/m_hh_z = g1 - 2*g2

    Args:
        gamma1: Luttinger parameter gamma1 (dimensionless)
        gamma2: Luttinger parameter gamma2 (dimensionless)

    Returns:
        The effective mass (in kg)
    """
    return electron_mass * (gamma1 - 2 * gamma2)


@CalculableParameters.register_calculable(
    description="Heavy hole effective mass along the [110] direction"
)
def eff_mass_hh_110(
    gamma1: Union[Quantity, float],
    gamma2: Union[Quantity, float],
    gamma3: Union[Quantity, float],
) -> Parameter:
    """Heavy hole effective mass along the [110] direction.

    Provided by Eq. 2.16 of Vurgaftman et al. JAP, 2001:

        m0/m_hh_111 = 1/2 * (2*g1 - g2 - 3*g3)

    Args:
        gamma1: Luttinger parameter gamma1 (dimensionless)
        gamma2: Luttinger parameter gamma2 (dimensionless)
        gamma3: Luttinger parameter gamma3 (dimensionless)

    Returns:
        The effective mass (in kg)
    """
    return electron_mass * (1 / 2 * (2 * gamma1 - gamma2 - 3 * gamma3))


@CalculableParameters.register_calculable(
    description="Heavy hole effective mass along the [111] direction"
)
def eff_mass_hh_111(
    gamma1: Union[Quantity, float], gamma3: Union[Quantity, float]
) -> Parameter:
    """Heavy hole effective mass along the [111] direction.

    Provided by Eq. 2.16 of Vurgaftman et al. JAP, 2001:

        m0/m_hh_110 = g1 - 2*g3

    Args:
        gamma1: Luttinger parameter gamma1 (dimensionless)
        gamma3: Luttinger parameter gamma3 (dimensionless)

    Returns:
        The effective mass (in kg)
    """
    return electron_mass * (gamma1 - 2 * gamma3)


@CalculableParameters.register_calculable(
    description="Light hole effective mass along the z [100] direction"
)
def eff_mass_lh_z(
    gamma1: Union[Quantity, float], gamma2: Union[Quantity, float]
) -> Parameter:
    """Light hole effective mass along the z direction [100].

    Provided by Eq. 2.17 of Vurgaftman et al. JAP, 2001:

        m0/m_hh_z = g1 + 2*g2

    Args:
        gamma1: Luttinger parameter gamma1 (dimensionless)
        gamma2: Luttinger parameter gamma2 (dimensionless)

    Returns:
        The effective mass (in kg)
    """
    return electron_mass * (gamma1 + 2 * gamma2)


@CalculableParameters.register_calculable(
    description="Light hole effective mass along the [110] direction"
)
def eff_mass_lh_110(
    gamma1: Union[Quantity, float],
    gamma2: Union[Quantity, float],
    gamma3: Union[Quantity, float],
) -> Parameter:
    """Light hole effective mass along the [110] direction.

    Provided by Eq. 2.17 of Vurgaftman et al. JAP, 2001:

        m0/m_hh_111 = 1/2 * (2*g1 + g2 + 3*g3)

    Args:
        gamma1: Luttinger parameter gamma1 (dimensionless)
        gamma2: Luttinger parameter gamma2 (dimensionless)
        gamma3: Luttinger parameter gamma3 (dimensionless)

    Returns:
        The effective mass (in kg)
    """
    return electron_mass * (1 / 2 * (2 * gamma1 + gamma2 + 3 * gamma3))


@CalculableParameters.register_calculable(
    description="Light hole effective mass along the [111] direction"
)
def eff_mass_lh_111(
    gamma1: Union[Quantity, float], gamma3: Union[Quantity, float]
) -> Parameter:
    """Light hole effective mass along the [111] direction.

    Provided by Eq. 2.17 of Vurgaftman et al. JAP, 2001:

        m0/m_hh_110 = g1 + 2*g3

    Args:
        gamma1: Luttinger parameter gamma1 (dimensionless)
        gamma3: Luttinger parameter gamma3 (dimensionless)

    Returns:
        The effective mass (in kg)
    """
    return electron_mass * (gamma1 + 2 * gamma3)


@CalculableParameters.register_calculable(description="Electron effective mass")
def eff_mass_electron(
    f: Union[Quantity, float],
    interband_matrix_element: Quantity,
    band_gap: Quantity,
    spin_orbit_splitting: Quantity,
) -> Parameter:
    """Electron effective mass

    Provided by Eq. 2.15 of Vurgaftman et al. JAP, 2001:

        m0/m_e = 1 + 2*F + Ep*(Eg + 2*Delta_so/3) / (Eg*(Eg+Delta_so))


    Args:
        f: Kane parameter (dimensionless)
        interband_matrix_element: Interband matrix element (eV)
        band_gap: Bandgap of the material (eV)
        spin_orbit_splitting: Splitting of the split-off band (eV)

    Returns:
        The effective mass (in kg)
    """

    mr = (
        1
        + 2 * f
        + interband_matrix_element
        * (band_gap + 2 * spin_orbit_splitting / 3)
        / (band_gap * (band_gap + spin_orbit_splitting))
    )
    return electron_mass * mr


@CalculableParameters.register_calculable(description="Absolute permittivity")
def permittivity(relative_permittivity: Union[Quantity, float]) -> Parameter:
    """Absolute permittivity of the material.

        epsilon = epsilon_r * epsilon_0

    Args:
        relative_permittivity: Relative permittivity (dimensionless)

    Returns:
        The absolute permittivity
    """
    return vacuum_permittivity * relative_permittivity


if __name__ == "__main__":
    print(
        CalculableParameters().get_parameter("GaAs", "band_gap", T=Quantity(298, "K"))
    )
    print(CalculableParameters().list_arguments("eff_mass_electron"))
    print(CalculableParameters().permittivity(9.1))
