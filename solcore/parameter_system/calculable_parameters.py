"""

"""
from ..constants import q


def eg(T: float, eg0: float, alpha: float, beta: float) -> float:
    """Energy gap as a function of temperature

    Calculate the energy gap for temperature T using the formula:

        Eg = Eg0 + alpha * T**2 / (T + beta)

    Args:
        T: Temperature (in Kelvin)
        eg0: Energy gap at T=0 (in eV)
        alpha: Proportionality constant (in eV/K)
        beta: Offset constant (in K)

    Returns:
        The gap at the chosen temperature (in eV)
    """
    return eg0 - alpha * T ** 2 / (T + beta)


def eg_gamma(
    T: float, eg0_gamma: float, alpha_gamma: float, beta_gamma: float
) -> float:
    """Energy gap at the Gamma point as a function of temperature.

    Calculate the energy gap for temperature T using the formula:

        Eg = Eg0 + alpha * T**2 / (T + beta)

    Args:
        T: Temperature (in Kelvin)
        eg0_gamma: Energy gap at T=0 (in eV)
        alpha_gamma: Proportionality constant (in eV/K)
        beta_gamma: Offset constant (in K)

    Returns:
        The gap at the chosen temperature (in eV)
    """
    return eg(T, eg0_gamma, alpha_gamma, beta_gamma)


def eg_x(T: float, eg0_x: float, alpha_x: float, beta_x: float) -> float:
    """Energy gap at the X point as a function of temperature.

    Calculate the energy gap for temperature T using the formula:

        Eg = Eg0 + alpha * T**2 / (T + beta)

    Args:
        T: Temperature (in Kelvin)
        eg0_x: Energy gap at T=0 (in eV)
        alpha_x: Proportionality constant (in eV/K)
        beta_x: Offset constant (in K)

    Returns:
        The gap at the chosen temperature (in eV)
    """
    return eg(T, eg0_x, alpha_x, beta_x)


def eg_l(T: float, eg0_l: float, alpha_l: float, beta_l: float) -> float:
    """Energy gap at the L point as a function of temperature.

    Calculate the energy gap for temperature T using the formula:

        Eg = Eg0 + alpha * T**2 / (T + beta)

    Args:
        T: Temperature (in Kelvin)
        eg0_l: Energy gap at T=0 (in eV)
        alpha_l: Proportionality constant (in eV/K)
        beta_l: Offset constant (in K)

    Returns:
        The gap at the chosen temperature (in eV)
    """
    return eg(T, eg0_l, alpha_l, beta_l)


def band_gap(eg_gamma: float, eg_x: float, eg_l: float) -> float:
    """Band gap energy, taken as the minimum of the Gamma, X and L bands.

    Args:
        eg_gamma: Energy gap at the gamma point (in eV)
        eg_x: Energy gap at the x point (in eV)
        eg_l: Energy gap at the l point (in eV)

    Returns:
        The band gap (in eV)
    """
    return min(eg_gamma, eg_x, eg_l)


def lowest_band(band_gap: float, eg_gamma: float, eg_x: float, eg_l: float) -> str:
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


"""

[Final Calculables] # parameters that are calculated at the end, when all the bowing is done
m0: get("eff_mass_split_off")*(get("gamma1")-get("interband_matrix_element")*get("spin_orbit_splitting")/(3*get("band_gap")*(get("band_gap")+get("spin_orbit_splitting"))))
eff_mass_hh_z: get('m0')/(get('gamma1')-2*get('gamma2'))
eff_mass_hh_110: get('m0')/(get('gamma1')-0.5*get('gamma2')-1.5*get('gamma3'))
eff_mass_hh_111: get('m0')/(get('gamma1')-2*get('gamma3'))
eff_mass_lh_z: get('m0')/(get('gamma1')+2*get('gamma2'))
eff_mass_lh_110: get('m0')/(get('gamma1')+0.5*get('gamma2')+1.5*get('gamma3'))
eff_mass_lh_111: get('m0')/(get('gamma1')+2*get('gamma3'))
eff_mass_electron: get('m0')/(1+2*get('F')+get('interband_matrix_element')*(get('band_gap')+2*get('spin_orbit_splitting')/3)/(get('band_gap')*(get('band_gap')+get('spin_orbit_splitting'))))
permittivity = 8.854187817e-12*get('relative_permittivity')





"""
