import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt

from solcore.constants import kb, q, c, h
from scipy.interpolate import interp1d


def iv_detailed_balance(junction, options):
    """ Calculates the IV curve of a junction of kind "DB". This is a detailed balanced calculation of the IV curve of a PN junction characterized by a certain eqe, temperature and chemical potential (voltage). Normally, the eqe will be calculated based on a a given absorption edge and absorptance level resulting in a "top hat" type of eqe, but it can also be provided externally.

    :param junction: A Junction object of kind "DB"
    :param options: Other arguments for the calculation of the IV.
    :return:
    """
    T = options.T
    light = options.light_iv
    mode = options.db_mode
    Ta = options.T_ambient
    E_range = options.E_range
    wl = options.wavelength

    try:
        Eg = junction.Eg
        A = min(junction.A, 1)
        n = junction.n

        if not hasattr(junction, 'eqe'):
            qe_detailed_balance(junction, wl)

        eqe = junction.eqe

        R_shunt = min(junction.R_shunt, 1e14) if hasattr(junction, 'R_shunt') else 1e14

    except AttributeError as err:
        raise AttributeError(
            'ERROR calculating the IV for the detailed balance Junction kind. Junction is missing one essential argument: {}'.format(
                err))

    j01 = q * 2 * np.pi * n ** 2 / (h / q) ** 3 / c ** 2

    if mode == 'top_hat':
        def emis(v, T):
            b = kb * T / q
            out = A * b * np.exp(-(Eg - v) / b) * (Eg ** 2 + (2 * Eg * b) + (2 * b ** 2))
            return out

    elif mode == 'boltzmann':
        def emis(v, T):
            b = kb * T / q
            out = quad(lambda x: eqe(1240 / x) * x ** 2 * np.exp(-(x - v) / b), max(Eg - E_range, 0), Eg + E_range)
            return out[0]

    else:
        def emis(v, T):
            b = kb * T / q
            out = quad(lambda x: eqe(1240 / x) * max(x ** 2 / (np.exp((x - v) / b) - 1), 0), max(Eg - E_range, 0),
                       Eg + E_range, limit=100)
            return out[0]

    # Now we calculate the generation current, made of thermal generation and photogeneration
    jthermal = j01 * emis(0, Ta)
    if light:
        wl, ph = options.light_source.spectrum(x=wl, output_units='photon_flux_per_nm')
        jsc = q * np.trapz(eqe(wl) * ph, wl) + jthermal
    else:
        jsc = jthermal

    # We loop over the (internal) voltages to get the (internal) current. The actual voltage applied is always smaller than Eg.
    junction.voltage = options.internal_voltages
    volt = np.where(junction.voltage < Eg - 0.001, junction.voltage, Eg - 0.001)

    # In the case of "top_hat", emis is vectorized so we use it as it is to speed up things
    # For the others, we have to manually loop over the voltages
    try:
        current = emis(volt, T)
    except (TypeError, ValueError):
        current = np.zeros_like(volt)
        for i, vv in enumerate(volt):
            current[i] = emis(vv, T)

    junction.current = j01 * current + volt / R_shunt - jsc

    junction.iv = interp1d(junction.voltage, junction.current, kind='linear', bounds_error=False, assume_sorted=True,
                           fill_value=(junction.current[0], junction.current[-1]))


def qe_detailed_balance(junction, wl):
    """ Calculates the EQE and the IQE in the detailed balanced limit... which is equal to the absorptance and the absorbed fraction of the light since, by definition of detailed balanced, the carrier collection is perfect.

    :param junction: Junction object of kind "DB"
    :return: None
    """
    absorptance_detailed_balance(junction)

    # This is true just in the detailed balance limit
    junction.iqe = junction.absorptance

    try:
        z = np.linspace(0, junction.width, 1001)
        all_abs = junction.absorbed(z)
        abs_vs_wl = np.trapz(all_abs, z, axis=0)
        junction.eqe = interp1d(wl, abs_vs_wl, bounds_error=False, fill_value=(0, 0))

    except AttributeError:
        junction.eqe = junction.absorptance


def absorptance_detailed_balance(junction):
    """ Calculates the absorptance of the junction in the detailed balanced case if it has not been calculated before. It uses a "top hat" absorptance with the absorption edge at Eg and magnitude A.

    Note that while the junction.absorptance represents the potential fraction of light that can be absorbed at a given energy, the junction.absorbed represents the actual fraction of input light that is absorbed. For detailed balanced calculations, these two are equivalent to the IQE and EQE respectively because carrier collection is perfect, but that will not be the case, in general.

    :param junction: Junction object of kind "DB"
    """
    try:
        if not hasattr(junction, 'absorptance'):
            wlg = 1240 / junction.Eg
            A = min(junction.A, 1)

            def absorptance(wl):
                return A * (wl < wlg)

            junction.absorptance = absorptance

    except AttributeError as err:
        raise AttributeError(
            'ERROR calculating absorptance for the detailed balance Junction kind.\nJunction is missing one essential argument: {}'.format(
                err))
