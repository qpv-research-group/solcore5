import math
from scipy.interpolate import interp1d
import numpy as np

from solcore.constants import kb, q, c, h
from solcore import si
from solcore.state import State

db_options = State()
db_options.db_mode = 'boltzmann'


def iv_detailed_balance(junction, options):
    """ Calculates the IV curve of a junction of kind "DB". This is a detailed balanced calculation of the IV curve of a PN junction characterized by a certain eqe, temperature and chemical potential (voltage). Normally, the eqe will be calculated based on a a given absorption edge and absorptance level resulting in a "top hat" type of eqe, but it can also be provided externally. By default, the solver uses the Boltzmann aproximation, although the full Planck equation might be used, also.

    :param junction: A Junction object of kind "DB"
    :param options: Other arguments for the calculation of the IV.
    :return:
    """
    T = options.T
    light = options.light_iv
    mode = options.db_mode
    Ta = options.T_ambient
    wl = options.wavelength

    try:
        Eg = junction.Eg
        n = junction.n

        if not hasattr(junction, 'eqe'):
            qe_detailed_balance(junction, wl)

        eqe = junction.eqe
        As = surface_integral(junction, wl)

        R_shunt = min(junction.R_shunt, 1e14) if hasattr(junction, 'R_shunt') else 1e14

    except AttributeError as err:
        raise AttributeError(
            'ERROR calculating the IV for the detailed balance Junction kind. Junction is missing one essential argument: {}'.format(
                err))

    j01 = 2 * q * n ** 2 * c
    # We loop over the (internal) voltages to get the (internal) current. The actual voltage applied is always smaller than Eg.
    junction.voltage = options.internal_voltages

    if mode == 'boltzmann':
        # If Boltzmann, the calculation can be vectorized as the voltages is outside the integral. No limit to V
        volt = junction.voltage

        def emis(v, T):
            boltz = j01 * (As / wl ** 4) * np.exp(- (h * c / wl) / (kb * T))
            junction.j01 = np.trapz(boltz, wl)
            out = np.trapz(boltz, wl) * np.exp((q * v) / (kb * T))
            return out

    else:
        # If the full Planck equation is used, we need to loop over the voltages. We limit them to the bandgap
        volt = np.where(junction.voltage < Eg - 0.001, junction.voltage, Eg - 0.001)

        def emis(v, T):
            out = []
            for vv in v:
                fermi = (As / wl ** 4) / (np.exp((h * c / wl - q * vv) / (kb * T)) - 1)
                out.append(np.trapz(fermi, wl))
            return j01*np.array(out)

    # Now we calculate the generation current, made of thermal generation and photogeneration
    jthermal = emis(np.array([0]), Ta)[0]
    if light:
        wl, ph = options.light_source.spectrum(x=wl, output_units='photon_flux_per_m')
        jsc = q * np.trapz(eqe(wl) * ph, wl) + jthermal
    else:
        jsc = jthermal

    junction.current = emis(volt, T) + volt / R_shunt - jsc

    junction.iv = interp1d(junction.voltage, junction.current, kind='linear', bounds_error=False, assume_sorted=True,
                           fill_value=(junction.current[0], junction.current[-1]))


def qe_detailed_balance(junction, wl):
    """ Calculates the EQE and the IQE in the detailed balanced limit... which is equal to the absorptance and the absorbed fraction of the light since, by definition of detailed balanced, the carrier collection is perfect.

    :param junction: Junction object of kind "DB"
    :param wl: wavelength in m
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

    junction.qe = State({'WL': wl, 'IQE': junction.iqe(wl), 'EQE': junction.eqe(wl)})


def absorptance_detailed_balance(junction):
    """ Calculates the absorptance of the junction in the detailed balanced case. If it has not been calculated before, provided as input, it uses a "top hat" absorptance with the absorption edge at Eg and magnitude A.

    Note that while the junction.absorptance represents the potential fraction of light that can be absorbed at a given energy, the junction.absorbed represents the actual fraction of input light that is absorbed. For detailed balanced calculations, these two are equivalent to the IQE and EQE respectively because carrier collection is perfect, but that will not be the case, in general.

    :param junction: Junction object of kind "DB"
    """
    try:
        if not hasattr(junction, 'absorptance'):
            wlg = si(1240 / junction.Eg, 'nm')
            A = min(junction.A, 1)

            def absorptance(wl):
                return A * (wl < wlg)

            junction.absorptance = absorptance

    except AttributeError as err:
        raise AttributeError(
            'ERROR calculating absorptance for the detailed balance Junction kind.\nJunction is missing one essential argument: {}'.format(
                err))


def surface_integral(junction, wl):
    """ Calculates the surface integral of the absorptivity

    :param junction:
    :param options:
    :return:
    """
    Rn = np.minimum(0.999999, junction.reflected(wl))

    # Out of the normal incidence reflection we calculate an effective refractive index
    neff = (1 + np.sqrt(Rn)) / (1 - np.sqrt(Rn))

    # And the corresponding critical angle
    cos_tc = np.sqrt(1 - 1 / neff ** 2)
    steps = int(math.ceil(20 / (1 - min(cos_tc))))
    cos_t = np.linspace(0.001, 1, steps)

    # And now, we back calculate the reflection at any internal angle
    R = np.ones((len(wl), steps))

    for i in range(len(wl)):
        idx = np.where(cos_t > cos_tc[i])

        Rs = np.ones(steps)
        Rp = np.ones(steps)

        cos_i = np.sqrt(1 - neff[i] ** 2 * (1 - cos_t[idx] ** 2))

        Rs[idx] = ((cos_i - neff[i] * cos_t[idx]) / (cos_i + neff[i] * cos_t[idx])) ** 2
        Rp[idx] = ((cos_t[idx] - neff[i] * cos_i) / (cos_t[idx] + neff[i] * cos_i)) ** 2

        R[i, :] = (Rs + Rp) / 2

    # B = exp(-alpha*w) but it is not necessary to calculate alpha and w explicitly
    B = 1 - junction.absorptance(wl)
    A_back = 0.5 * (1 + R[:, -1] * B) * (1 - B)
    A_front = 0.5 * (1 - R[:, -1]) * (1 - B)

    for i in range(1, steps - 1):
        t = cos_t[i]

        # The integral of the absorpton/emission across the back surface
        A_back = A_back + (1 + R[:, i] * B ** (1 / t)) * (1 - B ** (1 / t)) * t

        # The integral of the absorpton/emission across the front surface
        A_front = A_front + (1 - R[:, i]) * (1 - B ** (1 / t)) * t

    A_back = A_back / (steps + 1)
    A_front = A_front / (steps + 1)

    return 2 * np.pi * (A_back + A_front)


if __name__ == '__main__':
    import numpy as np
    import matplotlib.pyplot as plt

    from solcore.state import State
    from solcore import material

    junction = State()
    options = State()

    wl = np.linspace(300, 900, 300) * 1e-9
    GaAs = material('GaAs')(T=300)

    n = GaAs.n(wl)
    cos_tc = np.sqrt(1 - 1 / n ** 2)
    Rn = ((n - 1) / (n + 1)) ** 2
    wlg = 800 * 1e-9
    A = 1

    V = 0
    T = 300

    options.wavelength = wl
    junction.reflected = interp1d(wl, Rn)
    junction.absorptance = lambda x: A * (x < wlg)

    #
    As = surface_integral(junction, options)

    jE = 2 * q * n ** 2 * c * (1 / wl ** 4) * As / (np.exp((h * c / wl - q * V) / kb / T) - 1)

    plt.plot(wl, jE)
    plt.show()


    #
    # all_cos_t = np.linspace(min(cos_tc), 1, 10)
    # Rs = np.zeros_like(n)
    # Rp = np.zeros_like(n)
    # Rnew = np.zeros_like(n)
    # for cos_t in all_cos_t:
    #     for i in range(len(wl)):
    #         if (cos_t > cos_tc[i]):
    #             cos_i = np.sqrt(1 - n[i] ** 2 * (1 - cos_t ** 2))
    #             Rs[i] = ((cos_i - n[i] * cos_t) / (cos_i + n[i] * cos_t)) ** 2
    #             Rp[i] = ((cos_t - n[i] * cos_i) / (cos_t + n[i] * cos_i)) ** 2
    #         else:
    #             Rs[i] = 1
    #             Rp[i] = 1
    #
    #         Rnew[i] = myR(wl[i], cos_t)
    #
    #     plt.plot(wl, (Rs + Rp) / 2, 'o')
    #     plt.plot(wl, Rnew, 'k')
    #
    # plt.show()
