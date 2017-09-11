import numpy as np

from solcore.constants import kb, q


def iv_2diode(junction, options):
    """ Calculates the IV curve of a junction object using the 2-diode equation. All parameters needed for the calculation need to be included in the junction object. Series resistance is included at solar cell level, not at junction level.

    The junction is then updated with an "iv" function that calculates the IV curve at any voltage. """

    T = options.T
    light = options.light_iv

    try:
        j01 = junction.j01
        j02 = junction.j02
        n1 = junction.n1
        n2 = junction.n2
        R_shunt = min(junction.R_shunt, 1e14) if hasattr(junction, 'R_shunt') else 1e14

        # If the saturation currents correspond to a different temperature, we update them for the current temperature.
        if hasattr(junction, 'Tref') and (T != junction.Tref):
            assert hasattr(junction, 'Eg'), 'ERROR: The bandgap for each junction (Eg) must be provided if the working ' \
                                        'temperature (T) is different that the reference temperature (Tref). '

            Eg = junction.Eg
            Tref = junction.Tref
            kB = kb / q
            j01 = j01 * (T / Tref) ** 3 * np.exp(-Eg / (n1 * kB) * (1 / T - 1 / Tref))
            j02 = j02 * (T / Tref) ** (5. / 3.) * np.exp(-Eg / (n2 * kB) * (1 / T - 1 / Tref))

        if not light:
            jsc = 0
        else:
            if hasattr(junction, 'jsc'):
                jsc = junction.jsc
            elif hasattr(junction, 'eqe'):
                eqe = junction.eqe
                wl, ph = options.light_source.spectrum(output_units='photon_flux_per_nm')
                jsc = q * np.trapz(eqe(wl) * ph, wl)
            else:
                jsc = 0

    except AttributeError as err:
        raise AttributeError('ERROR in 2-diode equation. Junction is missing one essential argument.') from err

    junction.voltage = options.internal_voltages

    def iv(v):
        out = j01 * (np.exp(q * v / (n1 * kb * T)) - 1) + j02 * (np.exp(q * v / (n2 * kb * T)) - 1) + v / R_shunt - jsc
        return np.minimum(out, 1e8)

    junction.current = iv(junction.voltage)
    junction.iv = iv

