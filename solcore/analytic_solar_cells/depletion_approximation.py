import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import solve_bvp

from solcore.constants import kb, q, vacuum_permittivity, pi, hbar, electron_mass
from solcore.science_tracker import science_reference
from solcore.state import State


def iv_depletion(junction, options):
    """ Calculates the IV curve of a junction object using the depletion approximation as described in J. Nelson, “The Physics of Solar Cells”, Imperial College Press (2003). The junction is then updated with an "iv" function that calculates the IV curve at any voltage.

    :param junction: A junction object.
    :param options: Solver options.
    :return: None.
    """

    science_reference('Depletion approximation',
                      'J. Nelson, “The Physics of Solar Cells”, Imperial College Press (2003).')

    junction.voltage = options.internal_voltages
    T = options.T

    # First we have to figure out if we are talking about a PN, NP, PIN or NIP junction
    sn = 0 if not hasattr(junction, "sn") else junction.sn
    sp = 0 if not hasattr(junction, "sp") else junction.sp

    # We search for the emitter and check if it is n-type or p-type
    idx = 0
    pn_or_np = 'pn'

    for layer in junction:
        if layer.role is not 'emitter':
            idx += 1
        else:
            Na = 0
            Nd = 0
            if hasattr(layer.material, 'Na'): Na = layer.material.Na
            if hasattr(layer.material, 'Nd'): Nd = layer.material.Nd
            if Na < Nd:
                pn_or_np = "np"
                nRegion = junction[idx]
            else:
                pRegion = junction[idx]

            id_top = idx

            break

    # Now we check for an intrinsic region and, if there is, for the base.
    if junction[idx + 1].role is 'intrinsic':
        iRegion = junction[idx + 1]

        if junction[idx + 2].role is 'base':
            if pn_or_np == "pn":
                nRegion = junction[idx + 2]
                nidx = idx + 2
            else:
                pRegion = junction[idx + 2]
                pidx = idx + 2

            id_bottom = idx + 2

        else:
            raise RuntimeError(
                'ERROR processing junctions: A layer following the "intrinsic" layer must be defined as '
                '"base".')

    # If there is no intrinsic region, we check directly the base
    elif junction[idx + 1].role is 'base':
        if pn_or_np == "pn":
            nRegion = junction[idx + 1]
            nidx = idx + 1
        else:
            pRegion = junction[idx + 1]
            pidx = idx + 1
        iRegion = None

        id_bottom = idx + 1

    else:
        raise RuntimeError(
            'ERROR processing junctions: A layer following the "emitter" must be defined as "intrinsic"'
            'or "base".')

    # With all regions identified, it's time to start doing calculations
    kbT = kb * T
    Egap = nRegion.material.band_gap
    R_shunt = min(junction.R_shunt, 1e14) if hasattr(junction, 'R_shunt') else 1e14

    xp = pRegion.width
    xn = nRegion.width
    xi = 0 if iRegion is None else iRegion.width

    # Now we have to get all the material parameters needed for the calculation
    if hasattr(junction, "dielectric_constant"):
        es = junction.dielectric_constant
    else:
        es = nRegion.material.permittivity * vacuum_permittivity  # equal for n and p.  I hope.

    # For the diffusion lenght, subscript n and p refer to the carriers, electrons and holes
    if hasattr(junction, "ln"):
        ln = junction.ln
    else:
        ln = pRegion.material.electron_diffusion_length

    if hasattr(junction, "lp"):
        lp = junction.lp
    else:
        lp = nRegion.material.hole_diffusion_length

    # For the diffusion coefficient, n and p refer to the regions, n side and p side. Yeah, it's confusing...
    if hasattr(junction, "mup"):
        dp = junction.mup * kbT / q
    else:
        dp = pRegion.material.electron_mobility * kbT / q

    if hasattr(junction, "mun"):
        dn = junction.mun * kbT / q
    else:
        dn = nRegion.material.hole_mobility * kbT / q

    # Effective masses and effective density of states
    mEff_h = nRegion.material.eff_mass_hh_z * electron_mass
    mEff_e = pRegion.material.eff_mass_electron * electron_mass

    Nv = 2 * (mEff_h * kb * T / (2 * pi * hbar ** 2)) ** 1.5  # Jenny p58
    Nc = 2 * (mEff_e * kb * T / (2 * pi * hbar ** 2)) ** 1.5
    niSquared = Nc * Nv * np.exp(-Egap / (kb * T))
    ni = np.sqrt(niSquared)

    Na = pRegion.material.Na
    Nd = nRegion.material.Nd
    Vbi = (kbT / q) * np.log(Nd * Na / niSquared) if not hasattr(junction, "Vbi") else junction.Vbi  # Jenny p146

    # And now we account for the possible applied voltage, which can be, at most, equal to Vbi
    V = np.where(junction.voltage < Vbi - 0.001, junction.voltage, Vbi - 0.001)

    # It's time to calculate the depletion widths
    if not hasattr(junction, "wp") or not hasattr(junction, "wn"):

        if hasattr(junction, "depletion_approximation") and junction.depletion_approximation == "one-sided abrupt":
            print("using one-sided abrupt junction approximation for depletion width")
            science_reference("Sze abrupt junction approximation",
                              "Sze: The Physics of Semiconductor Devices, 2nd edition, John Wiley & Sons, Inc (2007)")
            wp = np.sqrt(2 * es * (Vbi - V) / (q * Na))
            wn = np.sqrt(2 * es * (Vbi - V) / (q * Nd))

        else:
            wn = (-xi + np.sqrt(xi ** 2 + 2. * es * (Vbi - V) / q * (1 / Na + 1 / Nd))) / (1 + Nd / Na)
            wp = (-xi + np.sqrt(xi ** 2 + 2. * es * (Vbi - V) / q * (1 / Na + 1 / Nd))) / (1 + Na / Nd)

    wn = wn if not hasattr(junction, "wn") else junction.wn
    wp = wp if not hasattr(junction, "wp") else junction.wp

    w = wn + wp + xi

    # Now it is time to calculate currents
    if pn_or_np == "pn":
        l_top, l_bottom = ln, lp
        x_top, x_bottom = xp, xn
        w_top, w_bottom = wp, wn
        s_top, s_bottom = sp, sn
        d_top, d_bottom = dp, dn
        min_top, min_bot = niSquared / Na, niSquared / Nd
    else:
        l_bottom, l_top = ln, lp
        x_bottom, x_top = xp, xn
        w_bottom, w_top = wp, wn
        s_bottom, s_top = sp, sn
        d_bottom, d_top = dp, dn
        min_bot, min_top = niSquared / Na, niSquared / Nd

    JtopDark = get_j_top(x_top, w_top, l_top, s_top, d_top, V, min_top, T)
    JbotDark = get_j_bot(x_bottom, w_bottom, l_bottom, s_bottom, d_bottom, V, min_bot, T)

    # hereby we define the subscripts to refer to the layer in which the current is generated:
    if pn_or_np == "pn":
        JnDark, JpDark = JbotDark, JtopDark
    else:
        JpDark, JnDark = JbotDark, JtopDark

    # These might not be the right lifetimes. Actually, they are not as they include all recombination processes, not
    # just SRH recombination, which is what the equation in Jenny, p159 refers to. Let´ leave them, for now.
    lifetime_n = ln ** 2 / dn
    lifetime_p = lp ** 2 / dp  # Jenny p163

    # Here we use the full version of the SRH recombination term as calculated by Sah et al. Works for positive bias
    # and moderately negative ones.
    science_reference('SRH current term.',
                      'C. T. Sah, R. N. Noyce, and W. Shockley, “Carrier Generation and Recombination in P-N Junctions and P-N Junction Characteristics,” presented at the Proceedings of the IRE, 1957, vol. 45, no. 9, pp. 1228–1243.')
    Jrec = get_Jsrh(ni, V, Vbi, lifetime_p, lifetime_n, w, kbT)

    J_sc_top = 0
    J_sc_bot = 0
    J_sc_scr = 0

    if options.light_iv:

        widths = []
        for layer in junction:
            widths.append(layer.width)

        cum_widths = np.cumsum([0] + widths)

        g = junction.absorbed
        wl = options.wavelength
        wl_sp, ph = options.light_source.spectrum(output_units='photon_flux_per_m', x=wl)
        id_v0 = np.argmin(abs(V))

        # The contribution from the Emitter (top side).
        xa = cum_widths[id_top]
        xb = cum_widths[id_top + 1] - w_top[id_v0]
        deriv = get_J_sc_diffusion(xa, xb, g, d_top, l_top, min_top, s_top, wl, ph, side='top')
        J_sc_top = q * d_top * abs(deriv)

        # The contribution from the Base (bottom side).
        xa = cum_widths[id_bottom] + w_bottom[id_v0]
        xb = cum_widths[id_bottom + 1]
        deriv = get_J_sc_diffusion(xa, xb, g, d_bottom, l_bottom, min_bot, s_bottom, wl, ph, side='bottom')
        J_sc_bot = q * d_bottom * abs(deriv)

        # The contribution from the SCR (includes the intrinsic region, if present).
        xa = cum_widths[id_top + 1] - w_top[id_v0]
        xb = cum_widths[id_bottom] + w_bottom[id_v0]
        J_sc_scr = q * get_J_sc_SCR(xa, xb, g, wl, ph)

    # And, finally, we output the currents
    junction.current = Jrec + JnDark + JpDark + V / R_shunt - J_sc_top - J_sc_bot - J_sc_scr
    junction.iv = interp1d(junction.voltage, junction.current, kind='linear', bounds_error=False, assume_sorted=True,
                           fill_value=(junction.current[0], junction.current[-1]))
    junction.region_currents = State({"Jn_dif": JnDark, "Jp_dif": JpDark, "Jscr_srh": Jrec,
                                      "J_sc_top": J_sc_top, "J_sc_bot": J_sc_bot, "J_sc_scr": J_sc_scr})


def get_j_top(x, w, l, s, d, V, minority, T):
    # We calculate some fractions
    harg = (x - w) / l
    sinh_harg = np.sinh(harg)
    cosh_harg = np.cosh(harg)
    lsod = (l * s) / d

    # And then we are ready to calculate the different currents
    # Missing the voltage dependent part of these equations.
    # They should be 6.34 and 6.39, not 6.62 and 6.63

    J_top_dark = (q * d * minority / l) * (np.exp(q * V / kb / T) - 1) * \
                 ((lsod * cosh_harg + sinh_harg) / (lsod * sinh_harg + cosh_harg))

    return J_top_dark


def get_j_bot(x, w, l, s, d, V, minority, T):
    # We calculate some fractions
    harg = (x - w) / l
    cosh_harg = np.cosh(harg)
    sinh_harg = np.sinh(harg)
    lsod = (l * s) / d

    # And then we are ready to calculate the different currents
    # Missing the voltage dependent part of these equations.
    # They should be 6.34 and 6.39, not 6.62 and 6.63

    J_bottom_dark = (q * d * minority / l) * (np.exp(q * V / kb / T) - 1) * \
                    ((lsod * cosh_harg + sinh_harg) / (lsod * sinh_harg + cosh_harg))

    return J_bottom_dark


def get_Jsrh(ni, V, Vbi, tp, tn, w, kbT, dEt=0):
    science_reference('SRH current term.',
                      'C. T. Sah, R. N. Noyce, and W. Shockley, “Carrier Generation and Recombination in P-N Junctions and P-N Junction Characteristics,” presented at the Proceedings of the IRE, 1957, vol. 45, no. 9, pp. 1228–1243.')

    out = np.zeros(V.shape)

    m = V >= -1120 * kbT / q
    out[m] = forward(ni, V[m], Vbi, tp, tn, w[m], kbT)

    return out


def forward(ni, V, Vbi, tp, tn, w, kbT, dEt=0):
    """ Equation 27 of Sah's paper. Strictly speaking, it is not valid for intermediate negative bias. """
    J0 = 2 * q * ni * w / np.sqrt(tn * tp)
    f = factor(V, Vbi, tp, tn, kbT, dEt)
    out = J0 * np.sinh(q * V / (2 * kbT)) / (q * (Vbi - V) / kbT) * f

    return out


def factor(V, Vbi, tp, tn, kbT, dEt=0):
    """ The integral of Eq. 27 in Sah's paper. While it is coninuum (in principle) it has to be done in two parts.
    (or three) """
    trap = q * dEt / kbT + np.log(tp / tn) / 2
    b = np.exp(-q * V / kbT / 2) * np.cosh(trap)

    z1 = np.sqrt(tp / tn) * np.exp(-q * (Vbi - V) / kbT / 2)
    z2 = np.sqrt(tp / tn) * np.exp(q * (Vbi - V) / kbT / 2)

    out = np.zeros(b.shape)

    # For b values < 1
    m = b < 1
    top = np.arctan((b[m] + z2[m]) / np.sqrt(1 - b[m] ** 2))
    bot = np.arctan((b[m] + z1[m]) / np.sqrt(1 - b[m] ** 2))

    out[m] = (top - bot) / np.sqrt(1 - b[m] ** 2)

    # For b values >= 1 and <=1e7
    m = (b >= 1) * (b <= 1e7)
    xx = b[m]
    yy = np.sqrt(xx ** 2 - 1)

    top = np.log(abs(z2[m] - yy + xx) / abs(z2[m] + yy + xx))
    bot = np.log(abs(z1[m] - yy + xx) / abs(z1[m] + yy + xx))

    out[m] = (top - bot) / (2 * yy)

    # For b values > 1e7
    m = b >= 1e7

    top = np.log(z2[m]) - np.log(z2[m] + 2 * b[m])
    bot = np.log(z1[m]) - np.log(z1[m] + 2 * b[m])

    out[m] = (top - bot) / (2 * b[m])

    return out


def get_J_sc_diffusion(xa, xb, g, D, L, y0, S, wl, ph, side='top'):
    zz = np.linspace(xa, xb, 1001)
    gg = g(zz) * ph
    g_vs_z = np.trapz(gg, wl, axis=1)

    A = lambda x: np.interp(x, zz, g_vs_z) / D + y0 / L ** 2

    def fun(x, y):
        out1 = y[1]
        out2 = y[0] / L ** 2 - A(x)
        return np.vstack((out1, out2))

    if side == 'top':
        def bc(ya, yb):
            left = ya[1] - S / D * (ya[0] - y0)
            right = yb[0]
            return np.array([left, right])
    else:
        def bc(ya, yb):
            left = ya[0]
            right = yb[1] - S / D * (yb[0] - y0)
            return np.array([left, right])

    guess = y0 * np.ones((2, zz.size))
    guess[1] = np.zeros_like(guess[0])

    solution = solve_bvp(fun, bc, zz, guess)

    if side == 'top':
        out = solution.y[1][-1]
    else:
        out = solution.y[1][0]

    return out


def get_J_sc_SCR(xa, xb, g, wl, ph):
    zz = np.linspace(xa, xb, 1001)
    gg = g(zz) * ph
    out = np.trapz(np.trapz(gg, wl, axis=1), zz)

    return out


def qe_depletion(junction, options):
    """ Calculates the QE curve of a junction object using the depletion approximation as described in J. Nelson, “The Physics of Solar Cells”, Imperial College Press (2003). The junction is then updated with an "iqe" and several "eqe" functions that calculates the QE curve at any wavelength.

    :param junction: A junction object.
    :param options: Solver options.
    :return: None.
    """

    science_reference('Depletion approximation',
                      'J. Nelson, “The Physics of Solar Cells”, Imperial College Press (2003).')

    T = options.T

    # First we have to figure out if we are talking about a PN, NP, PIN or NIP junction
    sn = 0 if not hasattr(junction, "sn") else junction.sn
    sp = 0 if not hasattr(junction, "sp") else junction.sp

    # We search for the emitter and check if it is n-type or p-type
    idx = 0
    pn_or_np = 'pn'

    for layer in junction:
        if layer.role is not 'emitter':
            idx += 1
        else:
            Na = 0
            Nd = 0
            if hasattr(layer.material, 'Na'): Na = layer.material.Na
            if hasattr(layer.material, 'Nd'): Nd = layer.material.Nd
            if Na < Nd:
                pn_or_np = "np"
                nRegion = junction[idx]
            else:
                pRegion = junction[idx]

            id_top = idx

            break

    # Now we check for an intrinsic region and, if there is, for the base.
    if junction[idx + 1].role is 'intrinsic':
        iRegion = junction[idx + 1]

        if junction[idx + 2].role is 'base':
            if pn_or_np == "pn":
                nRegion = junction[idx + 2]
                nidx = idx + 2
            else:
                pRegion = junction[idx + 2]
                pidx = idx + 2

            id_bottom = idx + 2

        else:
            raise RuntimeError(
                'ERROR processing junctions: A layer following the "intrinsic" layer must be defined as '
                '"base".')

    # If there is no intrinsic region, we check directly the base
    elif junction[idx + 1].role is 'base':
        if pn_or_np == "pn":
            nRegion = junction[idx + 1]
            nidx = idx + 1
        else:
            pRegion = junction[idx + 1]
            pidx = idx + 1
        iRegion = None

        id_bottom = idx + 1

    else:
        raise RuntimeError(
            'ERROR processing junctions: A layer following the "emitter" must be defined as "intrinsic"'
            'or "base".')

    # With all regions identified, it's time to start doing calculations
    kbT = kb * T
    Egap = nRegion.material.band_gap
    R_shunt = min(junction.R_shunt, 1e14) if hasattr(junction, 'R_shunt') else 1e14

    xp = pRegion.width
    xn = nRegion.width
    xi = 0 if iRegion is None else iRegion.width

    # Now we have to get all the material parameters needed for the calculation
    if hasattr(junction, "dielectric_constant"):
        es = junction.dielectric_constant
    else:
        es = nRegion.material.permittivity * vacuum_permittivity  # equal for n and p.  I hope.

    # For the diffusion lenght, subscript n and p refer to the carriers, electrons and holes
    if hasattr(junction, "ln"):
        ln = junction.ln
    else:
        ln = pRegion.material.electron_diffusion_length

    if hasattr(junction, "lp"):
        lp = junction.lp
    else:
        lp = nRegion.material.hole_diffusion_length

    # For the diffusion coefficient, n and p refer to the regions, n side and p side. Yeah, it's confusing...
    if hasattr(junction, "mup"):
        dp = junction.mup * kbT / q
    else:
        dp = pRegion.material.electron_mobility * kbT / q

    if hasattr(junction, "mun"):
        dn = junction.mun * kbT / q
    else:
        dn = nRegion.material.hole_mobility * kbT / q

    # Effective masses and effective density of states
    mEff_h = nRegion.material.eff_mass_hh_z * electron_mass
    mEff_e = pRegion.material.eff_mass_electron * electron_mass

    Nv = 2 * (mEff_h * kb * T / (2 * pi * hbar ** 2)) ** 1.5  # Jenny p58
    Nc = 2 * (mEff_e * kb * T / (2 * pi * hbar ** 2)) ** 1.5
    niSquared = Nc * Nv * np.exp(-Egap / (kb * T))
    ni = np.sqrt(niSquared)

    Na = pRegion.material.Na
    Nd = nRegion.material.Nd
    Vbi = (kbT / q) * np.log(Nd * Na / niSquared) if not hasattr(junction, "Vbi") else junction.Vbi  # Jenny p146

    # It's time to calculate the depletion widths
    if not hasattr(junction, "wp") or not hasattr(junction, "wn"):

        if hasattr(junction, "depletion_approximation") and junction.depletion_approximation == "one-sided abrupt":
            print("using one-sided abrupt junction approximation for depletion width")
            science_reference("Sze abrupt junction approximation",
                              "Sze: The Physics of Semiconductor Devices, 2nd edition, John Wiley & Sons, Inc (2007)")
            wp = np.sqrt(2 * es * Vbi / (q * Na))
            wn = np.sqrt(2 * es * Vbi / (q * Nd))

        else:
            wn = (-xi + np.sqrt(xi ** 2 + 2. * es * Vbi / q * (1 / Na + 1 / Nd))) / (1 + Nd / Na)
            wp = (-xi + np.sqrt(xi ** 2 + 2. * es * Vbi / q * (1 / Na + 1 / Nd))) / (1 + Na / Nd)

    wn = wn if not hasattr(junction, "wn") else junction.wn
    wp = wp if not hasattr(junction, "wp") else junction.wp

    w = wn + wp + xi

    # Now it is time to calculate currents
    if pn_or_np == "pn":
        l_top, l_bottom = ln, lp
        x_top, x_bottom = xp, xn
        w_top, w_bottom = wp, wn
        s_top, s_bottom = sp, sn
        d_top, d_bottom = dp, dn
        min_top, min_bot = niSquared / Na, niSquared / Nd
    else:
        l_bottom, l_top = ln, lp
        x_bottom, x_top = xp, xn
        w_bottom, w_top = wp, wn
        s_bottom, s_top = sp, sn
        d_bottom, d_top = dp, dn
        min_bot, min_top = niSquared / Na, niSquared / Nd

    widths = []
    for layer in junction:
        widths.append(layer.width)

    cum_widths = np.cumsum([0] + widths)

    g = junction.absorbed
    wl = options.wavelength
    wl_sp, ph = options.light_source.spectrum(output_units='photon_flux_per_m', x=wl)

    # The contribution from the Emitter (top side).
    xa = cum_widths[id_top]
    xb = cum_widths[id_top + 1] - w_top
    deriv = get_J_sc_diffusion_vs_WL(xa, xb, g, d_top, l_top, min_top, s_top, wl, ph, side='top')
    j_sc_top = d_top * abs(deriv)

    # The contribution from the Base (bottom side).
    xa = cum_widths[id_bottom] + w_bottom
    xb = cum_widths[id_bottom + 1]
    deriv = get_J_sc_diffusion_vs_WL(xa, xb, g, d_bottom, l_bottom, min_bot, s_bottom, wl, ph, side='bottom')
    j_sc_bot = d_bottom * abs(deriv)

    # The contribution from the SCR (includes the intrinsic region, if present).
    xa = cum_widths[id_top + 1] - w_top
    xb = cum_widths[id_bottom] + w_bottom
    j_sc_scr = get_J_sc_SCR_vs_WL(xa, xb, g, wl, ph)

    # The total ligth absorbed, but not necesarily collected, is:
    xa = cum_widths[id_top]
    xb = cum_widths[id_bottom + 1]
    zz = np.linspace(xa, xb, 3001)
    gg = g(zz) * ph
    current_absorbed = np.trapz(gg, zz, axis=0)

    # Now, we put everything together
    j_sc = j_sc_top + j_sc_bot + j_sc_scr

    eqe = j_sc / ph
    eqe_emitter = j_sc_top / ph
    eqe_base = j_sc_bot / ph
    eqe_scr = j_sc_scr / ph

    junction.iqe = interp1d(wl, j_sc / current_absorbed)
    junction.eqe = interp1d(wl, eqe, kind='linear', bounds_error=False, assume_sorted=True,
                            fill_value=(eqe[0], eqe[-1]))
    junction.eqe_emitter = interp1d(wl, eqe_emitter, kind='linear', bounds_error=False, assume_sorted=True,
                                    fill_value=(eqe_emitter[0], eqe_emitter[-1]))
    junction.eqe_base = interp1d(wl, eqe_base, kind='linear', bounds_error=False, assume_sorted=True,
                                 fill_value=(eqe_base[0], eqe_base[-1]))
    junction.eqe_scr = interp1d(wl, eqe_scr, kind='linear', bounds_error=False, assume_sorted=True,
                                fill_value=(eqe_scr[0], eqe_scr[-1]))


def get_J_sc_SCR_vs_WL(xa, xb, g, wl, ph):
    zz = np.linspace(xa, xb, 1001)
    gg = g(zz) * ph
    out = np.trapz(gg, zz, axis=0)

    return out


def get_J_sc_diffusion_vs_WL(xa, xb, g, D, L, y0, S, wl, ph, side='top'):
    zz = np.linspace(xa, xb, 1001)
    gg = g(zz) * ph

    out = np.zeros_like(wl)

    for i in range(len(wl)):
        A = lambda x: np.interp(x, zz, gg[:, i]) / D + y0 / L ** 2

        def fun(x, y):
            out1 = y[1]
            out2 = y[0] / L ** 2 - A(x)
            return np.vstack((out1, out2))

        if side == 'top':
            def bc(ya, yb):
                left = ya[1] - S / D * (ya[0] - y0)
                right = yb[0]
                return np.array([left, right])
        else:
            def bc(ya, yb):
                left = ya[0]
                right = yb[1] - S / D * (yb[0] - y0)
                return np.array([left, right])

        guess = y0 * np.ones((2, zz.size))
        guess[1] = np.zeros_like(guess[0])

        solution = solve_bvp(fun, bc, zz, guess)

        if side == 'top':
            out[i] = solution.y[1][-1]
        else:
            out[i] = solution.y[1][0]

    return out
