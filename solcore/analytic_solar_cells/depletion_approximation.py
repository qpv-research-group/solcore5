import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import solve_bvp, quad_vec
from functools import partial

from solcore.constants import kb, q
from solcore.science_tracker import science_reference
from solcore.state import State
from solcore.light_source import LightSource

import warnings

da_options = State()
da_options.da_mode = "green"


def identify_layers(junction):
    # First we have to figure out if we are talking about a PN, NP, PIN or NIP junction
    # We search for the emitter and check if it is n-type or p-type
    idx = 0
    pn_or_np = "pn"
    # homojunction = True

    for layer in junction:
        if layer.role.lower() != "emitter":
            idx += 1
        else:
            Na = 0
            Nd = 0
            if hasattr(layer.material, "Na"):
                Na = layer.material.Na
            if hasattr(layer.material, "Nd"):
                Nd = layer.material.Nd
            if Na < Nd:
                pn_or_np = "np"
                nRegion = junction[idx]
            else:
                pRegion = junction[idx]

            id_top = idx

            break

    # Now we check for an intrinsic region and, if there is, for the base.
    if junction[idx + 1].role.lower() == "intrinsic":
        iRegion = junction[idx + 1]

        if junction[idx + 2].role.lower() == "base":
            if pn_or_np == "pn":
                nRegion = junction[idx + 2]

            else:
                pRegion = junction[idx + 2]

            id_bottom = idx + 2
            # homojunction = (
            #   homojunction
            #   and nRegion.material.material_string == pRegion.material.material_string
            # )
            # homojunction = (
            #   homojunction
            #   and nRegion.material.material_string == iRegion.material.material_string
            # )

        else:
            raise RuntimeError("ERROR processing junctions: A layer following the "
                               '"intrinsic" layer must be defined as "base".')

    # If there is no intrinsic region, we check directly the base
    elif junction[idx + 1].role.lower() == "base":
        if pn_or_np == "pn":
            nRegion = junction[idx + 1]

        else:
            pRegion = junction[idx + 1]

        iRegion = None

        id_bottom = idx + 1
        # homojunction = (
        #     homojunction
        #     and nRegion.material.material_string == pRegion.material.material_string
        # )

    else:
        raise RuntimeError('ERROR processing junctions: A layer following the '
                           '"emitter" must be defined as "intrinsic" or "base".')

    # We assert that we are really working with an homojunction
    # assert homojunction, "ERROR: The DA solver only works with homojunctions,
    # for now."

    return id_top, id_bottom, pRegion, nRegion, iRegion, pn_or_np


def identify_parameters(junction, T, pRegion, nRegion, iRegion):
    kbT = kb * T

    xp = pRegion.width
    xn = nRegion.width
    xi = 0 if iRegion is None else iRegion.width

    sn = 0 if not hasattr(junction, "sn") else junction.sn
    sp = 0 if not hasattr(junction, "sp") else junction.sp

    # Now we have to get all the material parameters needed for the calculation
    # labels for permittivity refer to the doping of the layer
    if hasattr(junction, "permittivity"):
        es_n = junction.permittivity
        es_p = junction.permittivity
    else:
        es_n = nRegion.material.permittivity
        es_p = pRegion.material.permittivity

    # For the diffusion length, subscript n and p refer to the carriers,
    # electrons and holes
    if hasattr(junction, "ln"):
        ln = junction.ln
    else:
        ln = pRegion.material.electron_diffusion_length

    if hasattr(junction, "lp"):
        lp = junction.lp
    else:
        lp = nRegion.material.hole_diffusion_length

    # For the diffusion coefficient, n and p refer to the regions,
    # n side and p side. Yeah, it's confusing...
    if hasattr(junction, "mup"):
        dp = junction.mup * kbT / q
    else:
        dp = pRegion.material.electron_mobility * kbT / q

    if hasattr(junction, "mun"):
        dn = junction.mun * kbT / q
    else:
        dn = nRegion.material.hole_mobility * kbT / q

    ni = nRegion.material.ni

    Na = pRegion.material.Na
    Nd = nRegion.material.Nd

    return xn, xp, xi, sn, sp, ln, lp, dn, dp, Nd, Na, ni, es_n, es_p


def iv_depletion(junction, options):
    """Calculates the IV curve of a junction object using the depletion approximation as
    described in J. Nelson, “The Physics of Solar Cells”, Imperial College Press (2003).
    The junction is then updated with an "iv" function that calculates the IV curve at
    any voltage.

    :param junction: A junction object.
    :param options: Solver options.
    :return: None.
    """

    science_reference(
        "Depletion approximation",
        "J. Nelson, “The Physics of Solar Cells”, Imperial College Press (2003).",
    )

    T = options.T
    kbT = kb * T

    id_top, id_bottom, pRegion, nRegion, iRegion, pn_or_np = identify_layers(junction)

    if pn_or_np == "pn":
        junction.voltage = options.internal_voltages
        J_sign = 1

    else:
        junction.voltage = -options.internal_voltages
        J_sign = -1

        if np.all(options.internal_voltages >= 0):
            warnings.warn('All voltages are positive, but junction has been identified as n-p, so the '
                          'open-circuit voltage (Voc) of the junction will be negative.', UserWarning)

        if "voltages" in options.keys():
            if np.all(options.voltages >= 0):
                warnings.warn('All voltages are positive, but junction has been identified as n-p, so the '
                              'open-circuit voltage (Voc) of the junction will be negative.', UserWarning)


    xn, xp, xi, sn, sp, ln, lp, dn, dp, Nd, Na, ni, es_n, es_p = identify_parameters(
        junction, T, pRegion, nRegion, iRegion
    )

    niSquared = ni**2

    Vbi = (
        (kbT / q) * np.log(Nd * Na / niSquared)
        if not hasattr(junction, "Vbi")
        else junction.Vbi
    )  # Jenny p146

    R_shunt = min(junction.R_shunt, 1e14) if hasattr(junction, "R_shunt") else 1e14

    # And now we account for the possible applied voltage,
    # which can be, at most, equal to Vbi
    V = np.where(junction.voltage < Vbi - 0.001, junction.voltage, Vbi - 0.001)

    wn, wp = get_depletion_widths(junction, es_n, es_p, Vbi, V, Na, Nd, xi)

    # if the depletion region is calculated to be wider than the width of the n/p
    # region itself, the whole region is depleted:

    wn[wn > xn] = xn
    wp[wp > xp] = xp

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

    JtopDark = get_j_dark(x_top, w_top, l_top, s_top, d_top, V, min_top, T)
    JbotDark = get_j_dark(
        x_bottom, w_bottom, l_bottom, s_bottom, d_bottom, V, min_bot, T
    )

    # hereby we define the subscripts to refer to the layer in which
    # the current is generated:
    if pn_or_np == "pn":
        JnDark, JpDark = JbotDark, JtopDark
    else:
        JpDark, JnDark = JbotDark, JtopDark

    # These might not be the right lifetimes. Actually, they are not as
    # they include all recombination processes, not just SRH recombination,
    # which is what the equation in Jenny, p159 refers to. Let´s leave them, for now.
    lifetime_n = ln**2 / dn
    lifetime_p = lp**2 / dp  # Jenny p163

    # Here we use the full version of the SRH recombination term as calculated by
    # Sah et al. Works for positive bias and moderately negative ones.

    science_reference(
        "SRH current term.",
        "C. T. Sah, R. N. Noyce, and W. Shockley, “Carrier Generation and "
        "Recombination in P-N Junctions and P-N Junction Characteristics,” presented "
        "at the Proceedings of the IRE, 1957, vol. 45, no. 9, pp. 1228–1243.",
    )
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
        wl_sp, ph = options.light_source.spectrum(
            output_units="photon_flux_per_m", x=wl
        )
        id_v0 = np.argmin(abs(V))

        # The contribution from the Emitter (top side).
        xa = cum_widths[id_top]
        xb = cum_widths[id_top + 1] - w_top[id_v0]

        if options.da_mode == "bvp":
            deriv = get_J_sc_diffusion(
                xa, xb, g, d_top, l_top, min_top, s_top, wl, ph, side="top"
            )
        else:
            xbb = xb - (xb - xa) / 1001.0

            zz = np.linspace(xa, xb, 1001, endpoint=False)
            gg = ph * g(zz)
            g_vs_z = np.trapezoid(gg, wl, axis=1)
            g_vs_z[np.isnan(g_vs_z)] = 0
            g_vs_z = interp1d(zz, g_vs_z, axis=0)

            deriv = get_J_sc_diffusion_green(
                xa, xbb, g_vs_z, d_top, l_top, s_top, 1, side="top"
            )

        J_sc_top = q * d_top * abs(deriv)

        # The contribution from the Base (bottom side).
        xa = cum_widths[id_bottom] + w_bottom[id_v0]
        xb = cum_widths[id_bottom + 1]
        if options.da_mode == "bvp":
            deriv = get_J_sc_diffusion(
                xa, xb, g, d_bottom, l_bottom, min_bot, s_bottom, wl, ph, side="bottom"
            )
        else:
            xbb = xb - (xb - xa) / 1001.0
            zz = np.linspace(xa, xb, 1001, endpoint=False)
            gg = ph * g(zz)
            g_vs_z = np.trapezoid(gg, wl, axis=1)
            g_vs_z[np.isnan(g_vs_z)] = 0
            g_vs_z = interp1d(zz, g_vs_z, axis=0)

            deriv = get_J_sc_diffusion_green(
                xa, xbb, g_vs_z, d_bottom, l_bottom, s_bottom, 1, side="bottom"
            )
            # photogeneration is included in g_vs_z, so set ph=1

        J_sc_bot = q * d_bottom * abs(deriv)

        # The contribution from the SCR (includes the intrinsic region, if present).
        xa = cum_widths[id_top + 1] - w_top[id_v0]
        xb = cum_widths[id_bottom] + w_bottom[id_v0]
        J_sc_scr = q * get_J_sc_SCR(xa, xb, g, wl, ph)

    # And, finally, we output the currents
    junction.voltage = J_sign*junction.voltage # this flips the voltage sign back for an n-p junction, or leaves it
    # unchanged for a p-n junction

    junction.current = J_sign*(
        Jrec + JnDark + JpDark + V / R_shunt - J_sc_top - J_sc_bot - J_sc_scr
    )
    junction.iv = interp1d(
        junction.voltage,
        junction.current,
        kind="linear",
        bounds_error=False,
        assume_sorted=True,
        fill_value=(junction.current[0], junction.current[-1]),
    )
    junction.region_currents = State(
        {
            "Jn_dif": J_sign*JnDark,
            "Jp_dif": J_sign*JpDark,
            "Jscr_srh": J_sign*Jrec,
            "J_sc_top": J_sign*J_sc_top,
            "J_sc_bot": J_sign*J_sc_bot,
            "J_sc_scr": J_sign*J_sc_scr,
        }
    )


def get_j_dark(x, w, L, s, d, V, minority, T):
    """
    :param x: width of top junction
    :param w: depletion width in top junction
    :param L: diffusion length
    :param s: surface recombination velocity
    :param d: diffusion coefficient
    :param V: voltage
    :param minority: minority carrier density
    :param T: Temperature

    :return: J_top_dark
    """

    harg = (x - w) / L

    # only do cosh and sinh if harg is < 200 (avoid inf/nan errors)
    harg[harg > 200] = 200

    sinh_harg = np.sinh(harg)
    cosh_harg = np.cosh(harg)
    lsod = (L * s) / d

    # And then we are ready to calculate the different currents
    # Missing the voltage dependent part of these equations.
    # They should be 6.34 and 6.39, not 6.62 and 6.63

    J_dark = (
        (q * d * minority / L)
        * (np.exp(q * V / kb / T) - 1)
        * ((lsod * cosh_harg + sinh_harg) / (lsod * sinh_harg + cosh_harg))
    )

    return J_dark


def get_Jsrh(ni, V, Vbi, tp, tn, w, kbT, dEt=0):
    science_reference(
        "SRH current term.",
        "C. T. Sah, R. N. Noyce, and W. Shockley, “Carrier Generation and Recombination"
        " in P-N Junctions and P-N Junction Characteristics,” presented at the "
        "Proceedings of the IRE, 1957, vol. 45, no. 9, pp. 1228–1243.",
    )

    out = np.zeros(V.shape)

    m = V >= -1120 * kbT / q
    out[m] = forward(ni, V[m], Vbi, tp, tn, w[m], kbT)

    return out


def forward(ni, V, Vbi, tp, tn, w, kbT, dEt=0):
    """Equation 27 of Sah's paper. Strictly speaking, it is not valid for
    intermediate negative bias."""

    J0 = 2 * q * ni * w / np.sqrt(tn * tp)
    f = factor(V, Vbi, tp, tn, kbT, dEt)
    out = J0 * np.sinh(q * V / (2 * kbT)) / (q * (Vbi - V) / kbT) * f

    return out


def factor(V, Vbi, tp, tn, kbT, dEt=0):
    """The integral of Eq. 27 in Sah's paper. While it is continuum (in principle) it
    has to be done in two parts (or three)"""
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
    yy = np.sqrt(xx**2 - 1)

    top = np.log(abs(z2[m] - yy + xx) / abs(z2[m] + yy + xx))
    bot = np.log(abs(z1[m] - yy + xx) / abs(z1[m] + yy + xx))

    out[m] = (top - bot) / (2 * yy)

    # For b values > 1e7
    m = b >= 1e7

    top = np.log(z2[m]) - np.log(z2[m] + 2 * b[m])
    bot = np.log(z1[m]) - np.log(z1[m] + 2 * b[m])

    out[m] = (top - bot) / (2 * b[m])

    return out


def get_J_sc_diffusion(xa, xb, g, D, L, y0, S, wl, ph, side="top"):
    """
    :param xa:
    :param xb:
    :param g:
    :param D:
    :param L:
    :param y0:
    :param S:
    :param wl:
    :param ph:
    :param side:

    :return: out
    """

    # see comments in get_J_sc_diffusion_vs_WL for details on how differential
    # equations are solved
    zz = np.linspace(xa, xb, 1001, endpoint=False)
    gg = g(zz) * ph

    g_vs_z = np.trapezoid(gg, wl, axis=1)

    g_vs_z[np.isnan(g_vs_z)] = 0

    def A(x):
        return np.interp(x, zz, g_vs_z) / D + y0 / L**2

    def fun(x, y):
        out1 = y[1]
        out2 = y[0] / L**2 - A(x)
        return np.vstack((out1, out2))

    if side == "top":

        def bc(ya, yb):
            left = ya[1] - S / D * (ya[0] - y0)
            right = yb[0] - y0
            return np.array([left, right])

    else:

        def bc(ya, yb):
            left = ya[0] - y0
            right = yb[1] + S / D * (yb[0] - y0)
            return np.array([left, right])

    guess = y0 * np.ones((2, zz.size))
    guess[1] = np.zeros_like(guess[0])

    solution = solve_bvp(fun, bc, zz, guess, max_nodes=2 * zz.shape[0])

    if side == "top":
        out = solution.y[1][-1]
    else:
        out = solution.y[1][0]

    if solution.status != 0:
        warnings.warn(
            "Depletion approximation (DA) I-V calculation: "
            "solve_bvp did not converge as expected",
            RuntimeWarning,
        )

    return out


def _conv_exp_top(x, xa, xb, g, L, phoD):
    """Convolution of the carrier generation rate with the approximate Green's function
    kernel at point x. To be used with the numerical integration routine to compute the
    minority carrier derivative on the top edge. This kernel approximates the original
    one when the diffusion length is 2 orders of magnitude higher than the junction
    width by assuming that sinh(x) = cosh(x) = .5 * exp(x).

    :param x: Coordinate in the junction (variable to be integrated).
    :param xa: Coordinate at the start the junction.
    :param xb: Coordinate at the end the junction.
    :param g: Carrier generation rate at point x (expected as function).
    :param L: Diffusion length.
    :param phoD: Light spectrum divided by the diffusion constant D.
    """
    xc = (xa - x) / L
    xv = np.array(
        [
            xa + xb - x,
        ]
    )
    Pkern = -np.exp(xc)
    Gx = g(xv) * phoD
    return Pkern * Gx


def _conv_exp_bottom(x, xa, g, L, phoD):
    """Convolution of the carrier generation rate with the approximate Green's function
    kernel at point x. To be used with the numerical integration routine to compute the
    minority carrier derivative on the bottom edge. This kernel approximates the
    original one when the diffusion length is 2 orders of magnitude higher than
    the junction width by assuming that sinh(x) = cosh(x) = .5 * exp(x).

    :param x: Coordinate in the junction (variable to be integrated).
    :param xa: Coordinate at the start the junction.
    :param g: Carrier generation rate at point x (expected as function).
    :param L: Diffusion length.
    :param phoD: Light spectrum divided by the diffusion constant D.
    """
    xc = (xa - x) / L
    xv = np.array(
        [
            x,
        ]
    )
    Pkern = np.exp(xc)
    Gx = g(xv) * phoD
    return Pkern * Gx


def _conv_green_top(x, xa, xb, g, L, phoD, crvel):
    """Convolution of the carrier generation rate with the Green's function kernel at
    point x. To be used with the numerical integration routine to compute the minority
    carrier derivative on the top edge.

    :param x: Coordinate in the junction (variable to be integrated).
    :param xa: Coordinate at the start the junction.
    :param xb: Coordinate at the end the junction.
    :param g: Carrier generation rate at point x (expected as function).
    :param L: Diffusion length.
    :param phoD: Light spectrum divided by the diffusion constant D.
    :param crvel: Coefficient computed as S / D * L, with S the surface recombination
         velocity.
    """
    xc = (xb - x) / L
    xv = np.array(
        [
            xa + xb - x,
        ]
    )
    Pkern = np.cosh(xc) + crvel * np.sinh(xc)
    Gx = g(xv) * phoD

    return Pkern * Gx


def _conv_green_bottom(x, xb, g, L, phoD, crvel):
    """Convolution of the carrier generation rate with the Green's function kernel at
    point x. To be used with the numerical integration routine to compute the minority
    carrier derivative on the bottom edge.

    :param x: Coordinate in the junction (variable to be integrated).
    :param xb: Coordinate at the end the junction.
    :param g: Carrier generation rate at point x (expected as function).
    :param L: Diffusion length.
    :param phoD: Light spectrum divided by the diffusion constant D.
    :param crvel: Coefficient computed as S / D * L, with S the surface recombination
        velocity.
    """
    xc = (xb - x) / L
    xv = np.array(
        [
            x,
        ]
    )
    Pkern = np.cosh(xc) - crvel * np.sinh(xc)
    Gx = g(xv) * phoD
    return Pkern * Gx


def get_J_sc_diffusion_green(xa, xb, g, D, L, S, ph, side="top"):
    """Computes the derivative of the minority carrier concentration at the edge of the
    junction by approximating the convolution integral resulting from applying the
    Green's function method to the drift-diffusion equation.

    :param xa: Coordinate at the start the junction.
    :param xb: Coordinate at the end the junction.
    :param g: Carrier generation rate at point x (expected as function).
    :param D: Diffusion constant.
    :param L: Diffusion length.
    :param y0: Carrier equilibrium density.
    :param S: Surface recombination velocity.
    :param ph: Light spectrum.
    :param side: String to indicate the edge of interest. Either 'top' or 'bottom'.

    :return: The derivative of the minority carrier concentration at the edge of the
        junction.
    """

    science_reference(
        "DA Green's function method.",
        "T. Vasileiou, J. M. Llorens, J. Buencuerpo, J. M. Ripalda, D. Izzo and "
        "L. Summerer, “Light absorption enhancement and radiation hardening for triple "
        "junction solar cell through bioinspired nanostructures,” "
        "Bioinspir. Biomim., vol. 16, no. 5, pp. 056010, 2021.",
    )

    xbL = (xb - xa) / L
    crvel = S / D * L
    ph_over_D = ph / D

    # if L too low in comparison to junction width, avoid nan's
    if xbL > 1.0e2:
        if side == "top":
            fun = partial(_conv_exp_top, xa=xa, xb=xb, g=g, L=L, phoD=ph_over_D)
        else:
            fun = partial(_conv_exp_bottom, xa=xa, g=g, L=L, phoD=ph_over_D)
        cp = 1.0

    else:
        if side == "top":
            cp = -np.cosh(xbL) - crvel * np.sinh(xbL)

            fun = partial(
                _conv_green_top, xa=xa, xb=xb, g=g, L=L, phoD=ph_over_D, crvel=crvel
            )
        else:
            cp = np.cosh(xbL) + crvel * np.sinh(xbL)

            fun = partial(
                _conv_green_bottom, xb=xb, g=g, L=L, phoD=ph_over_D, crvel=-crvel
            )

    out, err = quad_vec(fun, xa, xb, epsrel=1.0e-5)
    return out.squeeze() / cp


def get_J_sc_SCR(xa, xb, g, wl, ph):
    zz = np.linspace(xa, xb, 1001, endpoint=False)
    gg = g(zz) * ph
    out = np.trapezoid(np.trapezoid(gg, wl, axis=1), zz)

    return out


def qe_depletion(junction, options):
    """Calculates the QE curve of a junction object using the depletion approximation as
    described in J. Nelson, “The Physics of Solar Cells”, Imperial College Press (2003).
    The junction is then updated with an "iqe" and several "eqe" functions that
    calculates the QE curve at any wavelength.

    :param junction: A junction object.
    :param options: Solver options.
    :return: None.
    """

    science_reference(
        "Depletion approximation",
        "J. Nelson, “The Physics of Solar Cells”, Imperial College Press (2003).",
    )

    # First we have to figure out if we are talking about a PN, NP, PIN or NIP junction
    T = options.T
    kbT = kb * T

    id_top, id_bottom, pRegion, nRegion, iRegion, pn_or_np = identify_layers(junction)
    xn, xp, xi, sn, sp, ln, lp, dn, dp, Nd, Na, ni, es_n, es_p = identify_parameters(
        junction, T, pRegion, nRegion, iRegion
    )

    niSquared = ni**2

    Vbi = (
        (kbT / q) * np.log(Nd * Na / niSquared)
        if not hasattr(junction, "Vbi")
        else junction.Vbi
    )  # Jenny p146

    wn, wp = get_depletion_widths(junction, es_n, es_p, Vbi, 0, Na, Nd, xi)

    # if the depletion region is calculated to be wider than the width of the n/p
    # region itself, the whole region is depleted:

    if wn > xn:
        wn = xn

    if wp > xp:
        wp = xp

    # Now it is time to calculate currents
    if pn_or_np == "pn":
        l_top, l_bottom = ln, lp
        w_top, w_bottom = wp, wn
        s_top, s_bottom = sp, sn
        d_top, d_bottom = dp, dn
        min_top, min_bot = niSquared / Na, niSquared / Nd
    else:
        l_bottom, l_top = ln, lp
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
    wl_sp, ph = LightSource(source_type="black body", x=wl, T=6000).spectrum(
        output_units="photon_flux_per_m", x=wl
    )
    # wl_sp, ph = options.light_source.spectrum(output_units='photon_flux_per_m', x=wl)
    ph = 1e10 * np.ones_like(ph)
    # The contribution from the Emitter (top side).
    xa = cum_widths[id_top]
    xb = cum_widths[id_top + 1] - w_top

    if options.da_mode == "bvp":
        deriv = get_J_sc_diffusion_vs_WL(
            xa, xb, g, d_top, l_top, min_top, s_top, wl, ph, side="top"
        )
    else:
        xbb = xb - (xb - xa) / 1001.0
        deriv = get_J_sc_diffusion_green(
            xa, xbb, g, d_top, l_top, s_top, ph, side="top"
        )

    j_sc_top = d_top * abs(deriv)

    # The contribution from the Base (bottom side).
    xa = cum_widths[id_bottom] + w_bottom
    xb = cum_widths[id_bottom + 1]

    if options.da_mode == "bvp":
        deriv = get_J_sc_diffusion_vs_WL(
            xa, xb, g, d_bottom, l_bottom, min_bot, s_bottom, wl, ph, side="bottom"
        )
    else:
        xbb = xb - (xb - xa) / 1001.0
        deriv = get_J_sc_diffusion_green(
            xa, xbb, g, d_bottom, l_bottom, s_bottom, ph, side="bottom"
        )

    j_sc_bot = d_bottom * abs(deriv)

    # The contribution from the SCR (includes the intrinsic region, if present).
    xa = cum_widths[id_top + 1] - w_top
    xb = cum_widths[id_bottom] + w_bottom
    j_sc_scr = get_J_sc_SCR_vs_WL(xa, xb, g, wl, ph)

    # The total light absorbed, but not necessarily collected, is:
    xa = cum_widths[id_top]
    xb = cum_widths[id_bottom + 1]
    zz = np.linspace(xa, xb, 1001)
    gg = g(zz) * ph
    current_absorbed = np.trapezoid(gg, zz, axis=0)

    # why does this happen sometimes?
    # j_sc_top[j_sc_top < 0] = 0
    # j_sc_bot[j_sc_bot < 0] = 0
    # j_sc_scr[j_sc_scr < 0] = 0

    # Now, we put everything together
    j_sc = j_sc_top + j_sc_bot + j_sc_scr

    eqe = j_sc / ph
    eqe_emitter = j_sc_top / ph
    eqe_base = j_sc_bot / ph
    eqe_scr = j_sc_scr / ph

    iqe = np.divide(j_sc, current_absorbed,
                    out=np.zeros_like(j_sc), where=current_absorbed != 0)

    junction.iqe = interp1d(wl, iqe)

    junction.eqe = interp1d(
        wl,
        eqe,
        kind="linear",
        bounds_error=False,
        assume_sorted=True,
        fill_value=(eqe[0], eqe[-1]),
    )
    junction.eqe_emitter = interp1d(
        wl,
        eqe_emitter,
        kind="linear",
        bounds_error=False,
        assume_sorted=True,
        fill_value=(eqe_emitter[0], eqe_emitter[-1]),
    )
    junction.eqe_base = interp1d(
        wl,
        eqe_base,
        kind="linear",
        bounds_error=False,
        assume_sorted=True,
        fill_value=(eqe_base[0], eqe_base[-1]),
    )
    junction.eqe_scr = interp1d(
        wl,
        eqe_scr,
        kind="linear",
        bounds_error=False,
        assume_sorted=True,
        fill_value=(eqe_scr[0], eqe_scr[-1]),
    )
    junction.qe = State(
        {
            "WL": wl,
            "IQE": junction.iqe(wl),
            "EQE": junction.eqe(wl),
            "EQE_emitter": junction.eqe_emitter(wl),
            "EQE_base": junction.eqe_base(wl),
            "EQE_scr": junction.eqe_scr(wl),
        }
    )


def get_J_sc_SCR_vs_WL(xa, xb, g, wl, ph):
    zz = np.linspace(xa, xb, 1001, endpoint=False)
    gg = g(zz) * ph
    out = np.trapezoid(gg, zz, axis=0)

    return out


def get_J_sc_diffusion_vs_WL(xa, xb, g, D, L, y0, S, wl, ph, side="top"):
    zz = np.linspace(xa, xb, 1001, endpoint=False)
    # excluding the last point - depending on the mesh/floating point errors,
    # sometimes this is actually in the next layer

    gg = g(zz) * ph  # generation rate * photon flux

    out = np.zeros_like(wl)
    sol_success = np.zeros_like(wl)  # keep track of whether solve_bvp is converging

    for i in range(len(wl)):
        if np.all(
            gg[:, i] == 0
        ):  # no reason to solve anything if no generation at this wavelength
            out[i] = 0
            sol_success[i] = 0

        else:

            def A(x):
                return np.interp(x, zz, gg[:, i]) / D + y0 / L**2

            # generation and n0/p0 term in differential equation
            # (eq. 6.15 & 6.20 in Jenny Nelson, Physics of Solar Cells)

            def fun(x, y):
                # differential equation (eq. 6.15 & 6.20 in Jenny Nelson,
                # Physics of Solar Cells)
                # solve_bvp solves equation of form:
                # dy / dx = f(x, y, p), a <= x <= b
                # in this case y = [n or p, dn/dx or dp/dx]
                # y[0] = carrier concentration (n or p)
                # y[1] = carrier concentration gradient (dn/dx or dp/dx)
                out1 = y[1]  # by definition! dy/dx = dy/dx

                out2 = y[0] / L**2 - A(x)
                # actually solving the differential equation (6.15 & 6.20)

                return np.vstack((out1, out2))

            # boundary conditions for solve_bvp:
            if side == "top":

                def bc(ya, yb):
                    left = ya[1] - S / D * (ya[0] - y0)
                    # eq. 6.18 - b.c. at front of junction (surface recombination)

                    right = yb[0] - y0
                    # eq. 6.17 - b.c. edge of depletion region, top half of junction
                    # added - y0 (generally very small so makes almost no difference)
                    return np.array([left, right])

            else:

                def bc(ya, yb):
                    left = ya[0] - y0
                    # eq. 6.21 - b.c. edge of depletion region, bottom half of junction
                    # added - y0 (generally very small so makes almost no difference)

                    right = yb[1] + S / D * (yb[0] - y0)
                    # eq. 6.22 - b.c. at back of junction (surface recombination)
                    # changed sign! Current is going the other way
                    return np.array([left, right])

            guess = y0 * np.ones((2, zz.size))
            guess[1] = np.zeros_like(guess[0])

            solution = solve_bvp(fun, bc, zz, guess, max_nodes=2 * zz.shape[0])
            # increase max_nodes to avoid "too many mesh points" message

            sol_success[i] = solution.status

            if side == "top":
                out[i] = solution.y[1][-1]
                # current at edge of depletion region (top half of junction), eq. 6.33
            else:
                out[i] = solution.y[1][0]
                # current at edge of depletion region (bottom half of junction), eq 6.38

        # give a warning f any of the solution statuses are not 0 using warnings.warn:
        if np.any(sol_success != 0):
            warnings.warn(
                "Depletion approximation (DA) EQE calculation: "
                "solve_bvp did not converge as expected for some wavelengths",
                RuntimeWarning,
            )

    return out


def get_depletion_widths(junction, es_n, es_p, Vbi, V, Na, Nd, xi):
    if not hasattr(junction, "wp") or not hasattr(junction, "wn"):
        if (
            hasattr(junction, "depletion_approximation")
            and junction.depletion_approximation == "one-sided abrupt"
        ):
            print("using one-sided abrupt junction approximation for depletion width")
            one_sided = True
        else:
            one_sided = False

        if one_sided:
            science_reference(
                "Sze abrupt junction approximation",
                "Sze: The Physics of Semiconductor Devices, "
                "2nd edition, John Wiley & Sons, Inc (2007)",
            )
            wn = np.sqrt(2 * es_n * (Vbi - V) / (q * Nd))
            wp = np.sqrt(2 * es_p * (Vbi - V) / (q * Na))

        else:
            wn = (
                -xi + np.sqrt(xi**2 + 2.0 * es_n * (Vbi - V) / q * (1 / Na + 1 / Nd))
            ) / (1 + Nd / Na)
            wp = (
                -xi + np.sqrt(xi**2 + 2.0 * es_p * (Vbi - V) / q * (1 / Na + 1 / Nd))
            ) / (1 + Na / Nd)

    wn = wn if not hasattr(junction, "wn") else junction.wn
    wp = wp if not hasattr(junction, "wp") else junction.wp

    return wn, wp


# def iv_depletion_analytical(junction, options):
#     """Calculates the IV curve of a junction object using the analytical solutions to
#     the depletion approximation with Beer-Lambert absorption as described in e.g. J.
#     Nelson, “The Physics of Solar Cells”, Imperial College Press (2003). Note that the
#     equations used here are not exactly those from the book, as they contain typos.
#     The junction is then updated with an "iv" function that calculates the IV curve at
#     any voltage.
#
#     :param junction: A junction object.
#     :param options: Solver options.
#     :return: None.
#     """
#
#     pass


# def minority_carrier_top_analytical(x, L, alpha, xa, D, S, n0, A):
#     aL = alpha * L
#     a1 = -aL * np.exp(2 * alpha * (x + xa) + (2 * x / L)) + aL * np.exp(
#         2 * alpha * (x + xa) + (2 * xa / L)
#     )
#
#     a2 = np.exp(2 * alpha * x + alpha * xa + (xa / L)) - np.exp(
#         alpha * x + 2 * alpha * xa + (x / L)
#     )
#
#     a3 = np.exp(((aL + 1) * (2 * x + xa)) / L) - np.exp(((aL + 1) * (x + 2 * xa)) / L)
#
#     b1 = np.exp(((aL + 1) * (2 * x + xa)) / L) - np.exp(((aL + 1) * (x + 2 * xa)) / L)
#
#     b2 = np.exp(alpha * x + 2 * alpha * xa + (x / L)) - np.exp(
#         2 * alpha * x + alpha * xa + (xa / L)
#     )
#
#     b3 = -np.exp(2 * alpha * (x + xa) + (2 * x / L)) + np.exp(
#         2 * alpha * (x + xa) + (2 * xa / L)
#     )
#
#     numerator = (
#         A
#         * L**2
#         * np.exp(-2 * alpha * (x + xa) - (x / L))
#         * (D * (a1 + a2 + a3) + L * S * (b1 + b2 + b3))
#     )
#     denominator = (
#         D
#         * (alpha**2 * L**2 - 1)
#         * (D * (np.exp(2 * xa / L) + 1) + L * S * (np.exp(2 * xa / L) - 1))
#     )
#     result = numerator / denominator + n0
#
#     return result


# def minority_carrier_bottom_analytical(x, L, alpha, xb, D, S, p0, A):
#     aL = alpha * L
#     denominator = (
#         D
#         * (-1 + alpha**2 * L**2)
#         * (D * (1 + np.exp((2 * xb) / L)) + (-1 + np.exp((2 * xb) / L)) * L * S)
#     )
#
#     a1 = np.exp((2 * xb) / L + 2 * alpha * (x + xb)) - np.exp(
#         (((1 + aL) * (x + 2 * xb)) / L)
#     )
#
#     a2 = np.exp((2 * x) / L + 2 * alpha * (x + xb)) - np.exp(
#         alpha * x + x / L + 2 * alpha * xb
#     )
#
#     a3 = (
#         alpha * np.exp(2 * alpha * x + alpha * xb + xb / L) * L
#         - alpha * np.exp(((1 + aL) * (2 * x + xb)) / L) * L
#     )
#
#     b1 = np.exp(((1 + aL) * (2 * x + xb)) / L) - np.exp(((1 + aL) * (x + 2 * xb)) / L)
#
#     b2 = +np.exp(alpha * x + x / L + 2 * alpha * xb) - np.exp(
#         2 * alpha * x + alpha * xb + xb / L
#     )
#
#     b3 = np.exp((2 * xb) / L + 2 * alpha * (x + xb)) - np.exp(
#         (2 * x) / L + 2 * alpha * (x + xb)
#     )
#
#     numerator = (
#         A
#         * np.exp(-(x / L) - 2 * alpha * (x + xb))
#         * L**2
#         * (D * (a1 + a2 + a3) + (b1 + b2 + b3) * L * S)
#     )
#
#     p = p0 + numerator / denominator
#     return p
