import numpy as np
from solcore.interpolate import interp1d
from solcore.structure import Junction, Layer
from solcore.constants import kb, q, hbar, electron_mass, vacuum_permittivity, pi
from solcore.science_tracker import science_reference
from solcore import eVnm, siUnits, convert

fs = 6.8e-5
Ts = 5762.0



#################
### Old Functions

def spectral_response_all_junctions(solar_cell, incident_light=None, energy=None, V=0, verbose=False):
    """ Calculates the spectral response of any number of junctions using analytical diode equations as described in
    J. Nelson's book "The Physics of Solar Cells" (2003). All parameters must be in SI units. It only works for
    homojunctions.

    :param solar_cell: A structure object with one or more Junction objects and a few optional parameters. Each junction
     is made of a sequence of *Layers* (with a thickness, a material and a role) and the surface recombination
     velocities for electrons and holes. Optional attributes for the structure are:
        - shading: (optional) Shading loses.
        - reflectivity: (optional) Function that calculates the reflectivity as a function of energy (in J).
    :param incident_light: (optional) 2D array containing the energies and the spectral power density (in photons m-2 J-1).
    :param energy: (optional) energies at which to calculate the QE. If not provided, the range of the incident ligth
    is used and if that is not defined, a "sensible" range is used.
    :param V: The voltage at which to perform the calculations (in V)
    :param verbose: If the information about the calculation must be printed (default = FALSE).
    :return: A dictionary containing, as a function of energy:

        - junctions: A list containing the QE data for each junction
        - transmitted: The transmitted power density
        - transmitted_fraction: The fraction of transmitted power
        - passive_loss: The passive losses (light absorbed in the encapsulants, AR coatings, window layers, etc)
        - reflected: reflected power density
        - reflected_fraction: Fraction of reflected power
        - e: The energy values
    """
    science_reference("Nelson pin spectral response",
                      "Jenny: (Nelson. The Physics of Solar Cells. Imperial College Press (2003))")

    # Get the energy range and incident spectrum. If they are not inputs, we create some sensible values.
    if energy is None:
        if incident_light is not None:
            energy = incident_light[0]
            bs = np.copy(incident_light[1])
        else:
            energy = siUnits(np.linspace(0.5, 3.5, 450), 'eV')
            bs = np.ones_like(energy)
    else:
        if incident_light is not None:
            bs = np.interp(energy, incident_light[0], incident_light[1])
        else:
            bs = np.ones_like(energy)

    bs_initial = np.copy(bs)

    # We include the shadowing losses
    if hasattr(solar_cell, 'shading'):
        bs *= (1 - solar_cell.shading)

    # And the reflexion losses
    if hasattr(solar_cell, 'reflectivity') and solar_cell.reflectivity is not None:
        ref = solar_cell.reflectivity(energy)
        bs *= (1 - ref)
        reflected = ref * bs_initial
    else:
        reflected = np.zeros_like(bs)

    # And now we perform the calculation, each junction at a time
    qe_result = []
    passive_loss = np.ones_like(bs)

    for layer_index, layer_object in enumerate(solar_cell):

        # Attenuation due to absorption in the AR coatings or any layer in the front that is not part of the
        # junction
        if type(layer_object) is Layer:
            bs = bs * np.exp(-layer_object.material.alphaE(energy) * layer_object.width)
            passive_loss *= np.exp(-layer_object.material.alphaE(energy) * layer_object.width)

        # For each junction, we calculate the spectral response
        elif type(layer_object) is Junction:

            # If there are window layers, passively absorbing light above the emitter, we attenuate the intensity
            idx = 0
            for junction_layer_object in layer_object:
                if junction_layer_object.role != 'emitter':
                    bs = bs * np.exp(-junction_layer_object.material.alphaE(energy) * junction_layer_object.width)
                    passive_loss *= np.exp(-junction_layer_object.material.alphaE(energy) * junction_layer_object.width)
                    idx += 1
                else:
                    break

            output = calculate_junction_sr(layer_object, energy, bs, bs_initial, V, printParameters=verbose)
            qe_result.append(output)

            # And we reduce the amount of light reaching the next junction
            for junction_layer_object in layer_object[idx:]:
                bs *= np.exp(-junction_layer_object.material.alphaE(energy) * junction_layer_object.width)

        else:
            raise ValueError("Strange layer-like object discovered in structure stack: {}".format(type(layer_object)))

    return {"junctions": qe_result,
            "transmitted": bs,
            "transmitted_fraction": bs / bs_initial,
            "passive_loss": 1 - passive_loss,
            "reflected": reflected,
            "reflected_fraction": reflected / bs_initial,
            "e": energy}


def calculate_junction_sr(junc, energies, bs, bs_initial, V, printParameters=False):
    """ Calculates the total quantum efficiency, the QE splitted by regions, photocurrent and other parameters for a
    given junction at a given voltage.

    :param junc: The junction object
    :param energies: The energies at which to perform the calculation
    :param bs: The spectral power density reaching the junction
    :param bsInitial: The initial power density
    :param V: The voltage at which to perform the calculations (not implemented, yet)
    :param printParameters: If a list of all parameters must be printed
    :return:
    """

    # First we have to figure out if we are talking about a PN, NP, PIN or NIP junction
    sn = 0 if not hasattr(junc, "sn") else junc.sn
    sp = 0 if not hasattr(junc, "sp") else junc.sp

    # First we search for the emitter and check if it is n-type or p-type
    idx = 0
    pn_or_np = 'pn'
    for layer in junc:
        if layer.role is not 'emitter':
            idx += 1
        else:
            Na = 0
            Nd = 0
            if hasattr(layer.material, 'Na'): Na = layer.material.Na
            if hasattr(layer.material, 'Nd'): Nd = layer.material.Nd
            if Na < Nd:
                pn_or_np = "np"
                nRegion = junc[idx]
            else:
                pRegion = junc[idx]
            break

    # Now we check for an intrinsic region and, if there is, for the base.
    if junc[idx + 1].role is 'intrinsic':
        iRegion = junc[idx + 1]

        if junc[idx + 2].role is 'base':
            if pn_or_np == "pn":
                nRegion = junc[idx + 2]
            else:
                pRegion = junc[idx + 2]
        else:
            raise RuntimeError('ERROR processing junctions: A layer following the "intrinsic" layer must be defined as '
                               '"base".')

    # If there is no intrinsic region, we check directly the base
    elif junc[idx + 1].role is 'base':
        if pn_or_np == "pn":
            nRegion = junc[idx + 1]
        else:
            pRegion = junc[idx + 1]
        iRegion = None

    else:
        raise RuntimeError('ERROR processing junctions: A layer following the "emitter" must be defined as "intrinsic"'
                           'or "base".')

    # With all regions identified, it's time to start doing calculations
    T = nRegion.material.T
    kbT = kb * T
    Egap = nRegion.material.band_gap

    xp = pRegion.width
    xn = nRegion.width
    xi = 0 if iRegion is None else iRegion.width

    # Now we have to get all the material parameters needed for the calculation
    if hasattr(junc, "dielectric_constant"):
        es = junc.dielectric_constant
    else:
        es = nRegion.material.permittivity * vacuum_permittivity  # equal for n and p.  I hope.

    # For the diffusion lenght, subscript n and p refer to the carriers, electrons and holes
    if hasattr(junc, "ln"):
        ln = junc.ln
    else:
        ln = pRegion.material.electron_diffusion_length

    if hasattr(junc, "lp"):
        lp = junc.lp
    else:
        lp = nRegion.material.hole_diffusion_length

    # For the diffusion coefficient, n and p refer to the regions, n side and p side. Yeah, it's confusing...
    if hasattr(junc, "mup"):
        dp = junc.mup * kb * T / q
    else:
        dp = pRegion.material.electron_mobility * kb * T / q

    if hasattr(junc, "mun"):
        dn = junc.mun * kb * T / q
    else:
        dn = nRegion.material.hole_mobility * kb * T / q

    # Effective masses and effective density of states
    mEff_h = nRegion.material.eff_mass_hh_z * electron_mass
    mEff_e = pRegion.material.eff_mass_electron * electron_mass

    Nv = 2 * (mEff_h * kb * T / (2 * pi * hbar ** 2)) ** 1.5  # Jenny p58
    Nc = 2 * (mEff_e * kb * T / (2 * pi * hbar ** 2)) ** 1.5
    niSquared = Nc * Nv * np.exp(-Egap / (kb * T))
    ni = np.sqrt(niSquared)

    Na = pRegion.material.Na
    Nd = nRegion.material.Nd
    Vbi = (kb * T / q) * np.log(Nd * Na / niSquared) if not hasattr(junc, "Vbi") else junc.Vbi  # Jenny p146

    # And now we account for the possible applied voltage, which can be, at most, equal to Vbi
    V = min(Vbi, V)
    Vbi = Vbi - V

    # It's time to calculate the depletion widths
    if not hasattr(junc, "wp") or not hasattr(junc, "wn"):

        if hasattr(junc, "depletion_approximation") and junc.depletion_approximation == "one-sided abrupt":
            print("using one-sided abrupt junction approximation for depletion width")
            science_reference("Sze abrupt junction approximation",
                              "Sze: The Physics of Semiconductor Devices, 2nd edition, John Wiley & Sons, Inc (2007)")
            wp = np.sqrt(2 * es * Vbi / (q * Na))
            wn = np.sqrt(2 * es * Vbi / (q * Nd))

        else:
            wn = (-xi + np.sqrt(xi ** 2 + 2. * es * Vbi / q * (1 / Na + 1 / Nd))) / (1 + Nd / Na)
            wp = (-xi + np.sqrt(xi ** 2 + 2. * es * Vbi / q * (1 / Na + 1 / Nd))) / (1 + Na / Nd)

    wn = wn if not hasattr(junc, "wn") else junc.wn
    wp = wp if not hasattr(junc, "wp") else junc.wp

    # we have an array of alpha values that needs to be interpolated to the right energies
    alphaN = nRegion.material.alphaE(energies)  # create numpy array at right energies.
    alphaP = pRegion.material.alphaE(energies)  # create numpy array at right energies.
    alphaI = iRegion.material.alphaE(energies) if iRegion else 0
    depleted_width = wn + wp + xi
    bs_incident_on_top = bs

    # Now it is time to calculate currents
    if pn_or_np == "pn":
        bs_incident_on_bottom = bs * np.exp(-alphaP * xp - alphaN * wn - xi * alphaI)
        bs_incident_on_depleted = bs * np.exp(-alphaP * (xp - wp))
        alphaTop = alphaP
        alphaBottom = alphaN

        l_top, l_bottom = ln, lp
        x_top, x_bottom = xp, xn
        w_top, w_bottom = wp, wn
        s_top, s_bottom = sp, sn
        d_top, d_bottom = dp, dn
        min_top, min_bot = niSquared / Na, niSquared / Nd
    else:
        bs_incident_on_bottom = bs * np.exp(-alphaN * xn - alphaP * wp - xi * alphaI)
        bs_incident_on_depleted = bs * np.exp(-alphaN * (xn - wn))
        alphaTop = alphaN
        alphaBottom = alphaP

        l_bottom, l_top = ln, lp
        x_bottom, x_top = xp, xn
        w_bottom, w_top = wp, wn
        s_bottom, s_top = sp, sn
        d_bottom, d_top = dp, dn
        min_bot, min_top = niSquared / Na, niSquared / Nd

    j_top, JtopDark = get_j_top(x_top, w_top, l_top, s_top, d_top, alphaTop, bs_incident_on_top, V, min_top, T)
    j_bottom, JbotDark = get_j_bot(x_bottom, w_bottom, l_bottom, s_bottom, d_bottom, alphaBottom, bs_incident_on_bottom,
                                   V, min_bot, T)

    jgen = q * bs_incident_on_depleted * (1 - np.exp(-alphaI * xi - alphaN * wn - alphaP * wp))  # jgen. Jenny, p. 159

    # hereby we define the subscripts to refer to the layer in which the current is generated:
    if pn_or_np == "pn":
        jn, jp = j_bottom, j_top
        JnDark, JpDark = JbotDark, JtopDark
    else:
        jp, jn = j_bottom, j_top
        JpDark, JnDark = JbotDark, JtopDark

    # These might not be the right lifetimes. Actually, they are not as they include all recombination processes, not
    # just SRH recombination, which is what the equation in Jenny, p159 refers to. Let´ leave them, for now.
    lifetime_n = ln ** 2 / dn
    lifetime_p = lp ** 2 / dp  # Jenny p163

    # Jrec. Jenny, p.159. Note capital J. This does not need integrating over energies
    Jrec = q * ni * (wn + wp + xi) / np.sqrt(lifetime_n * lifetime_p) * np.sinh(q * V / (2 * kbT)) / (
        q * Vbi / kbT) * pi

    # jgen = q* bs*(1 - exp(-depleted_width*alpha))*exp(-(xn-wn)*alpha);
    nDepletionCharge = wn * Nd * q
    pDepletionCharge = wp * Na * q

    Vbi2 = (0.5 * (wn + wp) + xi) * pDepletionCharge / es

    good_indeces = np.isfinite(jn) * np.isfinite(jp) * np.isfinite(jgen)

    energies = energies[good_indeces]
    jn = jn[good_indeces]
    jp = jp[good_indeces]
    jgen = jgen[good_indeces]

    # jn[jn < 0] = 0
    # jp[jp < 0] = 0
    # jgen[jgen < 0] = 0

    bs_initial = bs_initial[good_indeces]

    Jn = np.trapz(y=jn, x=energies)
    Jp = np.trapz(y=jp, x=energies)
    Jgen = np.trapz(y=jgen, x=energies)

    if printParameters:
        jSum = list((jn + jp + jgen) / bs_initial / q)
        peakQE = max(jSum)
        BandgapEV = convert(Egap, "J", "eV")
        peakQEE = convert(energies[jSum.index(peakQE)], "J", "eV")

        parameterDictionary = {
            "typee": "PIN" if iRegion is not None else "PN",
            "Na": Na, "Nd": Nd,
            "mEff_e": mEff_e, "mEff_h": mEff_h,
            "dp": dp, "dn": dn, "T": T, "es": es, "Vbi": Vbi,
            "xpUm": convert(xp, "m", "um"),
            "xnUm": convert(xn, "m", "um"),
            "BandgapEV": BandgapEV,
            "BandgapNM": eVnm(BandgapEV),
            "pMatrialString": str(pRegion.material),
            "nMatrialString": str(nRegion.material),
            "relativeEffectiveMassE": mEff_e / electron_mass,
            "relativeEffectiveMassH": mEff_h / electron_mass,
            "wnNm": convert(wn, "m", "nm"),
            "wpNm": convert(wp, "m", "nm"),
            "NaCM": convert(Na, "m-3", "cm-3"),
            "NdCM": convert(Nd, "m-3", "cm-3"),
            "permittivity": es / vacuum_permittivity,
            "peakQE": peakQE * 100,
            "peakQEE": peakQEE,
            "peakQENM": eVnm(peakQEE),
            "lpum": convert(lp, "m", "um"),
            "lnum": convert(ln, "m", "um"),
            "nQ": convert(nDepletionCharge, "m-2", "cm-2"),
            "pQ": convert(pDepletionCharge, "m-2", "cm-2"),
            "pDepletionCharge": pDepletionCharge,
            "nDepletionCharge": nDepletionCharge,
            "fieldMax": pDepletionCharge / (es),
            "iRegionString": "",
            "Vbi2": Vbi2,
            "ni": ni
        }
        if iRegion is not None:
            parameterDictionary["iRegionString"] = """
| i region: {xiUm:.2f} um {iMatrialString}
|   field: {fieldMax:.3e} V/m""".format(**{"iMatrialString": str(iRegion.material), "xiUm": convert(xi, "m", "um"),
                                           "fieldMax": pDepletionCharge / (es), })

        print("""\n
| Calculating {typee} QE. Active Parameters:
|
| p region: {xpUm:.2f} um {pMatrialString}
|   Na = {Na:.3e} m-3 ({NaCM:.3e} cm-3)
|   minority carrier (electron) diffusion length: {lpum:.3f} um
|   minority carrier (electron) effective mass: {relativeEffectiveMassE:.3f} (relative) {mEff_e:.3e} kg (absolute)
|   minority carrier (electron) diffusivity = {dp:.3e} m2/s
|   depletion width: {wpNm:.3f} nm
|   charge in depletion region: {pDepletionCharge:.3e} C m-2 ({pQ:.3e} C cm-2){iRegionString}
| n region: {xnUm:.2f} um {nMatrialString}
|   Nd = {Nd:.3e} m-3 ({NdCM:.3e} cm-3)
|   minority carrier (hole) diffusion length: {lnum:.3f} um
|   minority carrier (hole) effective mass: {relativeEffectiveMassH:.3f} (relative) {mEff_h:.3e} kg (absolute)
|   minority carrier (hole) diffusivity = {dn:.3e} m2/s
|   depletion width: {wnNm:.3f} nm
|   charge in depletion region: {nDepletionCharge:.3e} C m-2 ({nQ:.3e} C cm-2)
| Bandgap: {BandgapEV:.3f} eV ({BandgapNM:.2f} nm)
| Temperature: {T:.2f} K
| permittivity: {permittivity:.3f} (relative) {es:.3e} A s V-1 m-1 (absolute)
| built-in Voltage: {Vbi:.3f} V
| peak field: {fieldMax:.3e} V/m
| ni: {ni:e}
| Result:
| \tPeak QE = {peakQE:.1f} % at {peakQEE:.3f} eV ({peakQENM:.2f} nm)""".format(**parameterDictionary))

    return {
        "qe_n": jn / q / bs_initial,
        "qe_p": jp / q / bs_initial,
        "qe_scr": jgen / q / bs_initial,
        "qe_tot": (jn + jp + jgen) / q / bs_initial,
        "Jn_sc": Jn,
        "Jp_sc": Jp,
        "Jscr_sc": Jgen,
        "Jn_dif": JnDark,
        "Jp_dif": JpDark,
        "Jscr_srh": Jrec,
        "J": (Jn + Jp + Jgen - Jrec - JnDark - JpDark),
        "e": energies,
        "Temporary locals dictionary for radiative efficiency": locals()
    }


def get_j_top(x, w, l, s, d, alpha, bs, V, minority, T):
    # We calculate some fractions
    harg = (x - w) / l
    sinh_harg = np.sinh(harg)
    cosh_harg = np.cosh(harg)
    lsod = (l * s) / d

    # And then we are ready to calculate the different currents
    # Missing the voltage dependent part of these equations.
    # They should be 6.34 and 6.39, not 6.62 and 6.63
    j_top_light = (q * bs * alpha * l) / (alpha ** 2 * l ** 2 - 1) * \
                  ((lsod + alpha * l - np.exp(-alpha * (x - w)) * (lsod * cosh_harg + sinh_harg)) /
                   (lsod * sinh_harg + cosh_harg)
                   - alpha * l * np.exp(-alpha * (x - w)))

    J_top_dark = (q * d * minority / l) * (np.exp(q * V / kb / T) - 1) * \
                 ((lsod * cosh_harg + sinh_harg) / (lsod * sinh_harg + cosh_harg))

    return j_top_light, J_top_dark


def get_j_bot(x, w, l, s, d, alpha, bs, V, minority, T):
    # We calculate some fractions
    harg = (x - w) / l
    cosh_harg = np.cosh(harg)
    sinh_harg = np.sinh(harg)
    lsod = (l * s) / d

    # And then we are ready to calculate the different currents
    # Missing the voltage dependent part of these equations.
    # They should be 6.34 and 6.39, not 6.62 and 6.63

    j_bottom_light = (q * bs * alpha * l) / (alpha ** 2 * l ** 2 - 1) * \
                     (l * alpha -
                      (lsod * cosh_harg + sinh_harg - (lsod - l * alpha) * np.exp(-alpha * (x - w))) /
                      (cosh_harg + lsod * sinh_harg))

    J_bottom_dark = (q * d * minority / l) * (np.exp(q * V / kb / T) - 1) * \
                    ((lsod * cosh_harg + sinh_harg) / (lsod * sinh_harg + cosh_harg))

    return j_bottom_light, J_bottom_dark





