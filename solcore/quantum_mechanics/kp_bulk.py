import numpy as np
from numpy.linalg import eig
from functools import lru_cache

import solcore
from solcore.graphing import *
from solcore import constants
from solcore.crystals import traverse_brillouin, brillouin_critical_points, kvector
from solcore.science_tracker import science_reference

electron_charge = constants.q
hbar = constants.hbar
m0 = constants.electron_mass
q = constants.q


@lru_cache(maxsize=1000, typed=False)
def eight_band_strain_hamiltonian(kx, ky, kz, Ev0, Ec0, exx, ezz, me_eff, gamma1, gamma2, gamma3, a0, Delta, ac, av, b,
                                  Ep):
    """Hamiltonian to calculate cb, hh, lh, so bands and include strain for the biaxial (along [001]) special case.
    
    See Hamiltonian in ref. but remove row 1 and 6 and column 1 and 6.
    http://prb.aps.org/pdf/PRB/v73/i12/e125348"""
    science_reference("k.p Hamiltonian", "Stanko Tomic et al., Electronic structure of InyGa1yAs1xNxGaAs(N) quantum "
                                         "dots by ten-band k.p theory. Phys. Rev. B 73, 125348 (2006)")

    av = abs(
        av)  # The sign in Vurgaftmann is negative while Chuang (and Tomic) assumes it is possitve. We make sure it truly is.

    Eg = Ec0 - Ev0
    sqrt2 = np.sqrt(2)
    sqrt3 = np.sqrt(3)
    sqrt6 = np.sqrt(6)
    sqrt3o2 = np.sqrt(3 / 2)

    m0 = constants.electron_mass

    # The Kane matrix element (Vurgaftman 'interband matrix element')
    P0 = np.sqrt(Ep * constants.hbar ** 2 / (2 * m0))  # as a momentum

    # "Modified" Luttinger parameters
    gc = 1 / (me_eff / constants.electron_mass) - (Ep / 3) * (2 / Eg + 1 / (Eg + Delta))
    g1 = gamma1 - Ep / (
            3 * Eg + Delta)  # There was an error here. It was written as: g1 = gamma1 - Ep/(3*(Eg + Delta))
    g2 = gamma2 - Ep / (6 * Eg + 2 * Delta)
    g3 = gamma3 - Ep / (6 * Eg + 2 * Delta)

    Ck = constants.hbar ** 2 / (2 * m0)

    # Luttinger-Kohn terms
    Ok = lambda kx, ky, kz: Ck * gc * (kx ** 2 + ky ** 2 + kz ** 2)
    Pk = lambda kx, ky, kz: Ck * g1 * (kx ** 2 + ky ** 2 + kz ** 2)
    Qk = lambda kx, ky, kz: Ck * g2 * (kx ** 2 + ky ** 2 - 2 * kz ** 2)
    Rk = lambda kx, ky, kz: Ck * sqrt3 * (g2 * (kx ** 2 - ky ** 2) - 2 * 1j * g3 * kx * ky)
    Sk = lambda kx, ky, kz: Ck * sqrt6 * g3 * (kx - 1j * ky) * kz
    Tk = lambda kx, ky, kz: 1 / sqrt6 * P0 * (kx + 1j * ky)
    Uk = lambda kx, ky, kz: 1 / sqrt3 * P0 * kz

    # Conjugate
    Rkc = lambda kx, ky, kz: Ck * sqrt3 * (g2 * (kx ** 2 - ky ** 2) + 2 * 1j * g3 * kx * ky)
    Skc = lambda kx, ky, kz: Ck * sqrt6 * g3 * (kx + 1j * ky) * kz
    Tkc = lambda kx, ky, kz: 1 / sqrt6 * P0 * (kx - 1j * ky)

    # Strain terms (biaxial approximation makes many of these function zero)
    Oe = lambda exx, eyy, ezz: ac * (exx + eyy + ezz)
    Pe = lambda exx, eyy, ezz: -av * (exx + eyy + ezz)
    Qe = lambda exx, eyy, ezz: -b / 2 * (exx + eyy - 2 * ezz)
    # Re = lambda exx,eyy,ezz: 0
    # Se = lambda exx,eyy,ezz: 0
    # Tk = lambda exx,eyy,ezz: 0
    # Uk = lambda exx,eyy,ezz: 0

    # Complete terms
    O = lambda kx, ky, kz, exx, eyy, ezz: Ok(kx, ky, kz) + Oe(exx, eyy, ezz)
    P = lambda kx, ky, kz, exx, eyy, ezz: Pk(kx, ky, kz) + Pe(exx, eyy, ezz)
    Q = lambda kx, ky, kz, exx, eyy, ezz: Qk(kx, ky, kz) + Qe(exx, eyy, ezz)
    R = lambda kx, ky, kz, exx, eyy, ezz: Rk(kx, ky, kz)  # Here the strain term is zero
    S = lambda kx, ky, kz, exx, eyy, ezz: Sk(kx, ky, kz)  # Here the strain term is zero
    T = lambda kx, ky, kz, exx, eyy, ezz: Tk(kx, ky, kz)  # Here the strain term is zero
    U = lambda kx, ky, kz, exx, eyy, ezz: Uk(kx, ky, kz)  # Here the strain term is zero

    # Conjugates
    Rc = lambda kx, ky, kz, exx, eyy, ezz: Rkc(kx, ky, kz)  # Here the strain term is zero
    Sc = lambda kx, ky, kz, exx, eyy, ezz: Skc(kx, ky, kz)  # Here the strain term is zero
    Tc = lambda kx, ky, kz, exx, eyy, ezz: Tkc(kx, ky, kz)  # Here the strain term is zero

    # Bands
    Ecb = lambda kx, ky, kz, exx, eyy, ezz: Ec0 + O(kx, ky, kz, exx, eyy, ezz)
    Ehh = lambda kx, ky, kz, exx, eyy, ezz: Ev0 - (P(kx, ky, kz, exx, eyy, ezz) + Q(kx, ky, kz, exx, eyy, ezz))
    Elh = lambda kx, ky, kz, exx, eyy, ezz: Ev0 - (P(kx, ky, kz, exx, eyy, ezz) - Q(kx, ky, kz, exx, eyy, ezz))
    Eso = lambda kx, ky, kz, exx, eyy, ezz: Ev0 - (P(kx, ky, kz, exx, eyy, ezz) + Delta)

    eyy = exx
    p = P(kx, ky, kz, exx, eyy, ezz)
    q = Q(kx, ky, kz, exx, eyy, ezz)
    r = R(kx, ky, kz, exx, eyy, ezz)
    s = S(kx, ky, kz, exx, eyy, ezz)
    t = T(kx, ky, kz, exx, eyy, ezz)
    u = U(kx, ky, kz, exx, eyy, ezz)
    rc = Rc(kx, ky, kz, exx, eyy, ezz)
    sc = Sc(kx, ky, kz, exx, eyy, ezz)
    tc = Tc(kx, ky, kz, exx, eyy, ezz)
    cb = Ecb(kx, ky, kz, exx, eyy, ezz)
    hh = Ehh(kx, ky, kz, exx, eyy, ezz)
    lh = Elh(kx, ky, kz, exx, eyy, ezz)
    so = Eso(kx, ky, kz, exx, eyy, ezz)

    H_ST = np.mat([
        [cb, -sqrt3 * t, sqrt2 * u, -u, 0, 0, -tc, -sqrt2 * tc],
        [-sqrt3 * tc, hh, sqrt2 * s, -s, 0, 0, -r, -sqrt2 * r],
        [sqrt2 * u, sqrt2 * sc, lh, -sqrt2 * q, tc, r, 0, sqrt3 * s],
        [-u, -sc, -sqrt2 * q, so, sqrt2 * tc, sqrt2 * r, -sqrt3 * s, 0],
        [0, 0, t, sqrt2 * t, cb, -sqrt3 * tc, sqrt2 * u, -u],
        [0, 0, rc, sqrt2 * rc, -sqrt3 * t, hh, sqrt2 * sc, -sc],
        [-t, -rc, 0, -sqrt3 * sc, sqrt2 * u, sqrt2 * s, lh, -sqrt2 * q],
        [-sqrt2 * t, -sqrt2 * rc, sqrt3 * sc, 0, -u, -s, -sqrt2 * q, so]
    ])

    #     View hamiltonian, good for debugging
    #     if True:
    #         import pylab
    #         pylab.imshow(numpy.real(H_ST))
    #         pylab.show()

    H_ST = H_ST.transpose()
    E, Psi = eig(H_ST)  # do all the math

    bands = np.array((sorted(E.real)))
    return bands


def effective_mass_energy(k, m_eff, E0):
    return E0 + constants.hbar ** 2 / (2 * m_eff) * k ** 2


def kp_bands(material, host_material, kx=0, ky=0, kz=0, return_so=False, graph=False, fit_effective_mass=False,
             effective_mass_direction="X", additional_traverse=None):
    """

    :param material:
    :param host_material:
    :param kx:
    :param ky:
    :param kz:
    :param return_so:
    :param graph:
    :param fit_effective_mass:
    :param effective_mass_direction:
    :param additional_traverse:
    :return:
    """
    Delta = material.spin_orbit_splitting
    Eg = material.band_gap
    Ev0 = material.valence_band_offset
    Ec0 = material.valence_band_offset + Eg
    a0 = material.lattice_constant
    g1 = material.gamma1
    g2 = material.gamma2
    g3 = material.gamma3
    ac = material.a_c
    av = material.a_v
    b = material.b
    Ep = material.interband_matrix_element
    # print (ac, av, b, Ep, a0)
    me_eff = material.eff_mass_electron_Gamma * constants.electron_mass

    # Strain parameters
    exx = - (material.lattice_constant - host_material.lattice_constant) / material.lattice_constant
    ezz = - 2 * material.c12 / material.c11 * exx  # why this?
    # print("suspect stuff around here")
    zMax = np.pi / (a0) / 4
    kx = 0 * zMax
    ky = 0 * zMax

    result = eight_band_strain_hamiltonian(kx, ky, kz, Ev0, Ec0, exx, ezz, me_eff, g1, g2, g3, a0, Delta, ac, av, b, Ep)
    so, lh, hh, c = result[::2]
    if material.lattice_constant < host_material.lattice_constant:
        # print("CAUTION: BLINDLY SWAPPING HH,LH LABELS BASED ON LATTICE CONSTANT")
        lh, hh = hh, lh

    if fit_effective_mass:
        fitpoint_x = 0.05

        critical_points = brillouin_critical_points(material.lattice_constant)
        critical_gamma = critical_points["Gamma"]
        critical_L = critical_points[effective_mass_direction]

        fitpoint = critical_gamma + (critical_L - critical_gamma) * fitpoint_x
        k_cGx, k_cGy, k_cGz = critical_gamma
        k_fx, k_fy, k_fz = fitpoint
        bands_at_gamma = eight_band_strain_hamiltonian(k_cGx, k_cGy, k_cGz, Ev0, Ec0, exx, ezz, me_eff, g1, g2, g3, a0,
                                                       Delta, ac, av, b, Ep)
        bands_at_fitpoint = eight_band_strain_hamiltonian(k_fx, k_fy, k_fz, Ev0, Ec0, exx, ezz, me_eff, g1, g2, g3, a0,
                                                          Delta, ac, av, b, Ep)
        k_max = np.linalg.norm(fitpoint)

        effective_masses = []
        e0s = []
        for band in range(len(result)):
            fpE = bands_at_fitpoint[band]
            oE = bands_at_gamma[band]

            deltaE = fpE - oE
            m_eff_fit = constants.hbar ** 2 / (2 * deltaE) * k_max ** 2
            effective_masses.append(m_eff_fit)
            e0s.append(oE)

        m_eff_so, m_eff_lh, m_eff_hh, m_eff_c = effective_masses[::2]
        if material.lattice_constant < host_material.lattice_constant:
            # print("CAUTION: BLINDLY SWAPPING EFF MASS HH,LH LABELS BASED ON LATTICE CONSTANT")
            m_eff_lh, m_eff_hh = m_eff_hh, m_eff_lh

    if graph:
        traverse_order = ("Gamma", effective_mass_direction)
        if additional_traverse is not None:
            traverse_order = traverse_order + additional_traverse
        allk, x_graph, xticks = traverse_brillouin(a0, traverse_order=traverse_order, steps=300)  # ()

        ticks = [x[0] for x in xticks]
        labels = [x[1] for x in xticks]

        bands = []
        for kx, ky, kz in allk:
            bands.append(
                eight_band_strain_hamiltonian(kx, ky, kz, Ev0, Ec0, exx, ezz, me_eff, g1, g2, g3, a0, Delta, ac, av, b,
                                              Ep))

        g = []

        for i, b in enumerate(zip(*bands)):
            # print (asUnit(array(b),"eV"))
            if fit_effective_mass:
                g.append(GraphData(x_graph, effective_mass_energy(x_graph / fitpoint_x * k_max, effective_masses[i],
                                                                  e0s[i]) / electron_charge, color="red", linewidth=2,
                                   alpha=0.7))

            g.append(GraphData(x_graph, np.array(b) / electron_charge, color="black"))
        a = Graph(g, xticks=ticks, xticklabels=labels, ylabel="$E$ (eV)", xlabel="$k$", palette="hue wheel",
                  ylim=(-2, 1))
        # a = Graph(g, xticks = xticks, ylabel="$E$ (eV)", palette="hue wheel", ylim = (-2,1))
        # a = Graph(g,  ylabel="$E$ (eV)", palette="hue wheel", title="{} on {} -- Luttinger-Kohn w/Pikus-Bir(Stanko)".format(str(material), str(host_material)), ylim = (-2,1))
        a.draw()

    if fit_effective_mass:
        if return_so:
            return (c, hh, lh, so, np.abs(m_eff_c), np.abs(m_eff_hh), np.abs(m_eff_lh), np.abs(m_eff_so))
        else:
            return (c, hh, lh, np.abs(m_eff_c), np.abs(m_eff_hh), np.abs(m_eff_lh))

    if return_so:
        return (c, hh, lh, so)
    else:
        return (c, hh, lh)


def KPbands(material, host_material, return_edges_only=False, plot_result=False, t=0, p=np.pi / 2, fraction=0.2,
            points=50, vin=None):
    """ New version of the above function that produces either the band edges or the full bands in one specifc direction
    from Gamma.
    
    - Calculation of the full Brillouin zone is not allowed. After all, KP is only valid near k=0
    - The calculation of the effective masses is done externally by another function
    - The direction to calculate the bands is provided either with directory angles or a directory vector
    - The number of points of the bands and the fraction up to the border of the Brilluin zone in the given direction are also inputs
    - The default direction is X, as given by the angles theta (t) and phi (p)
    
    """

    Delta = material.spin_orbit_splitting
    Eg = material.band_gap
    Ev0 = material.valence_band_offset
    Ec0 = material.valence_band_offset + Eg
    a0 = material.lattice_constant
    g1 = material.gamma1
    g2 = material.gamma2
    g3 = material.gamma3
    ac = material.a_c
    av = material.a_v
    b = material.b
    Ep = material.interband_matrix_element
    me_eff = material.eff_mass_electron_Gamma * constants.electron_mass

    # Strain parameters
    exx = (host_material.lattice_constant - a0) / a0
    ezz = - 2 * material.c12 / material.c11 * exx

    if return_edges_only:
        # If we only want the band edges, we stop here
        result = eight_band_strain_hamiltonian(0, 0, 0, Ev0, Ec0, exx, ezz, me_eff, g1, g2, g3, a0, Delta, ac, av, b,
                                               Ep)
        so, lh, hh, c = result[::2]

        if a0 < host_material.lattice_constant:
            # We need to figure out a safer way of doing this
            # print("CAUTION: BLINDLY SWAPPING HH,LH LABELS BASED ON LATTICE CONSTANT")
            lh, hh = hh, lh

        return (c, hh, lh, so)

    allk = kvector(a=a0, t=t, p=p, fraction=fraction, points=points, vin=vin)

    bands = []
    for kx, ky, kz in allk:
        bands.append(
            eight_band_strain_hamiltonian(kx, ky, kz, Ev0, Ec0, exx, ezz, me_eff, g1, g2, g3, a0, Delta, ac, av, b, Ep))

    output = [np.array([np.linalg.norm(vect) for vect in allk])]
    for i, band in enumerate(zip(*bands)):
        output.append(np.array(band))

    if plot_result:
        import matplotlib.pyplot as plt
        for i in range(1, len(output)):
            plt.plot(output[0] * a0 / 2 / np.pi, output[i] / q, 'k')
        plt.xlabel('k (2*pi/a)')
        plt.ylabel('Energy (eV)')
        plt.show()

    output = np.array(output)
    return output


def fit_effective_masses(bands, material, host_material, plot_result=False, dk=0.026):
    from scipy.optimize import curve_fit

    masses = []
    a0 = material.lattice_constant

    # General maximum k value. The actual maximum value depends on the direction,
    # but the difference is generally not too big
    kmax = 2 * np.pi / a0

    # We give more weight in the fit to those k points near k=0
    weigth = np.exp(- bands[0] ** 2 / (dk * kmax) ** 2)
    for i in range(1, len(bands), 2):  # Bands are duplicated (at least for now) so we only fit half of them
        band = bands[i]

        # The initial guess for the effective mass (forward finite differences at k = 0)
        p0 = (bands[0][2] - bands[0][1]) * (bands[0][1] - bands[0][0]) / (
                    band[2] - 2 * band[1] + band[0]) * hbar ** 2 / m0

        # We consider only points with weight > 0.01 to speed up the algorithm and prevent convergence errors
        popt, pcov = curve_fit(parabolic, bands[0][weigth > 0.01], band[weigth > 0.01] - band[0], p0=p0,
                               sigma=weigth[weigth > 0.01])

        # Standard error in the calculation of the mass
        # perr = np.sqrt(np.diag(pcov))
        masses.append(popt[0])
        fit = band[0] + parabolic(bands[0], popt[0])

        if plot_result:
            import matplotlib.pyplot as plt
            plt.plot(bands[0] / kmax, band / q, 'k')
            plt.plot(bands[0] / kmax, fit / q, 'r')

    if a0 < host_material.lattice_constant:
        # We need to figure out a safer way of doing this
        # print("CAUTION: BLINDLY SWAPPING HH,LH LABELS BASED ON LATTICE CONSTANT")
        masses[1], masses[2] = masses[2], masses[1]

    if plot_result:
        plt.xlabel('k (2*pi/a)')
        plt.ylabel('Energy (eV)')
        plt.ylim(-2, 2)
        plt.show()

    masses = abs(np.array(masses)) * m0
    return tuple(masses[::-1].tolist())


def parabolic(k, mr):
    return hbar ** 2 / (2 * m0 * mr) * k ** 2


def average_inplane_effective_mass(material, host_material, averaging_points=10, dk=0.026):
    """ Calculates the average in-plane effective mass for the four bands, assuming always the XY plane. """

    # p = phi = Polar angle, always constant for the XY plane
    # t = theta = Azimuth. Given the symmetry of the problem, we need to average only 45˚ of the k space
    p = np.pi / 2
    t = np.linspace(0, np.pi / 4, averaging_points)

    mc = 0
    mhh = 0
    mlh = 0
    mso = 0

    # We calculate the effective mass on each direction and average them
    for tt in t:
        bands = KPbands(material, host_material, p=p, t=tt)
        result = fit_effective_masses(bands, material, host_material, dk=dk)

        mc = mc + result[0]
        mhh = mhh + result[1]
        mlh = mlh + result[2]
        mso = mso + result[3]

    mc = mc / averaging_points
    mhh = mhh / averaging_points
    mlh = mlh / averaging_points
    mso = mso / averaging_points

    return mc, mhh, mlh, mso


def kp8x8_bulk(material, host_material, averaging_points=10, dk=0.026):
    """ Combines the above functions to provide a compact result of the band edges and the effective masses. """

    mc, mhh, mlh, mso = average_inplane_effective_mass(material, host_material, averaging_points=averaging_points,
                                                       dk=dk)
    c, hh, lh, so = KPbands(material, host_material, return_edges_only=True)

    return c, hh, lh, so, mc, mhh, mlh, mso


if __name__ == "__main__":
    import solcore

    # Material parameters
    GaAs = solcore.material("GaAs")(T=300)
    InGaAs2 = solcore.material("InGaAs")(In=0.28, T=300)
    InGaAs3 = solcore.material("GaInAs")(In=0.2, T=300)
    GaAsP = solcore.material("GaInP")(In=0.7, T=300)

    bands = kp_bands(GaAs, InGaAs3, kx=0, ky=0, kz=0, graph=True, fit_effective_mass=True, effective_mass_direction="L",
                     return_so=True)

    # bands = KPbands(GaAs, GaAsP, fraction=0.2, plot_result=False)
    # result = fit_effective_masses(bands, GaAs, GaAsP, plot_result=False)
    #
    # print(result)
