import numpy as np
import sys
from scipy.sparse import dia_matrix
from scipy.sparse.linalg import eigs
from solcore.constants import q, hbar, electron_mass as m0

from operator import itemgetter

sort_simultaneous = lambda *lists: [list(x) for x in zip(*sorted(zip(*lists), key=itemgetter(0)))]


def alfa(Eband, m):
    return (1 - m) ** 2 * q / Eband


def kp4x4(mat, host_material):
    """ Calculates the bandedges and bulk effective masses of a material under biaxial strain at the gamma point.
    SO band is ignored.

    This function provides the position at the zone center of the C, HH and LH bands of a material under biaxial strain.
    Both, the paralell and transverse effective masses (with respect the Z axis) are provided. SO band is ignored,
    therefore the masses are independent of strain. """

    Delta = mat.spin_orbit_splitting
    Eg = mat.band_gap
    Ev0 = mat.valence_band_offset
    Ec0 = mat.valence_band_offset + Eg
    a0 = mat.lattice_constant
    g1 = mat.gamma1
    g2 = mat.gamma2
    g3 = mat.gamma3
    ac = mat.a_c
    av = abs(mat.a_v)
    b = mat.b
    me = mat.eff_mass_electron_Gamma

    # Strain parameters
    exx = (host_material.lattice_constant - a0) / a0
    ezz = - 2 * mat.c12 / mat.c11 * exx

    Oe = ac * (exx + exx + ezz)
    Pe = -av * (exx + exx + ezz)
    Qe = -b / 2 * (exx + exx - 2 * ezz)

    Ec = Ec0 + Oe
    Evhh = Ev0 - Pe - Qe
    Evlh = Ev0 - Pe + Qe

    # Effective masses paralell to the Z axis
    mhh_p = 1. / (g1 - 2 * g2)
    mlh_p = 1. / (g1 + 2 * g2)

    # Effective masses transverse to the Z axis, that is, in the plane of the strain
    mhh_t = 1. / (g1 + g2)
    mlh_t = 1. / (g1 - g2)

    return Ec, Evhh, Evlh, me * m0, mhh_p * m0, mlh_p * m0, mhh_t * m0, mlh_t * m0


def kp6x6(mat, host_material):
    """ Calculates the bandedges and bulk effective masses of a material under biaxial strain at the gamma point
    including the SO band.

    This function provides the position at the zone center of the C, HH, LH and SO bands of a material under biaxial
    strain. Both, the paralell and transverse effective masses (with respect the Z axis) are provided."""

    Delta = mat.spin_orbit_splitting
    Eg = mat.band_gap
    Ev0 = mat.valence_band_offset
    Ec0 = mat.valence_band_offset + Eg
    a0 = mat.lattice_constant
    g1 = mat.gamma1
    g2 = mat.gamma2
    g3 = mat.gamma3
    ac = mat.a_c
    av = abs(mat.a_v)
    b = mat.b
    me = mat.eff_mass_electron_Gamma

    # Strain parameters
    exx = (host_material.lattice_constant - a0) / a0
    ezz = - 2 * mat.c12 / mat.c11 * exx

    Oe = ac * (exx + exx + ezz)
    Pe = -av * (exx + exx + ezz)
    Qe = -b / 2 * (exx + exx - 2 * ezz)

    x = Qe / Delta
    if x == 0: x = 1e-10
    com = lambda s: x - 1 + s * np.sqrt(1 + 2 * x + 9 * x ** 2)
    f = lambda s: (2 * x * (1 + 1.5 * com(s)) + 6 * x ** 2) / (
        0.75 * com(s) ** 2 + com(s) - 3 * x ** 2)  # s=+1 for LH and s=-1 for SO

    Ec = Ec0 + Oe
    Evhh = Ev0 - Pe - Qe
    Evlh = Ev0 - Pe + 0.5 * (Qe - Delta + np.sqrt(Delta ** 2 + 2 * Delta * Qe + 9 * Qe ** 2))
    Evso = Ev0 - Pe + 0.5 * (Qe - Delta - np.sqrt(Delta ** 2 + 2 * Delta * Qe + 9 * Qe ** 2))

    # Effective masses paralell to the Z axis
    mhh_p = 1. / (g1 - 2 * g2)
    mlh_p = 1. / (g1 + 2 * f(+1) * g2)
    mso_p = 1. / (g1 + 2 * f(-1) * g2)

    # Effective masses transverse to the Z axis, that is, in the plane of the strain
    mhh_t = 1. / (g1 + g2)
    mlh_t = 1. / (g1 - f(+1) * g2)
    mso_t = 1. / (g1 - f(-1) * g2)

    return Ec, Evhh, Evlh, Evso, me * m0, mhh_p * m0, mlh_p * m0, mso_p * m0, mhh_t * m0, mlh_t * m0, mso_t * m0


def fill_hamiltonian_holes_4x4(A1, A2, B1, B2, C, D, delta, block):
    # Number of meshpoints in z direction. The number of equations to solve is 2N
    N = len(A1)

    # Values of the diagonals, depending on the equation: 1, 2, 3 or 4
    ll14 = lambda i: (A2[i + 1] - 4 * A2[i] - A2[i - 1]) / (4 * delta ** 2)
    l13 = lambda i: D[i] / (2 * delta)
    d14 = lambda i: 2 * A2[i] / delta ** 2 + A1[i]
    u13 = lambda i: C[i] - (D[i + 1] - D[i - 1]) / (4 * delta)
    uu14 = lambda i: - (A2[i + 1] + 4 * A2[i] - A2[i - 1]) / (4 * delta ** 2)
    uuu13 = lambda i: - D[i] / (2 * delta)

    lll24 = lambda i: - D[i] / (2 * delta)
    ll23 = lambda i: (B2[i + 1] - 4 * B2[i] - B2[i - 1]) / (4 * delta ** 2)
    l24 = lambda i: C[i] + (D[i + 1] - D[i - 1]) / (4 * delta)
    d23 = lambda i: 2 * B2[i] / delta ** 2 + B1[i]
    u24 = lambda i: D[i] / (2 * delta)
    uu23 = lambda i: - (B2[i + 1] + 4 * B2[i] - B2[i - 1]) / (4 * delta ** 2)

    lll = np.zeros(2 * N)
    ll = np.zeros(2 * N)
    l = np.zeros(2 * N)
    d = np.zeros(2 * N)
    u = np.zeros(2 * N)
    uu = np.zeros(2 * N)
    uuu = np.zeros(2 * N)

    if block is "U":
        # We write the equations for the UPPER 2x2 block of the Hamiltonian, for g1 and g2

        # We fill the first two equations separatelly
        d[0] = d14(0)
        u[0] = C[0] - (D[1] - D[0]) / (4 * delta)
        uu[0] = - (A2[1] + 3 * A2[0]) / (4 * delta ** 2)
        uuu[0] = uuu13(0)

        l[1] = C[0] + (D[1] - D[0]) / (4 * delta)
        d[1] = d23(0)
        u[1] = u24(0)
        uu[1] = - (B2[1] + 3 * B2[0]) / (4 * delta ** 2)

        for i in range(1, N - 1):
            # The first equation, for g1
            k = 2 * i

            ll[k] = ll14(i)
            l[k] = l13(i)
            d[k] = d14(i)
            u[k] = u13(i)
            uu[k] = uu14(i)
            uuu[k] = uuu13(i)

            # The second equation, for g2
            k = k + 1

            lll[k] = lll24(i)
            ll[k] = ll23(i)
            l[k] = l24(i)
            d[k] = d23(i)
            u[k] = u24(i)
            uu[k] = uu23(i)

        # We fill the last two equations separatelly
        ll[2 * (N - 1)] = (- 3 * A2[N - 1] - A2[N - 2]) / (4 * delta ** 2)
        l[2 * (N - 1)] = l13(N - 1)
        d[2 * (N - 1)] = d14(N - 1)
        u[2 * (N - 1)] = C[N - 1] - (D[N - 1] - D[N - 2]) / (4 * delta)

        lll[2 * (N - 1) + 1] = lll24(N - 1)
        ll[2 * (N - 1) + 1] = (- 3 * B2[N - 1] - B2[N - 2]) / (4 * delta ** 2)
        l[2 * (N - 1) + 1] = C[N - 1] + (D[N - 1] - D[N - 2]) / (4 * delta)
        d[2 * (N - 1) + 1] = d23(N - 1)

    else:
        # We write the equations for the LOWER 2x2 block of the Hamiltonian, for g3 and g4

        # We fill the first two equations separatelly
        d[0] = d23(0)
        u[0] = C[0] - (D[1] - D[0]) / (4 * delta)
        uu[0] = - (B2[1] + 3 * B2[0]) / (4 * delta ** 2)
        uuu[0] = uuu13(0)

        l[1] = C[0] + (D[1] - D[0]) / (4 * delta)
        d[1] = d14(0)
        u[1] = u24(0)
        uu[1] = - (A2[1] + 3 * A2[0]) / (4 * delta ** 2)

        for i in range(1, N - 1):
            # The first equation, for g3
            k = 2 * i

            ll[k] = ll23(i)
            l[k] = l13(i)
            d[k] = d23(i)
            u[k] = u13(i)
            uu[k] = uu23(i)
            uuu[k] = uuu13(i)

            # The second equation, for g4
            k = k + 1

            lll[k] = lll24(i)
            ll[k] = ll14(i)
            l[k] = l24(i)
            d[k] = d14(i)
            u[k] = u24(i)
            uu[k] = uu14(i)

        # We fill the last two equations separatelly

        ll[2 * (N - 1)] = (- 3 * B2[N - 1] - B2[N - 2]) / (4 * delta ** 2)
        l[2 * (N - 1)] = l13(N - 1)
        d[2 * (N - 1)] = d23(N - 1)
        u[2 * (N - 1)] = C[N - 1] - (D[N - 1] - D[N - 2]) / (4 * delta)

        lll[2 * (N - 1) + 1] = lll24(N - 1)
        ll[2 * (N - 1) + 1] = (- 3 * A2[N - 1] - A2[N - 2]) / (4 * delta ** 2)
        l[2 * (N - 1) + 1] = C[N - 1] + (D[N - 1] - D[N - 2]) / (4 * delta)
        d[2 * (N - 1) + 1] = d14(N - 1)

    lll = np.roll(lll, -3)
    ll = np.roll(ll, -2)
    l = np.roll(l, -1)
    u = np.roll(u, 1)
    uu = np.roll(uu, 2)
    uuu = np.roll(uuu, 3)

    offsets = [-3, -2, -1, 0, 1, 2, 3]
    diagonals = np.array([lll, ll, l, d, u, uu, uuu])

    output = dia_matrix((diagonals, offsets), shape=(2 * N, 2 * N))

    # View hamiltonian, good for debugging
    # import pylab
    # pylab.imshow(np.real(output.todense()))
    # pylab.show()

    return output


def solve_holes_QW_at_kt_4x4(kt, z, fhh, flh, g1, g2, g3, num=(10, 10), quasiconfined=0.0, symmetric=False):
    # Normalization factors
    N = len(z)
    L = max(z)
    delta = (z[1] - z[0]) / L
    kt = kt * L
    E0 = hbar ** 2 / (m0 * L ** 2)

    # We invert and shift the potential so that it is possitive with the minimum value around zero
    offset = min(np.amin(fhh), np.amin(flh))
    fhh = (offset - fhh) / E0
    flh = (offset - flh) / E0

    # Fill the components of the equations based on the effective masses, potential and k value
    A1 = 0.5 * (g1 + g2) * kt ** 2 + fhh
    B1 = 0.5 * (g1 - g2) * kt ** 2 + flh
    A2 = 0.5 * (g1 - 2 * g2)
    B2 = 0.5 * (g1 + 2 * g2)
    C = np.sqrt(3) / 2 * 0.5 * (g2 + g3) * kt ** 2
    D = np.sqrt(3) * g3 * kt

    # Using the shift-invert mode will help to speed up the algorithm by calculating the eigenvalues only inside the QW or very close
    # See explanation: https://docs.scipy.org/doc/scipy-0.15.1/reference/tutorial/arpack.html
    sigma = min(np.amin(fhh), np.amin(
        flh))  # The bottom of the potential (valence band edge) is the target energy for the eigenvalues

    # Allow for quasiconfined levels to go through. They can be discarded later with the filter, if necessary
    potmax = max(np.amax(fhh), np.amax(flh)) + quasiconfined * q / E0

    # Solutions for upper hamiltonian
    H_U = fill_hamiltonian_holes_4x4(A1, A2, B1, B2, C, D, delta, "U")
    E_U, Psi_U = eigs(H_U, k=num[0], which='LR', sigma=sigma)

    Psi_g1 = Psi_U[range(0, len(Psi_U), 2)]
    Psi_g2 = Psi_U[range(1, len(Psi_U), 2)]

    if quasiconfined < 0.0:
        potmax = max(E_U)

    confined_levels = [i for i, e in enumerate(E_U) if (e < potmax and i < num[1])]
    E_U = E_U[confined_levels].real
    idx_sorted = np.argsort(E_U)
    E_U = E_U[idx_sorted]
    Psi_g1 = np.array(Psi_g1[:, confined_levels]).transpose().real[idx_sorted]
    Psi_g2 = np.array(Psi_g2[:, confined_levels]).transpose().real[idx_sorted]

    norm1 = [(np.sqrt(np.trapz(p * p, x=z))) for p in Psi_g1]
    norm2 = [(np.sqrt(np.trapz(p * p, x=z))) for p in Psi_g2]
    # Normalise the wavefunction except if the norm is very small. That's the case for some E at k near 0 only
    for i in range(len(Psi_g1)):
        if norm1[i] > 0.01 * max(norm1):
            Psi_g1[i] = Psi_g1[i] / norm1[i]
        else:
            Psi_g1[i] = Psi_g1[i] * 0.0

        if norm2[i] > 0.01 * max(norm1):
            Psi_g2[i] = Psi_g2[i] / norm2[i]
        else:
            Psi_g2[i] = Psi_g2[i] * 0.0

    E_U = offset - E_U * E0

    # If the structure is not symmetric, we have to solve the lower hamiltonian, too
    if not symmetric:
        H_L = fill_hamiltonian_holes_4x4(A1, A2, B1, B2, C, D, delta, "L")
        E_L, Psi_L = eigs(H_L, k=num[0], which='LR', sigma=sigma)

        Psi_g3 = Psi_L[range(0, len(Psi_L), 2)]
        Psi_g4 = Psi_L[range(1, len(Psi_L), 2)]

        if quasiconfined < 0.0:
            potmax = max(E_L)

        confined_levels = [i for i, e in enumerate(E_L) if (e < potmax and i < num[1])]
        E_L = E_L[confined_levels].real
        Psi_g3 = np.array(Psi_g3[:, confined_levels]).transpose().real
        Psi_g4 = np.array(Psi_g4[:, confined_levels]).transpose().real

        norm3 = [(np.sqrt(np.trapz(p * p, x=z))) for p in Psi_g3]
        norm4 = [(np.sqrt(np.trapz(p * p, x=z))) for p in Psi_g4]
        for i in range(len(Psi_g3)):
            if norm3[i] > 0.01 * max(norm3):
                Psi_g3[i] = Psi_g3[i] / norm3[i]
            else:
                Psi_g3[i] = Psi_g3[i] * 0.0

            if norm4[i] > 0.01 * max(norm4):
                Psi_g4[i] = Psi_g4[i] / norm4[i]
            else:
                Psi_g4[i] = Psi_g4[i] * 0.0

        E_L = offset - E_L * E0

        return {"E_U": E_U,
                "E_L": E_L,
                "Psi_g1": Psi_g1,
                "Psi_g2": Psi_g2,
                "Psi_g3": Psi_g3,
                "Psi_g4": Psi_g4}

    else:
        # If the structure is symmetric, we simply duplicate the energies and reverse the wavefunctions
        return {"E_U": E_U,
                "Psi_g1": Psi_g1,
                "Psi_g2": Psi_g2}


def fill_hamiltonian_holes_6x6(A1, A2, B1, B2, C, D, delta, block):
    print('Fully coupled 6x6 kp solver not implemented yet. Only the band edges can be calculated. ')
    sys.exit(-1)


def solve_holes_QW_at_kt_6x6(kt, z, fhh, flh, g1, g2, g3, num=(10, 10), quasiconfined=0.0, symmetric=False):
    print('Fully coupled 6x6 kp solver not implemented yet. Only the band edges can be calculated. ')
    sys.exit(-1)


def solve_electrons_QW_at_kt_parabolic(kt, z, fe, me, num=(10, 10), quasiconfined=0.0):
    """ Returns eignvalue and eigenvectors of an arbitrary potential.

    A tridiagonal matrix is constructed by writing the variable effective
    mass Schrodinger equation over a series of mesh points. The eigenvalues
    of the matrix correspond to the allowed energy levels of the system.

    The previous solver, eig, has been replaced by the spare matrix version, eigs, that is faster to compute


    "Varible effective mass Schordinger equation and tridiagonal solution method.",
    "Frensley, W. R. (1991). \
    Numerical evaluation of resonant states. \
    Superlattices and Microstructures, 11(3), 347350. \
    doi:10.1016/0749-6036(92)90396-M")

    """

    # Normalization factors
    L = max(z)
    dz = np.gradient(z / L)
    kt = kt * L
    E0 = hbar ** 2 / (m0 * L ** 2)

    # We shift the potential so that it is possitive with the minimum value around zero
    offset = np.amin(fe)
    fe = (fe - offset) / E0

    N = len(fe)
    m = me / m0

    # Vectorise effective mass differences to avoid a loop
    m = np.insert(m, (0, len(m)), (m[0], m[-1]))
    m_a = m[0:-2]  # m_(j-1) from Frensley, W. R. (1991)
    m_b = m[1:-1]  # m_j
    m_c = m[2:]  # m_(j+1)

    # These are the interior diagonals of equation 18 in the above ref.
    axis = 0.25 / dz ** 2 * (
    1 / m_a + 2 / m_b + 1 / m_c) + fe + 0.5 * kt ** 2 / m_b  # d_j from Frensley, W. R. (1991), modified to include kt and normalised variables
    upper = 0.25 / dz ** 2 * (1 / m_a + 1 / m_b)  # s_(j+1)
    lower = 0.25 / dz ** 2 * (1 / m_b + 1 / m_c)  # s_j

    index = (-1, 0, 1)
    diagonals = (-lower, axis, -upper)

    H = dia_matrix((diagonals, index), shape=(N, N))

    sigma = 0.0  # The bottom of the potential (conduction band) is the target energy for the eigenvalues

    # The heavy numerical calculation
    E, Psi = eigs(H, k=num[0], which='LR', sigma=sigma)

    # Allow for quasi confined levels to go through.
    potmax = max(fe) + quasiconfined * q / E0
    if quasiconfined < 0.0:
        potmax = max(E)

    confined_levels = [i for i, e in enumerate(E) if (e < potmax and i < num[1])]
    E = E[confined_levels].real
    Psi = np.array(Psi[:, confined_levels]).transpose()
    Psi = [(p / np.sqrt(np.trapz(p * p.conjugate(), x=z))) for p in Psi]

    E = E * E0 + offset

    # The solutions are degenerate: same energy
    return {"E_e": E,
            "Psi_e": Psi,
            "z": z}


def band_sorting(bands, symmetric=False):
    """ Sort the bands in C, HH and LH according to their character at k=0
    :param bands: A dictionary containing the bandstructure as calculated by solve_electrons_QW_at_kt_parabolic and
            solve_holes_QW_at_kt_4x4
    :return: A dictionary with the same input information but arranged in a different way, including labels
    """

    print("Organising the bands into E, HH and LH")

    Ee = []
    Ehh = []
    Elh = []
    psi_e = []
    psi_hh = []
    psi_lh = []

    # We vectorise the energies and the wavefunctions
    num_e = len(bands[0][1]["E_e"])
    num_hu = len(bands[0][1]["E_U"])
    kt = np.zeros(len(bands))
    electrons = np.zeros((len(bands), num_e))
    holes_U = np.zeros((len(bands), num_hu))
    Psi_e = np.zeros((len(bands), num_e, len(bands[0][1]["z"])))
    Psi_g1 = np.zeros((len(bands), num_hu, len(bands[0][1]["z"])))
    Psi_g2 = np.zeros((len(bands), num_hu, len(bands[0][1]["z"])))

    for i in range(len(bands)):
        kt[i] = bands[i][0]

        # Electrons
        for j in range(num_e):
            try:
                electrons[i, j] = bands[i][1]["E_e"][j]
                Psi_e[i, j, :] = bands[i][1]["Psi_e"][j].real
            except IndexError:
                electrons[i, j] = np.inf
        # Holes
        for j in range(num_hu):
            try:
                holes_U[i, j] = bands[i][1]["E_U"][j]
                Psi_g1[i, j, :] = bands[i][1]["Psi_g1"][j].real
                Psi_g2[i, j, :] = bands[i][1]["Psi_g2"][j].real
            except IndexError:
                holes_U[i, j] = np.inf

    # Conduction bands are easy to organise... the are all the same!!
    # For the hole bands, we say it is HH if the norm of Psi_g1 at k=0 is 1, and LH otherwise
    # This sorting is not very reliable as just outside k=0 the band can be strongly (mostly) LH-like
    for i in range(num_e):
        Ee.append(electrons[:, i])
        psi_e.append(Psi_e[:, i, :])

    norm = [(np.sqrt(np.trapz(p * p, x=bands[0][1]["z"]))) for p in Psi_g1[0, :]]
    for i in range(num_hu):
        if norm[i] > 0.5:
            Ehh.append(holes_U[:, i])
            psi_hh.append(Psi_g1[:, i, :])
        else:
            Elh.append(holes_U[:, i])
            psi_lh.append(Psi_g2[:, i, :])

    # We repeat the process if we are using the non-symmetric mode
    if not symmetric:
        num_hl = len(bands[0][1]["E_L"])
        holes_L = np.zeros((len(bands), num_hl))
        Psi_g3 = np.zeros((len(bands), num_hl, len(bands[0][1]["z"])))
        Psi_g4 = np.zeros((len(bands), num_hl, len(bands[0][1]["z"])))

        for i in range(len(bands)):

            for j in range(num_hl):
                try:
                    holes_L[i, j] = bands[i][1]["E_L"][j]
                    Psi_g3[i, j, :] = bands[i][1]["Psi_g3"][j].real
                    Psi_g4[i, j, :] = bands[i][1]["Psi_g4"][j].real
                except IndexError:
                    holes_L[i, j] = np.inf

                    # Plot the wavefunctions vs kt and z. Good for debugging
                    # import matplotlib.pyplot as plt
                    # import sys
                    # K, Z = np.meshgrid(kt, bands[0][1]["z"])
                    # plt.contourf(Z, K, Psi_e[:, 0, :].T**2, 50)
                    # plt.show()
                    # sys.exit()

        norm = [(np.sqrt(np.trapz(p * p, x=bands[0][1]["z"]))) for p in Psi_g4[0, :]]
        for i in range(len(norm)):
            if norm[i] > 0.5:
                Ehh.append(holes_L[:, i])
                psi_hh.append(Psi_g4[:, i, :])
            else:
                Elh.append(holes_L[:, i])
                psi_lh.append(Psi_g3[:, i, :])

        Ehh = np.array(Ehh)
        Elh = np.array(Elh)
        psi_hh = np.array(psi_hh)
        psi_lh = np.array(psi_lh)

    # Finally, we convert the lists containing the useful stuff into arrays
    Ee = np.array(Ee)
    Ehh = np.array(Ehh)
    Elh = np.array(Elh)
    psi_e = np.array(psi_e)
    psi_hh = np.array(psi_hh)
    psi_lh = np.array(psi_lh)

    return {"kt": kt,           # k points
            "Ee": Ee,           # Electron energy levels vs K
            "psi_e": psi_e,     # Electron wavefunctions vs X and K
            "Ehh": Ehh,         # Heavy hole energy levels vs K
            "psi_hh": psi_hh,   # Heavy hole wavefunctions vs X and K
            "Elh": Elh,         # Light hole energy levels vs K
            "psi_lh": psi_lh,   # Light hole wavefunctions vs X and K
            "symmetric": symmetric}  # If the structure is symmetric


def solve_bandstructure_QW(structure, num=10, kpoints=40, krange=5e9, quasiconfined=0.0, symmetric=False,
                           plot_bands=False):
    allk_t = np.linspace(0, krange, kpoints)
    bands = []

    z = structure["x"]
    fe = structure["Ve"]
    me = structure["me"]
    fhh = structure["Vhh"]
    flh = structure["Vlh"]
    g1 = structure["g1"]
    g2 = structure["g2"]
    g3 = structure["g3"]

    # First we solve at kt = 0
    new_kpoint = solve_electrons_QW_at_kt_parabolic(0.0, z, fe, me, num=(num, num), quasiconfined=quasiconfined)
    new_kpoint.update(solve_holes_QW_at_kt_4x4(0.0, z, fhh, flh, g1, g2, g3, num=(num, num), symmetric=symmetric,
                                               quasiconfined=quasiconfined))
    bands.append([0.0, new_kpoint])

    # We search only for the solutions corresponding to energy levels confined at k=0
    num_e = len(new_kpoint['E_e'])
    num_h = len(new_kpoint['E_U'])

    for kt in allk_t[1:]:
        new_kpoint = solve_electrons_QW_at_kt_parabolic(kt, z, fe, me, num=(num, num_e), quasiconfined=-0.1)
        new_kpoint.update(solve_holes_QW_at_kt_4x4(kt, z, fhh, flh, g1, g2, g3, num=(num, num_h), symmetric=symmetric,
                                                   quasiconfined=0.1))

        bands.append([kt, new_kpoint])

    print("Calculation finished!!")

    bands = band_sorting(bands, symmetric=symmetric)

    if plot_bands:
        import matplotlib.pyplot as plt

        fig, (ax1, ax2) = plt.subplots(2, 1)

        # Electrons
        ax1.plot(bands["kt"]*1e-9, bands["Ee"][0] / q, 'b', label='Ee')
        for i in range(1, len(bands["Ee"])):
            ax1.plot(bands["kt"]*1e-9, bands["Ee"][i] / q, 'b')

        # HH
        ax2.plot(bands["kt"]*1e-9, bands["Ehh"][0] / q, 'r', label='Ehh')
        for i in range(1, len(bands["Ehh"])):
            ax2.plot(bands["kt"]*1e-9, bands["Ehh"][i] / q, 'r')

        # LH
        ax2.plot(bands["kt"]*1e-9, bands["Elh"][0] / q, 'g', label='Elh')
        for i in range(1, len(bands["Elh"])):
            ax2.plot(bands["kt"]*1e-9, bands["Elh"][i] / q, 'g')

        ax1.set_ylabel('Energy (eV)')
        ax2.set_ylabel('Energy (eV)')
        plt.xlabel('k (nm$^{-1}$)')
        plt.xlim(0, max(bands["kt"]*1e-9))
        plt.setp(ax1.get_xticklabels(), visible=False)
        plt.legend()
        plt.tight_layout()
        plt.show()

    bands2 = structure.copy()
    bands2.update(bands)

    return bands2
