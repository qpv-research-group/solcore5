from numpy import ones, array, sqrt, trapz
from scipy.sparse import dia_matrix
from scipy.sparse.linalg import eigs
from solcore.science_tracker import science_reference
from solcore.constants import *
import numpy as np

from operator import itemgetter

sort_simultaneous = lambda *lists: [list(x) for x in zip(*sorted(zip(*lists), key=itemgetter(0)))]

from . import structure_utilities
from .graphics import Graph, GraphData


def tridiag_euler(V, z, m, periodic=False, num_eigenvalues=10, quasiconfined=0):
    """
    Returns eignvalue and eigenvectors of an arbitrary potential.
    
    A tridiagonal matrix is constructed by writing the variable effective
    mass Schrodinger equation over a series of mesh points. The eigenvalues
    of the matrix correspond to the allowed energy levels of the system.
    
    The previous solver, eig, has been replaced by the spare matrix version, eigs, that is faster to compute
    """

    science_reference("Varible effective mass Schordinger equation and tridiagonal solution method.",
                      "Frensley, W. R. (1991). \
                      Numerical evaluation of resonant states. \
                      Superlattices and Microstructures, 11(3), 347350. \
                      doi:10.1016/0749-6036(92)90396-M")

    N = len(V)
    dz = np.gradient(z)
    m = m * ones(N)

    # Vectorise effective mass differences to avoid a loop
    m = np.insert(m, (0, len(m)), (m[0], m[-1]))
    m_a = m[0:-2]  # m_(j-1) from Frensley, W. R. (1991)
    m_b = m[1:-1]  # m_j
    m_c = m[2:]  # m_(j+1)

    # These are the interior diagonals of equation 18 in the above ref.
    axis = hbar ** 2 / (4 * dz ** 2) * (1 / m_a + 2 / m_b + 1 / m_c) + V  # d_j from Frensley, W. R. (1991)
    upper = hbar ** 2 / (4 * dz ** 2) * (1 / m_a + 1 / m_b)  # s_(j+1)
    lower = hbar ** 2 / (4 * dz ** 2) * (1 / m_b + 1 / m_c)  # s_j

    if periodic:
        TopRight = np.zeros(len(axis))
        BottomLeft = np.zeros(len(axis))
        TopRight[-1] = -lower[-1]  # The last point becomes the one before the first one
        BottomLeft[0] = -upper[0]  # The first point becomes the one after the last one

        index = (-N + 1, -1, 0, 1, N - 1)
        diagonals = (BottomLeft, -lower, axis, -upper, TopRight)

    else:
        index = (-1, 0, 1)
        diagonals = (-lower, axis, -upper)

    H = dia_matrix((diagonals, index), shape=(N, N))

    # H[0,-1] = -lower[-1]    # The last point becomes the one before the first one
    # H[-1,0] = -upper[0]     # The first point becomes the one after the last one

    sigma = np.min(V)  # The top of the potential (valence band) is the target energy for the eigenvalues

    # The heavy numerical calculation
    E, Psi = eigs(H, k=num_eigenvalues, which='LR', sigma=sigma)

    # Allow for quasi confined levels to go through. They can be discarded later with the filter
    confined_levels = [i for i, e in enumerate(E) if e < quasiconfined * q]
    E, Psi = E[confined_levels].real, array(Psi[:, confined_levels]).transpose()
    Psi = [(p / sqrt(trapz(p * p, x=z))).real for p in Psi]

    E, Psi = sort_simultaneous(E, Psi)

    return E, Psi


def schroedinger_solve(x, V, m, num_eigenvalues=10, periodic=False, offset=0, electron=True, quasiconfined=0):
    """Returns normalised wavefuctions from the potential profile.

    Arguments:

        x -- spatial grid
        V -- potential
        m -- effective mass

    Keywords:

        electron -- whether the wavefunctions describe electrons or holes (default: True)
        num_eigenvalues -- Number of eigenvalues to calculate (default = 10)
        periodic   -- not to sure what this does (default: False)
    """

    if electron:
        v_shift = max(V) + offset
        E, psi = tridiag_euler(V - v_shift, x, m, num_eigenvalues=num_eigenvalues, periodic=periodic,
                               quasiconfined=quasiconfined)
        E = list() if len(E) == 0 else np.array(E) + v_shift

    else:
        v_shift = min(V) - offset
        E, psi = tridiag_euler(-V + v_shift, x, m, num_eigenvalues=num_eigenvalues, periodic=periodic,
                               quasiconfined=quasiconfined)
        E = list() if len(E) == 0 else -np.array(E) + v_shift

    return E, psi


def __potentials_to_wavefunctions_energies_internal(x, Ve, me, Vhh, mhh, Vlh, mlh, num_eigenvalues=10, periodic=False,
                                                    offset=0, filter_strength=0, structure=None, quasiconfined=0):
    """Returns normalised wavefuctions from the potential profile.

    Arguments:

        x -- spatial grid
        Ve -- electron potential
        me -- electron effective mass
        Vlh -- light hole potential
        mlh -- light hole effective mass
        Vhh -- heavy hole potential
        mhh -- heavy hole effective mass


    Keywords:

        num_eigenvalues -- Number of eigenvalues to calculate (default = 10)
        periodic   -- not to sure what this does (default: False)
    """

    Ee, psi_e = schroedinger_solve(x, Ve, me, num_eigenvalues, periodic, offset, electron=True,
                                   quasiconfined=quasiconfined)
    Ehh, psi_hh = schroedinger_solve(x, Vhh, mhh, num_eigenvalues, periodic, offset, electron=False,
                                     quasiconfined=quasiconfined)
    Elh, psi_lh = schroedinger_solve(x, Vlh, mlh, num_eigenvalues, periodic, offset, electron=False,
                                     quasiconfined=quasiconfined)

    if filter_strength != 0:
        assert structure is not None, "Need to provide structure to find well regions for filtering"
        print("filtering")
        Ee, psi_e = discard_unconfined(x, structure, Ee, psi_e, filter_strength)
        Ehh, psi_hh = discard_unconfined(x, structure, Ehh, psi_hh, filter_strength)
        Elh, psi_lh = discard_unconfined(x, structure, Elh, psi_lh, filter_strength)

    # Ee, psi_e = discard_unconfined_energy(Ee, psi_e, Ve, quasiconfined)
    # Ehh, psi_hh = discard_unconfined_energy(Ehh, psi_hh, Vhh, quasiconfined)
    # Elh, psi_lh = discard_unconfined_energy(Elh, psi_lh, Vlh, quasiconfined)

    # print(me/electron_mass, mhh/electron_mass, mlh/electron_mass)
    # print(me_plane, mhh_plane, mlh_plane)
    #
    # import sys
    # sys.exit()

    return x, Ee, psi_e, Ehh, psi_hh, Elh, psi_lh


def discard_unconfined_energy(E, psi, V, quasiconfined):
    before = len(E)
    maxE = min(abs(V[0]), max(abs(V))) + quasiconfined
    try:
        E, psi = zip(*[(E_i, psi_i)
                       for (E_i, psi_i) in zip(E, psi)
                       if abs(E_i) <= maxE])
    except ValueError as exception:
        print("Warning: wavefunction filter removed all states for this band, try reducing the filter strength.")
        return ([], [])
    print("Wavefunction filter removed %d state%s." % (before - len(E), "" if (before - len(E)) == 1 else "s"))
    return E, psi


def discard_unconfined(x, structure, E, psi, threshold=0.8):
    if threshold == 0:  # bypass the filter code, saving time
        return E, psi

    indx = structure_utilities.well_regions(x, structure)

    before = len(E)
    # for (E_i, psi_i) in zip(E, psi): 
    #     print (psi_i[indx])

    try:

        E, psi = zip(*[(E_i, psi_i)
                       for (E_i, psi_i) in zip(E, psi)
                       if np.trapz(psi_i[indx] ** 2, x=x[indx]) / np.trapz(psi_i ** 2, x=x) >= threshold])
    except ValueError as exception:
        print("Warning: wavefunction filter removed all states for this band, try reducing the filter strength.")
        return ([], [])
    print("Wavefunction filter removed %d state%s." % (before - len(E), "" if (before - len(E)) == 1 else "s"))
    return E, psi


def potentials_to_wavefunctions_energies(x, Ve, me, Vhh, mhh, Vlh, mlh, num_eigenvalues=10, periodic=False, offset=0,
                                         filter_strength=0, structure=None, quasiconfined=0, **kwargs):
    x, Ee, psi_e, Ehh, psi_hh, Elh, psi_lh = __potentials_to_wavefunctions_energies_internal(
        x, Ve, me, Vhh, mhh, Vlh, mlh, num_eigenvalues,
        periodic, offset,
        filter_strength, structure,
        quasiconfined)
    return {
        "x": x,
        "Ee": Ee,
        "psi_e": psi_e,
        "Ehh": Ehh,
        "psi_hh": psi_hh,
        "Elh": Elh,
        "psi_lh": psi_lh,
        "Ve": Ve,
        "me": me,
        "Vhh": Vhh,
        "Vlh": Vlh,
        "mhh": mhh,
        "mlh": mlh
    }


def graph(x, Ve, me, Vhh, mhh, Vlh, mlh, **kwargs):
    defaults = {
        "edit": lambda x, y: (x * 1e9, y / q),
        "xlabel": "Depth (nm)",
        "ylabel": "Energy (eV)",
    }
    defaults.update(kwargs)
    data = []

    normalise_psi = lambda p: p * q * 5e-6

    data.append(GraphData(x, Ve, color="black", linewidth=2))
    data.append(GraphData(x, Vlh, color="black", linewidth=2, dashes=[1, 1]))
    data.append(GraphData(x, Vhh, color="black", linewidth=2))

    g = Graph(data, **defaults)
    return g
