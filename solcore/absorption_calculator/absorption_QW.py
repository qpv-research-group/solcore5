from solcore.constants import *
from numpy import array
import numpy as np
from solcore.science_tracker import science_reference

H = lambda x: 0 if x < 0 else 1  ### blindly changed this !!
D = lambda x, width: 1 if (x >= -width) and (x <= width) else 0  # Delta function
L = lambda x, centre, hwhm: 1 / pi * (0.5 * hwhm) / (
    (x - centre) ** 2 + (0.5 * hwhm) ** 2)  # Lorenzian (area normalised to 1)
Gauss = lambda x, centre, hwhm: 1 / np.sqrt(2 * pi) / (0.5 * hwhm) * np.exp(
    -0.5 * (x - centre) ** 2 / (0.5 * hwhm) ** 2)  # Gaussian (area normalised to 1)


def exciton_rydberg_energy_2d(me, mh, eps_r):
    """
    :param me: electron effective mass (units: kg)
    :param mh: hole effective mass (units: kg)
    :param eps_r: dielectic constant (units: SI)
    :return: The exciton Rydberg energy
    """

    science_reference("Definition of the exciton Rydberg Energy for a quantum well.",
                      "S. L. Chuang, Physics of Optoelectonic Devices, Second Edition, p.554, Table 13.1")

    mr = me * mh / (me + mh)  # reduced effective mass
    if True:
        return mr * q ** 4 / (2 * hbar ** 2 * (4 * np.pi * eps_r * vacuum_permittivity) ** 2)
    Ry = 13.6 * 1.6e-19
    m0 = electron_mass
    Ry_eff = (mr / m0) * Ry / eps_r ** 2  # Effective Rydberg
    return Ry_eff


def exciton_bohr_radius(me, mh, eps):
    """
    :param me: electron effective mass (units: kg)
    :param mh: hole effective mass (units: kg)
    :param eps: dielectic constant (units: SI)
    :return: Exciton Borh radius
    """

    science_reference("Definition of the exciton bohr radius for a quantum well.",
                      "S. L. Chuang, Physics of Optoelectonic Devices, Second Edition, p.554, Table 13.1")
    mr = me * mh / (me + mh)  # reduced effective mass
    return hbar ** 2 / mr * (4 * np.pi * eps / q ** 2)


def alpha_c_hh_TE(E, z, E_e, E_hh, psi_e, psi_hh, well_width, me, mhh, Ep, nr):
    """     Absortion coefficient for incident light forming a transision between hh and c band of a quantum well

    NB. Assumes that valence band is a zero energy, might need to manually apply an offset.

    :param E: photon energy (units: J)
    :param z: mesh points along growth direction (units: m)
    :param E_e: Electron state energy (units: J)
    :param E_hh: Heavy hole state energy (units: J)
    :param psi_e: Electron envelope function (psi_e^2 must be normalised)
    :param psi_hh: Heavy hole envelope function (psi_hh^2 must be normalised)
    :param well_width: (units: m)
    :param me: electron effective mass (units: kg)
    :param mhh: heavy hole effective mass (units: kg)
    :param Ep: the Kane parameter "Optical dipole matrix elemet", sometimes "Momentum matrix element", e.g. Ep for GaAs ~ 28eV
    :param nr: refractive index
    :return:
    """

    m0 = electron_mass
    C0 = np.pi * q ** 2 / (nr * c * vacuum_permittivity * m0 ** 2 * E / hbar)
    mr = me * mhh / (me + mhh)  # reduced effective mass
    DOS_2D = mr / (np.pi * hbar ** 2 * well_width)
    Mbsq_3D = m0 / 6 * Ep
    Mbsq_2D = 3 / 2 * Mbsq_3D
    Ieh = np.trapz(psi_e * psi_hh, z) ** 2
    return C0 * Ieh * Mbsq_2D * DOS_2D * H(E - (E_e - E_hh))


def alpha_c_lh_TE(E, z, E_e, E_lh, psi_e, psi_lh, well_width, me, mlh, Ep, nr):
    """     Absortion coefficient for incident light forming a transision between hh and c band of a quantum well
    
    NB. Assumes that valence band is a zero energy, might need to manually apply an offset.
    
    :param E: photon energy (units: J)
    :param z: mesh points along growth direction (units: m)
    :param E_e: Electron state energy (units: J)
    :param E_hh: Heavy hole state energy (units: J)
    :param psi_e: Electron envelope function (psi_e^2 must be normalised)
    :param psi_hh: Heavy hole envelope function (psi_hh^2 must be normalised)
    :param well_width: (units: m)
    :param me: electron effective mass (units: kg)
    :param mhh: heavy hole effective mass (units: kg)
    :param Ep: the Kane parameter "Optical dipole matrix elemet", sometimes "Momentum matrix element", e.g. Ep for GaAs ~ 28eV
    :param nr: refractive index
    :return:
    """
    m0 = electron_mass
    C0 = np.pi * q ** 2 / (nr * c * vacuum_permittivity * m0 ** 2 * E / hbar)
    mr = me * mlh / (me + mlh)  # reduced effective mass
    DOS_2D = mr / (np.pi * hbar ** 2 * well_width)
    Mbsq_3D = m0 / 6 * Ep
    Mbsq_2D = 1 / 2 * Mbsq_3D
    Ieh = np.trapz(psi_e * psi_lh, z) ** 2
    return C0 * Ieh * Mbsq_2D * DOS_2D * H(E - (E_e - E_lh))


def alpha_exciton_ehh_TE(exciton_index, E, z, E_e, E_hh, psi_e, psi_hh, well_width, me, mhh, Ep, nr, eps,
                         hwhm=6e-3 * 1.6e-19, dimensionality=0.15, line_shape="Lorenzian"):
    """

    :param exciton_index:
    :param E:
    :param z:
    :param E_e:
    :param E_hh:
    :param psi_e:
    :param psi_hh:
    :param well_width:
    :param me:
    :param mhh:
    :param Ep:
    :param nr:
    :param eps:
    :param hwhm:
    :param dimensionality:
    :param line_shape:
    :return:
    """
    if exciton_index < 1:
        raise ValueError("The excition index must start on 1.")

    if dimensionality < 0 or dimensionality > 0.5:
        raise ValueError("The dimensionality parameter of the exciton must be a number between 0 and 0.5 (inclusive).")

    # Borrow constants from the bound-to-bound calculation
    m0 = electron_mass
    C0 = np.pi * q ** 2 / (nr * c * vacuum_permittivity * m0 ** 2 * E / hbar)
    mr = me * mhh / (me + mhh)  # reduced effective mass
    DOS_2D = mr / (np.pi * hbar ** 2 * well_width)
    Mbsq_3D = m0 / 6 * Ep
    Mbsq_2D = 3 / 2 * Mbsq_3D
    Ieh = np.trapz(psi_e * psi_hh, z) ** 2

    # Excitions are considered by modiftying the bulk absorption coefficient using an oscillator strength
    Ry_eff = exciton_rydberg_energy_2d(me, mhh, eps / vacuum_permittivity)
    En = -Ry_eff / ((exciton_index - dimensionality) ** 2)  # Exciton binding energy
    Et = (E_e - E_hh)
    # import pdb; pdb.set_trace()
    # NOTE TO MARKUS: Added the possiblity of using a Gauss lineshape
    if line_shape is "Gauss":
        shape = Gauss(E, Et + En, hwhm)
    else:
        shape = L(E, Et + En, hwhm)
    # NOTE TO MARKUS: It seems that there was a factor 2 missing in the oscilator strength as well as the exciton index. 
    # See: P C Klipstein and N Apsley 1986 J. Phys. C: Solid State Phys. 19 6461 doi:10.1088/0022-3719/19/32/020 
    oscillator_strength = (2 * Ry_eff / (exciton_index - dimensionality) ** 3)
    return C0 * Ieh * Mbsq_2D * DOS_2D * oscillator_strength * shape


def alpha_exciton_elh_TE(exciton_index, E, z, E_e, E_lh, psi_e, psi_lh, well_width, me, mlh, Ep, nr, eps,
                         hwhm=6e-3 * 1.6e-19, dimensionality=0.15, line_shape="Lorenzian"):
    """

    :param exciton_index:
    :param E:
    :param z:
    :param E_e:
    :param E_lh:
    :param psi_e:
    :param psi_lh:
    :param well_width:
    :param me:
    :param mlh:
    :param Ep:
    :param nr:
    :param eps:
    :param hwhm:
    :param dimensionality:
    :param line_shape:
    :return:
    """
    if exciton_index < 1:
        raise ValueError("The excition index must start on 1.")

    if dimensionality < 0 or dimensionality > 0.5:
        raise ValueError("The dimensionality parameter of the exciton must be a number between 0 and 0.5 (inclusive).")

    # Borrow constants from the bound-to-bound calculation
    m0 = electron_mass
    C0 = np.pi * q ** 2 / (nr * c * vacuum_permittivity * m0 ** 2 * E / hbar)
    mr = me * mlh / (me + mlh)  # reduced effective mass
    DOS_2D = mr / (np.pi * hbar ** 2 * well_width)
    Mbsq_3D = m0 / 6 * Ep
    Mbsq_2D = 1 / 2 * Mbsq_3D
    Ieh = np.trapz(psi_e * psi_lh, z) ** 2

    # Excitions are considered by modiftying the bulk absorption coefficient using an oscillator strength
    Ry_eff = exciton_rydberg_energy_2d(me, mlh, eps / vacuum_permittivity)
    En = -Ry_eff / ((exciton_index - dimensionality) ** 2)  # Exciton binding energy
    Et = (E_e - E_lh)
    # NOTE TO MARKUS: Added the possiblity of using a Gauss lineshape
    if line_shape is "Gauss":
        shape = Gauss(E, Et + En, hwhm)
    else:
        shape = L(E, Et + En, hwhm)
    # NOTE TO MARKUS: It seems that there was a factor 2 missing in the oscilator strength as well as the exciton index.
    # See: P C Klipstein and N Apsley 1986 J. Phys. C: Solid State Phys. 19 6461 doi:10.1088/0022-3719/19/32/020
    oscillator_strength = (2 * Ry_eff / (exciton_index - dimensionality) ** 3)
    return C0 * Ieh * Mbsq_2D * DOS_2D * oscillator_strength * shape


def sum_alpha_c_hh_TE(E, z, E_e, E_hh, psi_e, psi_hh, well_width, me, mh, Ep, nr):
    """

    :param E:
    :param z:
    :param E_e:
    :param E_hh:
    :param psi_e:
    :param psi_hh:
    :param well_width:
    :param me:
    :param mh:
    :param Ep:
    :param nr:
    :return:
    """
    alpha = 0
    for ee, pe in zip(E_e, psi_e):
        for eh, ph in zip(E_hh, psi_hh):
            alpha += alpha_c_hh_TE(E, z, ee, eh, pe, ph, well_width, me, mh, Ep, nr)
    return alpha


def sum_alpha_c_lh_TE(E, z, E_e, E_lh, psi_e, psi_lh, well_width, me, mh, Ep, nr):
    """

    :param E:
    :param z:
    :param E_e:
    :param E_lh:
    :param psi_e:
    :param psi_lh:
    :param well_width:
    :param me:
    :param mh:
    :param Ep:
    :param nr:
    :return:
    """
    alpha = 0
    for ee, pe in zip(E_e, psi_e):
        for eh, ph in zip(E_lh, psi_lh):
            alpha += alpha_c_lh_TE(E, z, ee, eh, pe, ph, well_width, me, mh, Ep, nr)
    return alpha


def sum_alpha_exciton_c_hh_TE(E, z, E_e, E_hh, psi_e, psi_hh, well_width, me, mh, Ep, nr, eps, hwhm=6e-3 * 1.6e-19,
                              dimensionality=0.5, line_shape="Lorenzian"):
    """

    :param E:
    :param z:
    :param E_e:
    :param E_hh:
    :param psi_e:
    :param psi_hh:
    :param well_width:
    :param me:
    :param mh:
    :param Ep:
    :param nr:
    :param eps:
    :param hwhm:
    :param dimensionality:
    :param line_shape:
    :return:
    """
    alpha = 0
    eh_pairs_with_zero_angular_momentum = zip(range(len(E_e)),
                                              range(len(E_hh)))  # l=0 for optically allowed transitions
    for exciton_index, eh_quantum_numbers in enumerate(eh_pairs_with_zero_angular_momentum):
        e_index = eh_quantum_numbers[0]
        hh_index = eh_quantum_numbers[1]
        ee = E_e[e_index]
        ehh = E_hh[hh_index]
        pe = psi_e[e_index]
        phh = psi_hh[e_index]
        alpha += alpha_exciton_ehh_TE(exciton_index + 1, E, z, ee, ehh, pe, phh, well_width, me, mh, Ep, nr, eps,
                                      hwhm=hwhm, dimensionality=dimensionality, line_shape=line_shape)
    return alpha


def sum_alpha_exciton_c_lh_TE(E, z, E_e, E_lh, psi_e, psi_lh, well_width, me, mlh, Ep, nr, eps, hwhm=6e-3 * 1.6e-19,
                              dimensionality=0.5, line_shape="Lorenzian"):
    """

    :param E:
    :param z:
    :param E_e:
    :param E_lh:
    :param psi_e:
    :param psi_lh:
    :param well_width:
    :param me:
    :param mlh:
    :param Ep:
    :param nr:
    :param eps:
    :param hwhm:
    :param dimensionality:
    :param line_shape:
    :return:
    """
    alpha = 0
    eh_pairs_with_zero_angular_momentum = zip(range(len(E_e)),
                                              range(len(E_lh)))  # l=0 for optically allowed transitions
    for exciton_index, eh_quantum_numbers in enumerate(eh_pairs_with_zero_angular_momentum):
        e_index = eh_quantum_numbers[0]
        lh_index = eh_quantum_numbers[1]
        ee = E_e[e_index]
        elh = E_lh[lh_index]
        pe = psi_e[e_index]
        plh = psi_lh[e_index]
        alpha += alpha_exciton_elh_TE(exciton_index + 1, E, z, ee, elh, pe, plh, well_width, me, mlh, Ep, nr, eps,
                                      hwhm=hwhm, dimensionality=dimensionality, line_shape=line_shape)
    return alpha


def calc_alpha(QM_result, well_width, kane_parameter=28 * 1.6e-19, refractive_index=3.5, hwhm=6e-3 * 1.6e-19,
               dimensionality=0.5, theta=0, eps=12.9 * vacuum_permittivity, espace=None, line_shape="Lorenzian"):
    """ Calculates the absorption coeficient of a quantum well structure assuming the parabolic approximation for the
    effective masses.

    :param QM_result: The output of the Schrodinger solver, incldued in the 'quantum_mechanics' package
    :param well_width: The well width
    :param kane_parameter: The Kane parameter
    :param refractive_index: Refractive (effective) index of the QW
    :param hwhm: Full width at half maximum of the excitonic lineshape
    :param dimensionality:
    :param theta:
    :param eps:
    :param espace:
    :param line_shape:
    :return:
    """
    results = QM_result

    # We pick the effective mass at the center of the structure as THE efective mass of the QW
    points = len(results["effective_masses"]["me"])
    me = results["effective_masses"]["me"][round(points / 2.0)]
    try:
        mhh = results["effective_masses"]["mhh"][round(points / 2.0)]
        mlh = results["effective_masses"]["mlh"][round(points / 2.0)]
    except KeyError:
        # If the calculation mode for the bands is kp4x4 or kp6x6, we have to use the effective mass transverse to the
        # growth direction (in the QW plane).
        mhh = results["effective_masses"]["mhh_t"][round(points / 2.0)]
        mlh = results["effective_masses"]["mlh_t"][round(points / 2.0)]
    #    print (results["effective_masses"]["me"][round(points/2.0)])
    #    print (results["effective_masses"]["mhh"][round(points/2.0)])
    #    print (results["effective_masses"]["mlh"][round(points/2.0)])

    result_c_hh_TE = []
    result_c_lh_TE = []
    result_ex_c_hh_TE = []
    result_ex_c_lh_TE = []
    for e in espace:
        c_hh = sum_alpha_c_hh_TE(e, results['x'], results['E']['Ee'], results['E']['Ehh'],
                                 results["wavefunctions"]['psi_e'], results["wavefunctions"]['psi_hh'], well_width, me,
                                 mhh, kane_parameter, refractive_index)
        c_lh = sum_alpha_c_lh_TE(e, results['x'], results['E']['Ee'], results['E']['Elh'],
                                 results["wavefunctions"]['psi_e'], results["wavefunctions"]['psi_lh'], well_width, me,
                                 mlh, kane_parameter, refractive_index)
        ex_c_hh = sum_alpha_exciton_c_hh_TE(e, results['x'], results['E']['Ee'], results['E']['Ehh'],
                                            results["wavefunctions"]['psi_e'], results["wavefunctions"]['psi_hh'],
                                            well_width, me, mhh, kane_parameter, refractive_index, eps, hwhm=hwhm,
                                            dimensionality=dimensionality, line_shape=line_shape)
        ex_c_lh = sum_alpha_exciton_c_lh_TE(e, results['x'], results['E']['Ee'], results['E']['Elh'],
                                            results["wavefunctions"]['psi_e'], results["wavefunctions"]['psi_lh'],
                                            well_width, me, mlh, kane_parameter, refractive_index, eps, hwhm=hwhm,
                                            dimensionality=dimensionality, line_shape=line_shape)
        result_c_hh_TE.append(c_hh)
        result_c_lh_TE.append(c_lh)
        result_ex_c_hh_TE.append(ex_c_hh)
        result_ex_c_lh_TE.append(ex_c_lh)

    result_c_hh_TE = array(result_c_hh_TE)
    result_c_lh_TE = array(result_c_lh_TE)
    result_c_lh_TM = 4 * array(result_c_lh_TE)
    result_ex_c_hh_TE = array(result_ex_c_hh_TE)
    result_ex_c_lh_TE = array(result_ex_c_lh_TE)
    result_ex_c_lh_TM = 4 * array(result_ex_c_lh_TE)
    result_TE_raw = result_c_hh_TE + result_c_lh_TE + result_ex_c_hh_TE + result_ex_c_lh_TE
    result_TE = result_TE_raw * (np.cos(theta) ** 2 + 0.5 * np.sin(theta) ** 2)
    result_TM_raw = result_c_lh_TM + result_ex_c_lh_TM
    result_TM = 0.5 * result_TM_raw * np.sin(theta) ** 2
    result_sum = result_TE + result_TM
    return espace, result_sum, result_TE, result_TM, result_TE_raw, result_TM_raw


def NonBlackBodyEmission(E, voltage=0, nr=3.5, T=300):
    """

    :param E:
    :param voltage:
    :param nr:
    :param T:
    :return:
    """
    A = 2 * nr ** 2 / h ** 3 / c ** 2
    B = (E - q * voltage) / kb / T

    bb = A * E ** 2 / (np.exp(B) - 1)

    return bb


def calc_emission(QM_result, well_width, voltage=0, theta=0):
    """

    :param QM_result:
    :param well_width:
    :param voltage:
    :param theta:
    :return:
    """
    alfa = QM_result["alpha"]
    aTE = alfa[4]
    aTM = alfa[5]

    result_TE = []
    result_TM = []
    result = []
    for index, e in enumerate(alfa[0]):
        bb = NonBlackBodyEmission(e)
        TE = bb * (1 - np.exp(-aTE[index] * well_width * (np.cos(theta) ** 2 + 0.5 * np.sin(theta) ** 2)))
        TM = bb * (1 - np.exp(-aTM[index] * well_width * 0.5 * np.sin(theta) ** 2))
        # All = bb*(1-np.exp(-(aTE[index]+aTM[index])*well_width))
        All = TE + TM
        result_TE.append(TE)
        result_TM.append(TM)
        result.append(All)

    result_TE = array(result_TE)
    result_TM = array(result_TM)
    result = array(result)

    return alfa[0], result, result_TE, result_TM
