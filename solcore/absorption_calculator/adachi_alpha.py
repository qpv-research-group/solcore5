import numpy as np
from solcore import get_parameter
from solcore.constants import hbar, pi, c
from solcore.science_tracker import science_reference


def create_adachi_alpha(material, Esteps=(1.42, 6, 3000), T=300, wl=None):
    """ Calculates the n, k and absorption coefficient of a material using Adachi's formalism of critical points.

    :param material: A solcore material
    :param Esteps: (1.42, 6, 3000) A tuple with the start, end and step energies in which calculating the optical data
    :param T: (300) Temeprature in kelvin
    :param wl: (None) Optional array indicating the wavelengths in which calculating the data
    :return: A tuple containing 4 arrays: (Energy, n, k, alpha)
    """

    science_reference("Adachi optical dispersion relations",
                      'S.Adachi, “Optical dispersion relations for GaP, GaAs, GaSb, InP, InAs, InSb, AlxGa1−xAs, '
                      'and In1−xGaxAsyP1−y” J.Appl.Phys., vol.66, no.12, pp.6030–12, 1989.')

    if type(material) != str:
        material = material.plain_string()

    e0 = get_parameter(material, "E0", T=T) + 0j
    Delta0 = get_parameter(material, "E0plusD0", T=T) + 0j - e0
    e1 = get_parameter(material, "E1", T=T) + 0j
    Delta1 = get_parameter(material, "E1plusD1", T=T) + 0j - e1
    e2 = get_parameter(material, "E2", T=T) + 0j
    egid = get_parameter(material, "Eg1d", T=T) + 0j
    a = get_parameter(material, "A", T=T) + 0j
    b1 = get_parameter(material, "B1", T=T) + 0j
    b11 = get_parameter(material, "B11", T=T) + 0j
    Gamma = get_parameter(material, "Gamma", T=T) + 0j  # , "eV",T=T)
    cc = get_parameter(material, "C", T=T) + 0j
    gamma = get_parameter(material, "gamma", T=T) + 0j
    d = get_parameter(material, "D", T=T) + 0j
    omegaphonon = 0
    b2 = b1
    b21 = b11

    a0 = get_parameter(material, "lattice_constant", T=T)

    b2 = 44 * (e1 + 2 * Delta1 / 3.) / (a0 * (e1 + Delta1) ** 2)
    b1 = 44 * (e1 + Delta1 / 3) / (a0 * e1 ** 2)
    b2 = 0
    b1 = 7
    b11 = 2 * b1
    b21 = 2 * b2

    other_params = {
        "$e_0$": e0.real + 0.005,
        "$e_1$": e1.real,
        "$e_1+\Delta_1$": (e1 + Delta1).real,
        "$e_2$": e2.real,
        "$e_0+\Delta_0$": (e0 + Delta0).real,
        "$e_{g,1d}$": egid.real
    }
    # print b1, b2, b11,b21

    # print other_params
    if wl is not None:
        E = 1240 * 1e-9 / wl[::-1]
    else:
        E = np.linspace(*Esteps)  # electron Volts

    omega = E * 1.6e-19 / hbar

    re = lambda x: x.real
    H = lambda x: (x.real >= 0)
    i = 1j

    XiSO = E / (e0 + Delta0)
    Xi0 = E / e0

    f_Xi_0 = (Xi0) ** -2 * (2 - (1 + Xi0) ** 0.5 - (1 - Xi0) ** 0.5 * H(1 - Xi0))
    f_Xi_s0 = (XiSO) ** -2 * (2 - (1 + XiSO) ** 0.5 - (1 - XiSO) ** 0.5 * H(1 - XiSO))  #

    # equation 7
    epsilon_2_e0 = (a / E ** 2) * (
        (E - e0) ** 0.5 * H((Xi0) - 1) + 0.5 * (E - e0 - Delta0) ** 0.5 * H((E / (e0 + Delta0)) - 1))  #
    epsilon_1_e0 = a * e0 ** -1.5 * (f_Xi_0 + 0.5 * (e0 / (e0 + Delta0)) ** 1.5 * f_Xi_s0)

    Xi1 = E / e1  #
    Xi2 = E / e2  #
    Xi1s = E / (e1 + Delta1)  #

    # equation 12
    epsilon_2_e1 = pi * Xi1 ** -2 * (b1 - b11 * (e1 - E) ** 0.5 * H(1 - Xi1)) * H(
        re(b1 - b11 * (e1 - E) ** 0.5 * H(1 - Xi1)))  #
    epsilon_2_e1_d1 = pi * Xi1s ** -2 * (b2 - b21 * (e1 + Delta1 - E) ** 0.5 * H(
        1 - Xi1s))  # where is this from? #*H(b2-b21*(e1+Delta1-E)**0.5*H(1-Xi1s)) #
    # from numpy import sin
    #    epsilon_2_3dcp = pi*b1*Xi1**-2*H(Xi1-1) + pi*b2*Xi1s**-2*H(Xi1s-1) # +10*sin(30*E)#
    # epsilon_1_3dcp = -b1*(Xi1+i*Gamma/e1)**-2*log(1-(Xi1+i*Gamma/e1)**2) # Ned Version
    # epsilon_1_3dcp = -b1*Xi1**-2*log(1-Xi1**2)-b2*Xi1s**-2*log(1-Xi1s**2) #adachi version, no damping
    Gamma = 0.06
    E_damped = E + i * Gamma
    Xi1damped = E_damped / e1
    Xi1sdamped = E_damped / (e1 + Delta1)

    # equation 16
    epsilon_1_3dcp = -b1 * Xi1damped ** -2 * np.log(1 - Xi1damped ** 2) - b2 * Xi1sdamped ** -2 * np.log(
        1 - Xi1sdamped ** 2)  # adachi version
    epsilon_2_3dcp = pi * b1 * Xi1damped ** -2 * H(Xi1 - 1) + pi * b2 * Xi1sdamped ** -2 * H(Xi1s - 1)  # +10*sin(30*E)#

    # equation 16
    epsilon_2_e2 = cc * Xi2 * gamma / ((1 - Xi2 ** 2) ** 2 + (Xi2 * gamma) ** 2) * H(
        E - e1)  # added last term for compatibility with adachi's graph, dammit.
    # # print log(epsilon_2_e2)
    epsilon_1_e2 = cc * (1 - Xi2 ** 2) / ((1 - Xi2 ** 2) ** 2 + (Xi2 * gamma) ** 2)  #
    # print egid
    XiG_plus = (egid + hbar * omegaphonon) / E
    XiG_minus = (egid - hbar * omegaphonon) / E
    XiCh = E / e1

    epsilon_2_indirect = d / E ** 2 * (E - egid + hbar * omegaphonon) ** 2 * H(1 - XiG_plus) * H(1 - XiCh)  #
    epsilon_2_indirect_b = d / E ** 2 * (E - egid - hbar * omegaphonon) ** 2 * H(1 - XiG_minus) * H(1 - XiCh)  #
    # epsilon_2_indirect_minus = d/E**2 *(E-egid - hbar*omegaphonon)**2 * H(1-XiG_plus) * H(1-XiCh)# markus made this up a bit
    # epsilon_2_indirect_b_minus = d/E**2 *(E-egid + hbar*omegaphonon)**2 * H(1-XiG_minus) * H(1-XiCh)#markus made this up a bit

    epsilon_2_indirect_total = epsilon_2_indirect + epsilon_2_indirect_b

    eps1 = epsilon_1_e0 + epsilon_1_3dcp + epsilon_1_e2
    eps2 = epsilon_2_e0 + epsilon_2_e1 + epsilon_2_e1_d1 + epsilon_2_e2 + epsilon_2_indirect_total

    eps2 = epsilon_2_e0 + epsilon_2_e1 + epsilon_2_e2 + epsilon_2_indirect_total

    E = E * 1.6e-19
    k = np.sqrt((np.sqrt(eps1 ** 2 + eps2 ** 2) - eps1) / 2).real
    n = np.sqrt((np.sqrt(eps1 ** 2 + eps2 ** 2) + eps1) / 2).real
    alpha_data = 2 * omega / c * k

    # If the input was in wavelengths, we reverse the output arrays
    if wl is not None:
        alpha_data = alpha_data[::-1]
        n = n[::-1]
        k = k[::-1]

    return E.real, abs(n + 1e-10), abs(k + 1e-10), abs(alpha_data + 1e-10)


if __name__ == '__main__':
    from solcore import material
    from solcore.graphing import Graph, GraphData

    material_name = "InGaAs"
    solcore_material = material(material_name)(T=300, In=0.1)

    E, nn, kk, adachi_alpha_data = create_adachi_alpha(solcore_material, T=300)

    idx = np.argmin(abs(E-2))

    solcore_alpha_data = solcore_material.alphaE(E)

    E /= 1.6e-19

    idx = np.argmin(abs(E - 2))
    print(round(E[idx]*1000), round(nn[idx]*1000), round(kk[idx]*1000))

    out = [round(E[idx]*1000), round(nn[idx]*1000), round(kk[idx]*1000)]

    data = [2000.0, 3629.0, 308.0]

    print(data==out)

    curves = [
        GraphData(E, adachi_alpha_data / 100, label="adachi-Generated Data", linewidth=2, color="red"),
        GraphData(E, solcore_alpha_data / 100, label="Default Solcore Data", linewidth=2, color="grey"),
    ]

    g = Graph(curves, yscale="log", legend="best", xlim=(0, 6), ylim=(1000, 1e7)).draw()

    curves = [
        GraphData(E, nn, label="n", linewidth=2, color="red"),
        GraphData(E, kk, label="k", linewidth=2, color="grey"),
    ]

    g = Graph(curves, yscale="log", legend="best", xlim=(0, 6), ylim=(1e-3, 10)).draw()
