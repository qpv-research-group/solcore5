# -*- coding: utf-8 -*-
"""
This is mostly based on Steven Byrnes' tmm package, with modifications to vectorize the code
over wavelengths (by Diego Alonso-Alvarez), and to include depth-dependent absorption calculations
in incoherent layers using the Beer-Lambert law (by Phoebe Pearce).

All credit of the algorithm, testing, etc. goes to Steven Byrnes. For more details, visit:

- https://arxiv.org/abs/1603.02720

- http://sjbyrnes.com

- For information see the docstring of each function, and also see manual.pdf

The most two important functions are:

coh_tmm(...) -- the transfer-matrix-method calculation in the coherent
case (i.e. thin films)

inc_tmm(...) -- the transfer-matrix-method calculation in the incoherent
case (i.e. films tens or hundreds of wavelengths thick, or whose
thickness is not very uniform.

These functions are all imported into the main package (tmm) namespace,
so you can call them with tmm.coh_tmm(...) etc.
"""


import scipy as sp
import numpy as np

import sys

EPSILON = sys.float_info.epsilon  # typical floating-point calculation error


def make_2x2_array(a, b, c, d, dtype=float):
    """
    Makes a 2x2 numpy array of [[a,b],[c,d]]

    Same as "numpy.array([[a,b],[c,d]], dtype=float)", but ten times faster
    """
    my_array = np.empty((len(a), 2, 2), dtype=dtype)
    my_array[:, 0, 0] = a
    my_array[:, 0, 1] = b
    my_array[:, 1, 0] = c
    my_array[:, 1, 1] = d
    return my_array


def snell(n_1, n_2, th_1):
    """
    return angle theta in layer 2 with refractive index n_2, assuming
    it has angle th_1 in layer with refractive index n_1. Use Snell's law. Note
    that "angles" may be complex!!
    """
    # Important that the arcsin here is scipy.arcsin, not numpy.arcsin!! (They
    # give different results e.g. for arcsin(2).)
    # Use real_if_close because e.g. arcsin(2 + 1e-17j) is very different from
    # arcsin(2) due to branch cut
    return sp.arcsin(np.real_if_close(n_1 * np.sin(th_1) / n_2))


def list_snell(n_list, th_0):
    """
    return list of angle theta in each layer based on angle th_0 in layer 0,
    using Snell's law. n_list is index of refraction of each layer. Note that
    "angles" may be complex!!
    """
    # Important that the arcsin here is scipy.arcsin, not numpy.arcsin!! (They
    # give different results e.g. for arcsin(2).)
    # Use real_if_close because e.g. arcsin(2 + 1e-17j) is very different from
    # arcsin(2) due to branch cut
    return sp.arcsin(np.real_if_close(n_list[0] * np.sin(th_0) / n_list))


def interface_r(polarization, n_i, n_f, th_i, th_f):
    """
    reflection amplitude (from Fresnel equations)

    polarization is either "s" or "p" for polarization

    n_i, n_f are (complex) refractive index for incident and final

    th_i, th_f are (complex) propegation angle for incident and final
    (in radians, where 0=normal). "th" stands for "theta".
    """
    if polarization == 's':
        # return 2 * n_i * np.cos(th_i) / (n_i * np.cos(th_i) + n_f * np.cos(th_f))
        return ((n_i * np.cos(th_i) - n_f * np.cos(th_f)) /
                (n_i * np.cos(th_i) + n_f * np.cos(th_f)))
    elif polarization == 'p':
        return ((n_f * np.cos(th_i) - n_i * np.cos(th_f)) /
                (n_f * np.cos(th_i) + n_i * np.cos(th_f)))
    else:
        raise ValueError("Polarization must be 's' or 'p'")


def interface_t(polarization, n_i, n_f, th_i, th_f):
    """
    transmission amplitude (frem Fresnel equations)

    polarization is either "s" or "p" for polarization

    n_i, n_f are (complex) refractive index for incident and final

    th_i, th_f are (complex) propegation angle for incident and final
    (in radians, where 0=normal). "th" stands for "theta".
    """
    if polarization == 's':
        return 2 * n_i * np.cos(th_i) / (n_i * np.cos(th_i) + n_f * np.cos(th_f))
    elif polarization == 'p':
        return 2 * n_i * np.cos(th_i) / (n_f * np.cos(th_i) + n_i * np.cos(th_f))
    else:
        raise ValueError("Polarization must be 's' or 'p'")


def R_from_r(r):
    """
    Calculate reflected power R, starting with reflection amplitude r.
    """
    return abs(r) ** 2


def T_from_t(pol, t, n_i, n_f, th_i, th_f):
    """
    Calculate transmitted power T, starting with transmission amplitude t.

    n_i,n_f are refractive indices of incident and final medium.

    th_i, th_f are (complex) propegation angles through incident & final medium
    (in radians, where 0=normal). "th" stands for "theta".

    In the case that n_i,n_f,th_i,th_f are real, formulas simplify to
    T=|t|^2 * (n_f cos(th_f)) / (n_i cos(th_i)).

    See manual for discussion of formulas
    """
    if (pol == 's'):
        return abs(t ** 2) * (((n_f * np.cos(th_f)).real) / (n_i * np.cos(th_i)).real)
    elif (pol == 'p'):
        return abs(t ** 2) * (((n_f * np.conj(np.cos(th_f))).real) /
                              (n_i * np.conj(np.cos(th_i))).real)
    else:
        raise ValueError("Polarization must be 's' or 'p'")


def power_entering_from_r(pol, r, n_i, th_i):
    """
    Calculate the power entering the first interface of the stack, starting with
    reflection amplitude r. Normally this equals 1-R, but in the unusual case
    that n_i is not real, it can be a bit different than 1-R. See manual.

    n_i is refractive index of incident medium.

    th_i is (complex) propegation angle through incident medium
    (in radians, where 0=normal). "th" stands for "theta".
    """
    if (pol == 's'):
        return ((n_i * np.cos(th_i) * (1 + np.conj(r)) * (1 - r)).real
                / (n_i * np.cos(th_i)).real)
    elif (pol == 'p'):
        return ((n_i * np.conj(np.cos(th_i)) * (1 + r) * (1 - np.conj(r))).real
                / (n_i * np.conj(np.cos(th_i))).real)
    else:
        raise ValueError("Polarization must be 's' or 'p'")


def interface_R(polarization, n_i, n_f, th_i, th_f):
    """
    Fraction of light intensity reflected at an interface.
    """
    r = interface_r(polarization, n_i, n_f, th_i, th_f)
    return R_from_r(r)


def interface_T(polarization, n_i, n_f, th_i, th_f):
    """
    Fraction of light intensity transmitted at an interface.
    """
    t = interface_t(polarization, n_i, n_f, th_i, th_f)
    return T_from_t(polarization, t, n_i, n_f, th_i, th_f)


def coh_tmm(pol, n_list, d_list, th_0, lam_vac):
    """
    This function is vectorized.
    Main "coherent transfer matrix method" calc. Given parameters of a stack,
    calculates everything you could ever want to know about how light
    propagates in it. (If performance is an issue, you can delete some of the
    calculations without affecting the rest.)

    pol is light polarization, "s" or "p".

    n_list is the list of refractive indices, in the order that the light would
    pass through them. The 0'th element of the list should be the semi-infinite
    medium from which the light enters, the last element should be the semi-
    infinite medium to which the light exits (if any exits).

    th_0 is the angle of incidence: 0 for normal, pi/2 for glancing.
    Remember, for a dissipative incoming medium (n_list[0] is not real), th_0
    should be complex so that n0 sin(th0) is real (intensity is constant as
    a function of lateral position).

    d_list is the list of layer thicknesses (front to back). Should correspond
    one-to-one with elements of n_list. First and last elements should be "inf".

    lam_vac is vacuum wavelength of the light.

    Outputs the following as a dictionary (see manual for details)

    * r--reflection amplitude
    * t--transmission amplitude
    * R--reflected wave power (as fraction of incident)
    * T--transmitted wave power (as fraction of incident)
    * power_entering--Power entering the first layer, usually (but not always)
      equal to 1-R (see manual).
    * vw_list-- n'th element is [v_n,w_n], the forward- and backward-traveling
      amplitudes, respectively, in the n'th medium just after interface with
      (n-1)st medium.
    * kz_list--normal component of complex angular wavenumber for
      forward-traveling wave in each layer.
    * th_list--(complex) propagation angle (in radians) in each layer
    * pol, n_list, d_list, th_0, lam_vac--same as input

    """
    # convert lists to numpy arrays if they're not already.
    n_list = np.array(n_list)
    d_list = np.array(d_list, dtype=float)[:, None]

    # input tests
    if hasattr(th_0, 'size') and th_0.size > 1 and th_0.size != lam_vac.size:
        raise ValueError('This function is not vectorized for angles; you need to run one angle calculation at a time.')
    if n_list.shape[0] != d_list.shape[0]:
        raise ValueError("Problem with n_list or d_list!")
    if (d_list[0] != np.inf) or (d_list[-1] != np.inf):
        raise ValueError('d_list must start and end with inf!')
    if any(abs((n_list[0] * np.sin(th_0)).imag) > 100 * EPSILON):
        raise ValueError('Error in n0 or th0!')
    if hasattr(th_0, 'size'):
        th_0 = np.array(th_0)
    num_layers = n_list.shape[0]
    num_wl = n_list.shape[1]

    # th_list is a list with, for each layer, the angle that the light travels
    # through the layer. Computed with Snell's law. Note that the "angles" may be
    # complex!
    th_list = list_snell(n_list, th_0)

    # kz is the z-component of (complex) angular wavevector for forward-moving
    # wave. Positive imaginary part means decaying.
    kz_list = 2 * np.pi * n_list * np.cos(th_list) / lam_vac

    # delta is the total phase accrued by traveling through a given layer.
    # ignore warning about inf multiplication
    olderr = sp.seterr(invalid='ignore')
    delta = kz_list * d_list
    sp.seterr(**olderr)

    # For a very opaque layer, reset delta to avoid divide-by-0 and similar
    # errors. The criterion imag(delta) > 35 corresponds to single-pass
    # transmission < 1e-30 --- small enough that the exact value doesn't
    # matter.
    # It DOES matter (for depth-dependent calculations!)
    delta[1:num_layers - 1, :] = np.where(delta[1:num_layers - 1, :].imag > 100, delta[1:num_layers - 1, :].real + 100j,
                                          delta[1:num_layers - 1, :])

    # t_list[i,j] and r_list[i,j] are transmission and reflection amplitudes,
    # respectively, coming from i, going to j. Only need to calculate this when
    # j=i+1. (2D array is overkill but helps avoid confusion.)
    t_list = np.zeros((num_wl, num_layers, num_layers), dtype=complex)
    r_list = np.zeros((num_wl, num_layers, num_layers), dtype=complex)

    for i in range(num_layers - 1):
        t_list[:, i, i + 1] = interface_t(pol, n_list[i], n_list[i + 1],
                                          th_list[i], th_list[i + 1])
        r_list[:, i, i + 1] = interface_r(pol, n_list[i], n_list[i + 1],
                                          th_list[i], th_list[i + 1])
    # At the interface between the (n-1)st and nth material, let v_n be the
    # amplitude of the wave on the nth side heading forwards (away from the
    # boundary), and let w_n be the amplitude on the nth side heading backwards
    # (towards the boundary). Then (v_n,w_n) = M_n (v_{n+1},w_{n+1}). M_n is
    # M_list[n]. M_0 and M_{num_layers-1} are not defined.
    # My M is a bit different than Sernelius's, but Mtilde is the same.

    M_list = np.zeros((num_layers, num_wl, 2, 2), dtype=complex)
    for i in range(1, num_layers - 1):
        A = make_2x2_array(np.exp(-1j * delta[i]), np.zeros_like(delta[i]), np.zeros_like(delta[i]), np.exp(1j * delta[i]),
                           dtype=complex)
        B = make_2x2_array(np.ones_like(delta[i]), r_list[:, i, i + 1], r_list[:, i, i + 1], np.ones_like(delta[i]),
                           dtype=complex)
        d = (1 / t_list[:, i, i + 1])

        M_list[i] = np.transpose(d * np.transpose(np.matmul(A, B)))  # , (1, 2, 0)), (2, 0, 1))

    Mtilde = make_2x2_array(np.ones_like(delta[i]), np.zeros_like(delta[i]), np.zeros_like(delta[i]),
                            np.ones_like(delta[i]), dtype=complex)
    for i in range(1, num_layers - 1):
        Mtilde = np.matmul(Mtilde, M_list[i])

    A = make_2x2_array(np.ones_like(delta[i]), r_list[:, 0, 1], r_list[:, 0, 1], np.ones_like(delta[i]),
                       dtype=complex)
    d = 1 / t_list[:, 0, 1]
    Mtilde = np.matmul(np.transpose(d * np.transpose(A, (1, 2, 0)), (2, 0, 1)), Mtilde)

    # Net complex transmission and reflection amplitudes
    r = Mtilde[:, 1, 0] / Mtilde[:, 0, 0]
    t = np.ones_like(Mtilde[:, 0, 0]) / Mtilde[:, 0, 0]

    # vw_list[n] = [v_n, w_n]. v_0 and w_0 are undefined because the 0th medium
    # has no left interface.
    vw_list = np.zeros((num_layers, num_wl, 2), dtype=complex)
    vw = np.zeros((num_wl, 2, 2), dtype=complex)
    I = np.identity(2)
    vw[:, 0, 0] = t
    vw[:, 0, 1] = t
    vw_list[-1] = vw[:, 0, :]
    for i in range(num_layers - 2, 0, -1):
        vw = np.matmul(M_list[i], vw)
        vw_list[i, :, :] = vw[:, :, 1]
    vw_list[-1, :, 1] = 0

    # Net transmitted and reflected power, as a proportion of the incoming light
    # power.
    R = R_from_r(r)
    T = T_from_t(pol, t, n_list[0], n_list[-1], th_0, th_list[-1])
    power_entering = power_entering_from_r(
        pol, r, n_list[0], th_0)

    return {'r': r, 't': t, 'R': R, 'T': T, 'power_entering': power_entering,
            'vw_list': vw_list, 'kz_list': kz_list, 'th_list': th_list,
            'pol': pol, 'n_list': n_list, 'd_list': d_list, 'th_0': th_0,
            'lam_vac': lam_vac}

    # return {'r': r, 't': t, 'R': R, 'T': T, 'power_entering': power_entering,
    #         'kz_list': kz_list, 'th_list': th_list,
    #         'pol': pol, 'n_list': n_list, 'd_list': d_list, 'th_0': th_0,
    #         'lam_vac': lam_vac}


def coh_tmm_reverse(pol, n_list, d_list, th_0, lam_vac):
    """
    Reverses the order of the stack then runs coh_tmm.
    """
    th_f = snell(n_list[0], n_list[-1], th_0)
    return coh_tmm(pol, n_list[::-1], d_list[::-1], th_f, lam_vac)


def ellips(n_list, d_list, th_0, lam_vac):
    """
    Calculates ellipsometric parameters, in radians.

    Warning: Conventions differ. You may need to subtract pi/2 or whatever.
    """

    s_data = coh_tmm('s', n_list, d_list, th_0, lam_vac)
    p_data = coh_tmm('p', n_list, d_list, th_0, lam_vac)
    rs = s_data['r']
    rp = p_data['r']
    return {'psi': np.arctan(abs(rp / rs)), 'Delta': np.angle(-rp / rs)}


def unpolarized_RT(n_list, d_list, th_0, lam_vac):
    """
    This function is vectorized.
    Calculates reflected and transmitted power for unpolarized light.
    """

    s_data = coh_tmm('s', n_list, d_list, th_0, lam_vac)
    p_data = coh_tmm('p', n_list, d_list, th_0, lam_vac)
    R = (s_data['R'] + p_data['R']) / 2.
    T = (s_data['T'] + p_data['T']) / 2.
    return {'R': R, 'T': T}


def position_resolved(layer, dist, coh_tmm_data):
    """
    This function is vectorized.
    Starting with output of coh_tmm(), calculate the Poynting vector
    and absorbed energy density a distance "dist" into layer number "layer"
    """

    vw = coh_tmm_data['vw_list'][layer]
    kz = coh_tmm_data['kz_list'][layer]
    th = coh_tmm_data['th_list'][layer]
    n = coh_tmm_data['n_list'][layer]
    n_0 = coh_tmm_data['n_list'][0]
    th_0 = coh_tmm_data['th_0']
    pol = coh_tmm_data['pol']

    # amplitude of forward-moving wave is Ef, backwards is Eb

    Ef = (vw.T[0] * np.exp(1j * kz.T * dist)).T
    Eb = (vw.T[1] * np.exp(-1j * kz.T * dist)).T


    # Poynting vector
    if (pol == 's'):
        poyn = ((n * np.cos(th) * np.conj(Ef + Eb) * (Ef - Eb)).real) / (n_0 * np.cos(th_0)).real
    elif (pol == 'p'):
        poyn = (((n * np.conj(np.cos(th)) * (Ef + Eb) * np.conj(Ef - Eb)).real)
                / (n_0 * np.conj(np.cos(th_0))).real)

    # absorbed energy density
    if (pol == 's'):
        absor = (n * np.cos(th) * kz * abs(Ef + Eb) ** 2).imag / (n_0 * np.cos(th_0)).real
    elif (pol == 'p'):
        absor = (n * np.conj(np.cos(th)) *
                 (kz * abs(Ef - Eb) ** 2 - np.conj(kz) * abs(Ef + Eb) ** 2)
                 ).imag / (n_0 * np.conj(np.cos(th_0))).real


    return ({'poyn': poyn.T, 'absor': absor.T})


def find_in_structure(d_list, dist):
    """
    This function is vectorized.
    d_list is list of thicknesses of layers, all of which are finite.

    dist is the distance from the front of the whole multilayer structure
    (i.e., from the start of layer 0.)

    Function returns [layer,z], where:

    layer is what number layer you're at.
    (For large enough dist, layer = len(d_list), even though d_list[layer]
    doesn't exist in that case.

    z is the distance into that layer.
    """
    if sum(d_list) == np.inf:
        raise ValueError('This function expects finite arguments')
    # layer = 0

    d_array = np.array([0] + d_list)
    cum_sum = d_array.cumsum()
    layer = np.searchsorted(cum_sum, dist, side='right')
    dist = dist - cum_sum[layer - 1]

    return layer, dist


def find_in_structure_with_inf(d_list, dist):
    """
    d_list is list of thicknesses of layers [inf, blah, blah, ..., blah, inf]

    dist is the distance from the front of the whole multilayer structure
    (i.e., frcom the start of layer 1.)

    Function returns [layer,z], where:

    layer is what number layer you're at,

    z is the distance into that layer.
    """
    layer, dist = find_in_structure(d_list[1:-1], dist)
    return layer, dist


def layer_starts(d_list):
    """
    Gives the location of the start of any given layer, relative to the front
    of the whole multilayer structure. (i.e. the start of layer 1)

    d_list is list of thicknesses of layers [inf, blah, blah, ..., blah, inf]

    """
    final_answer = np.zeros(len(d_list))
    final_answer[0] = -np.inf
    final_answer[1] = 0
    for i in range(2, len(d_list)):
        final_answer[i] = final_answer[i - 1] + d_list[i - 1]
    return final_answer


class absorp_analytic_fn:
    """
    This function (specifically, 'run') is vectorized.
    Absorption in a given layer is a pretty simple analytical function:
    The sum of four exponentials.

    a(z) = A1*exp(a1*z) + A2*exp(-a1*z)
           + A3*exp(1j*a3*z) + conj(A3)*exp(-1j*a3*z)

    where a(z) is absorption at depth z, with z=0 being the start of the layer,
    and A1,A2,a1,a3 are real numbers, with a1>0, a3>0, and A3 is complex.
    The class stores these five parameters, as well as d, the layer thickness.

    This gives absorption as a fraction of intensity coming towards the first
    layer of the stack.
    """

    def fill_in(self, coh_tmm_data, layer):
        """
        fill in the absorption analytic function starting from coh_tmm_data
        (the output of coh_tmm), for absorption in the layer with index
        "layer".
        """

        pol = coh_tmm_data['pol']
        v = coh_tmm_data['vw_list'][layer, :, 0]
        w = coh_tmm_data['vw_list'][layer, :, 1]
        kz = coh_tmm_data['kz_list'][layer]
        n = coh_tmm_data['n_list'][layer]
        n_0 = coh_tmm_data['n_list'][0]
        th_0 = coh_tmm_data['th_0']
        th = coh_tmm_data['th_list'][layer]
        self.d = coh_tmm_data['d_list'][layer]

        self.a1 = 2 * kz.imag
        self.a3 = 2 * kz.real

        self.a1[self.a1 < 1e-30] = 0
        self.a3[self.a3 < 1e-30] = 0

        if pol == 's':
            temp = (n * np.cos(th) * kz).imag / (n_0 * np.cos(th_0)).real
            self.A1 = temp * abs(w) ** 2
            self.A2 = temp * abs(v) ** 2
            self.A3 = temp * v * np.conj(w)
        else:  # pol=='p'
            temp = (2 * (kz.imag) * (n * np.cos(np.conj(th))).real /
                    (n_0 * np.conj(np.cos(th_0))).real)
            self.A1 = temp * abs(w) ** 2
            self.A2 = temp * abs(v) ** 2
            self.A3 = v * np.conj(w) * (-2 * (kz.real) * (n * np.cos(np.conj(th))).imag /
                                     (n_0 * np.conj(np.cos(th_0))).real)
        return self

    def copy(self):
        """
        Create copy of an absorp_analytic_fn object
        """
        a = absorp_analytic_fn()
        (a.A1, a.A2, a.A3, a.a1, a.a3, a.d) = (
            self.A1, self.A2, self.A3, self.a1, self.a3, self.d)
        return a

    def run(self, z):
        """
        Calculates absorption at a given depth z, where z=0 is the start of the
        layer.
        """

        if 'ndarray' in str(type(z)) and z.ndim > 0:
            part1 = self.A1[:, None] * np.exp(self.a1[:, None] * z[None, :])
            part2 = self.A2[:, None] * np.exp(-self.a1[:, None] * z[None, :])
            part3 = self.A3[:, None] * np.exp(1j * self.a3[:, None] * z[None, :])
            part4 = np.conj(self.A3[:, None]) * np.exp(-1j * self.a3[:, None] * z[None, :])

            part1[self.A1 < 1e-100, :] = 0
            return (part1 + part2 + part3 + part4)

        else:
            return (self.A1 * np.exp(self.a1 * z) + self.A2 * np.exp(-self.a1 * z)
                    + self.A3 * np.exp(1j * self.a3 * z) + np.conj(self.A3) * np.exp(-1j * self.a3 * z))



    def flip(self):
        """
        Flip the function front-to-back, to describe a(d-z) instead of a(z),
        where d is layer thickness.
        """
        expn = np.exp(self.a1 * self.d)
        #expn[expn > 1e100] = 1e100
        newA1 = self.A2 * np.exp(-self.a1 * self.d)
        newA1[self.A2 == 0] = 0
        newA2 = self.A1 * expn
        newA2[self.A1 == 0] = 0
        self.A1, self.A2 = newA1, newA2
        self.A3 = np.conj(self.A3 * np.exp(1j * self.a3 * self.d))
        return self

    def scale(self, factor):
        """
        multiplies the absorption at each point by "factor".
        """
        self.A1 *= factor
        self.A2 *= factor
        self.A3 *= factor
        self.A1[np.isnan(self.A1)] = 0
        self.A2[np.isnan(self.A2)] = 0
        self.A3[np.isnan(self.A3)] = 0
        return self

    def add(self, b):
        """
        adds another compatible absorption analytical function
        """

        if np.all(b.a1 != self.a1) or np.all(b.a3 != self.a3):
            raise ValueError('Incompatible absorption analytical functions!')

        self.A1 += b.A1
        self.A2 += b.A2
        self.A3 += b.A3
        return self


def absorp_in_each_layer(coh_tmm_data):
    """
    An array listing what proportion of light is absorbed in each layer.

    Assumes the final layer eventually absorbs all transmitted light.

    Assumes the initial layer eventually absorbs all reflected light.

    Entries of array should sum to 1.

    coh_tmm_data is output of coh_tmm()
    """
    num_layers = len(coh_tmm_data['d_list'])
    num_lam_vec = len(coh_tmm_data['lam_vac'])
    power_entering_each_layer = np.zeros((num_layers, num_lam_vec))
    power_entering_each_layer[0] = 1
    power_entering_each_layer[1] = coh_tmm_data['power_entering']
    power_entering_each_layer[-1] = coh_tmm_data['T']
    for i in range(2, num_layers - 1):
        power_entering_each_layer[i] = position_resolved(i, 0, coh_tmm_data)['poyn']
    final_answer = np.zeros((num_layers, num_lam_vec))
    final_answer[0:-1] = -np.diff(power_entering_each_layer, axis=0)
    final_answer[-1] = power_entering_each_layer[-1]


    final_answer[final_answer < 0] = 0

    return final_answer


def inc_group_layers(n_list, d_list, c_list):
    """
    Helper function for inc_tmm. Groups and sorts layer information.

    See coh_tmm for definitions of n_list, d_list.

    c_list is "coherency list". Each entry should be 'i' for incoherent or 'c'
    for 'coherent'.

    A "stack" is a group of one or more consecutive coherent layers. A "stack
    index" labels the stacks 0,1,2,.... The "within-stack index" counts the
    coherent layers within the stack 1,2,3... [index 0 is the incoherent layer
    before the stack starts]

    An "incoherent layer index" labels the incoherent layers 0,1,2,...

    An "alllayer index" labels all layers (all elements of d_list) 0,1,2,...

    Returns info about how the layers relate:

    * stack_d_list[i] = list of thicknesses of each coherent layer in the i'th
      stack, plus starting and ending with "inf"
    * stack_n_list[i] = list of refractive index of each coherent layer in the
      i'th stack, plus the two surrounding incoherent layers
    * all_from_inc[i] = j means that the layer with incoherent index i has
      alllayer index j
    * inc_from_all[i] = j means that the layer with alllayer index i has
      incoherent index j. If j = nan then the layer is coherent.
    * all_from_stack[i1][i2] = j means that the layer with stack index i1 and
      within-stack index i2 has alllayer index j
    * stack_from_all[i] = [j1 j2] means that the layer with alllayer index i is
      part of stack j1 with withinstack-index j2. If stack_from_all[i] = nan
      then the layer is incoherent
    * inc_from_stack[i] = j means that the i'th stack comes after the layer
      with incoherent index j, and before the layer with incoherent index j+1.
    * stack_from_inc[i] = j means that the layer with incoherent index i comes
      immediately after the j'th stack. If j=nan, it is not immediately
      following a stack.
    * num_stacks = number of stacks
    * num_inc_layers = number of incoherent layers
    * num_layers = number of layers total
    """

    if (d_list.ndim != 1):
        raise ValueError("Problem with n_list or d_list!")
    if (d_list[0] != np.inf) or (d_list[-1] != np.inf):
        raise ValueError('d_list must start and end with inf!')
    if (c_list[0] != 'i') or (c_list[-1] != 'i'):
        raise ValueError('c_list should start and end with "i"')
    if not len(n_list) == d_list.size == len(c_list):
        raise ValueError('List sizes do not match!')
    inc_index = 0
    stack_index = 0
    stack_d_list = []
    stack_n_list = []
    all_from_inc = []
    inc_from_all = []
    all_from_stack = []
    stack_from_all = []
    inc_from_stack = []
    stack_from_inc = []
    stack_in_progress = False
    for alllayer_index in range(len(n_list)):
        if c_list[alllayer_index] == 'c':  # coherent layer
            inc_from_all.append(np.nan)
            if not stack_in_progress:  # this layer is starting new stack
                stack_in_progress = True
                ongoing_stack_d_list = [np.inf, d_list[alllayer_index]]
                ongoing_stack_n_list = [n_list[alllayer_index - 1],
                                        n_list[alllayer_index]]
                stack_from_all.append([stack_index, 1])
                all_from_stack.append([alllayer_index - 1, alllayer_index])
                inc_from_stack.append(inc_index - 1)
                within_stack_index = 1
            else:  # another coherent layer in the same stack
                ongoing_stack_d_list.append(d_list[alllayer_index])
                ongoing_stack_n_list.append(n_list[alllayer_index])
                within_stack_index += 1
                stack_from_all.append([stack_index, within_stack_index])
                all_from_stack[-1].append(alllayer_index)
        elif c_list[alllayer_index] == 'i':  # incoherent layer
            stack_from_all.append(np.nan)
            inc_from_all.append(inc_index)
            all_from_inc.append(alllayer_index)
            if not stack_in_progress:  # previous layer was also incoherent
                stack_from_inc.append(np.nan)
            else:  # previous layer was coherent
                stack_in_progress = False
                stack_from_inc.append(stack_index)
                ongoing_stack_d_list.append(np.inf)
                stack_d_list.append(ongoing_stack_d_list)
                ongoing_stack_n_list.append(n_list[alllayer_index])
                stack_n_list.append(ongoing_stack_n_list)
                all_from_stack[-1].append(alllayer_index)
                stack_index += 1
            inc_index += 1
        else:
            raise ValueError("Error: c_list entries must be 'i' or 'c'!")
    return {'stack_d_list': stack_d_list,
            'stack_n_list': stack_n_list,
            'all_from_inc': all_from_inc,
            'inc_from_all': inc_from_all,
            'all_from_stack': all_from_stack,
            'stack_from_all': stack_from_all,
            'inc_from_stack': inc_from_stack,
            'stack_from_inc': stack_from_inc,
            'num_stacks': len(all_from_stack),
            'num_inc_layers': len(all_from_inc),
            'num_layers': len(n_list)}


def inc_tmm(pol, n_list, d_list, c_list, th_0, lam_vac):
    """
    This function is vectorized.
    Incoherent, or partly-incoherent-partly-coherent, transfer matrix method.

    See coh_tmm for definitions of pol, n_list, d_list, th_0, lam_vac.

    c_list is "coherency list". Each entry should be 'i' for incoherent or 'c'
    for 'coherent'.

    If an incoherent layer has real refractive index (no absorption), then its
    thickness doesn't affect the calculation results.

    See https://arxiv.org/abs/1603.02720 for physics background and some
    of the definitions.

    Outputs the following as a dictionary:

    * R--reflected wave power (as fraction of incident)
    * T--transmitted wave power (as fraction of incident)
    * VW_list-- n'th element is [V_n,W_n], the forward- and backward-traveling
      intensities, respectively, at the beginning of the n'th incoherent medium.
    * coh_tmm_data_list--n'th element is coh_tmm_data[n], the output of
      the coh_tmm program for the n'th "stack" (group of one or more
      consecutive coherent layers).
    * coh_tmm_bdata_list--n'th element is coh_tmm_bdata[n], the output of the
      coh_tmm program for the n'th stack, but with the layers of the stack
      in reverse order.
    * stackFB_list--n'th element is [F,B], where F is light traveling forward
      towards the n'th stack and B is light traveling backwards towards the n'th
      stack.
    * num_layers-- total number both coherent and incoherent.
    * power_entering_list--n'th element is the normalized Poynting vector
      crossing the interface into the n'th incoherent layer from the previous
      (coherent or incoherent) layer.
    * Plus, all the outputs of inc_group_layers

    """
    # Convert lists to numpy arrays if they're not already.
    n_list = np.array(n_list)
    d_list = np.array(d_list, dtype=float)

    # Input tests
    # =============================================================================
    #     if hasattr(n_list[0], 'size') and 'complex' not in str(type(n_list[0])):
    #         raise ValueError(
    #             'This function is not vectorized, yet; you need to run one calculation at a time (one wavelength, one angle).')
    # =============================================================================
    if any((np.real_if_close(n_list[0] * np.sin(th_0))).imag != 0):
        raise ValueError('Error in n0 or th0!')

    group_layers_data = inc_group_layers(n_list, d_list, c_list)
    num_inc_layers = group_layers_data['num_inc_layers']
    num_stacks = group_layers_data['num_stacks']
    stack_n_list = group_layers_data['stack_n_list']
    stack_d_list = group_layers_data['stack_d_list']
    all_from_stack = group_layers_data['all_from_stack']
    all_from_inc = group_layers_data['all_from_inc']
    all_from_stack = group_layers_data['all_from_stack']
    stack_from_inc = group_layers_data['stack_from_inc']
    inc_from_stack = group_layers_data['inc_from_stack']

    # th_list is a list with, for each layer, the angle that the light travels
    # through the layer. Computed with Snell's law. Note that the "angles" may be
    # complex!
    th_list = list_snell(n_list, th_0)

    # coh_tmm_data_list[i] is the output of coh_tmm for the i'th stack
    coh_tmm_data_list = []
    # coh_tmm_bdata_list[i] is the same stack as coh_tmm_data_list[i] but
    # with order of layers reversed

    coh_tmm_bdata_list = []
    for i in range(num_stacks):

        coh_tmm_data_list.append(coh_tmm(pol, stack_n_list[i],
                                         stack_d_list[i],
                                         th_list[all_from_stack[i][0]],
                                         lam_vac))
        coh_tmm_bdata_list.append(coh_tmm_reverse(pol, stack_n_list[i],
                                                  stack_d_list[i],
                                                  th_list[all_from_stack[i][0]],
                                                  lam_vac))

    # P_list[i] is fraction not absorbed in a single pass through i'th incoherent
    # layer.

    P_list = np.zeros((num_inc_layers, len(lam_vac)))
    for inc_index in range(1, num_inc_layers - 1):  # skip 0'th and last (infinite)
        i = all_from_inc[inc_index]
        P_list[inc_index] = np.exp(-4 * np.pi * d_list[i]
                                * (n_list[i] * np.cos(th_list[i])).imag / lam_vac)

    # For a very opaque layer, reset P to avoid divide-by-0 and similar
    # errors.
    P_list[P_list < 1e-30] = 1e-30

    # T_list[i,j] and R_list[i,j] are transmission and reflection powers,
    # respectively, coming from the i'th incoherent layer, going to the j'th
    # incoherent layer. Only need to calculate this when j=i+1 or j=i-1.
    # (2D array is overkill but helps avoid confusion.)
    # initialize these arrays
    T_list = np.zeros((num_inc_layers, num_inc_layers, len(lam_vac)))
    R_list = np.zeros((num_inc_layers, num_inc_layers, len(lam_vac)))
    for inc_index in range(num_inc_layers - 1):  # looking at interface i -> i+1
        alllayer_index = all_from_inc[inc_index]
        nextstack_index = stack_from_inc[inc_index + 1]
        if np.isnan(nextstack_index):  # next layer is incoherent
            R_list[inc_index, inc_index + 1] = (
                interface_R(pol, n_list[alllayer_index],
                            n_list[alllayer_index + 1],
                            th_list[alllayer_index],
                            th_list[alllayer_index + 1]))
            T_list[inc_index, inc_index + 1] = (
                interface_T(pol, n_list[alllayer_index],
                            n_list[alllayer_index + 1],
                            th_list[alllayer_index],
                            th_list[alllayer_index + 1]))
            R_list[inc_index + 1, inc_index] = (
                interface_R(pol, n_list[alllayer_index + 1],
                            n_list[alllayer_index],
                            th_list[alllayer_index + 1],
                            th_list[alllayer_index]))
            T_list[inc_index + 1, inc_index] = (
                interface_T(pol, n_list[alllayer_index + 1],
                            n_list[alllayer_index],
                            th_list[alllayer_index + 1],
                            th_list[alllayer_index]))
        else:  # next layer is coherent
            R_list[inc_index, inc_index + 1] = (
                coh_tmm_data_list[nextstack_index]['R'])
            T_list[inc_index, inc_index + 1] = (
                coh_tmm_data_list[nextstack_index]['T'])
            R_list[inc_index + 1, inc_index] = (
                coh_tmm_bdata_list[nextstack_index]['R'])
            T_list[inc_index + 1, inc_index] = (
                coh_tmm_bdata_list[nextstack_index]['T'])

    # L is the transfer matrix from the i'th to (i+1)st incoherent layer, see
    # manual

    # For a very opaque layer, reset T to avoid divide-by-0 and similar
    # errors.
    T_list[T_list < 1e-30] = 1e-30

    L_list = [np.nan]  # L_0 is not defined because 0'th layer has no beginning.
    Ltilde = (np.array([[np.ones_like(lam_vac), -R_list[1, 0]],
                     [R_list[0, 1],
                      T_list[1, 0] * T_list[0, 1] - R_list[1, 0] * R_list[0, 1]]])
              / T_list[0, 1]).transpose(2, 0, 1)

    # Ltilde = Ltilde.transpose(2,0,1)

    for i in range(1, num_inc_layers - 1):
        L = np.matmul(
            np.array([[1 / P_list[i], np.zeros(len(lam_vac))], [np.zeros(len(lam_vac)), P_list[i]]]).transpose(2, 0, 1),
            np.array([[np.ones_like(lam_vac), -R_list[i + 1, i]],
                   [R_list[i, i + 1],
                    T_list[i + 1, i] * T_list[i, i + 1] - R_list[i + 1, i] * R_list[i, i + 1]]]).transpose(2, 0, 1)
        ) / T_list[i, i + 1][:, None, None]
        L_list.append(L)
        Ltilde = np.matmul(Ltilde, L)

    T = 1 / Ltilde[:, 0, 0]
    R = Ltilde[:, 1, 0] / Ltilde[:, 0, 0]

    # VW_list[n] = [V_n, W_n], the forward- and backward-moving intensities
    # at the beginning of the n'th incoherent layer. VW_list[0] is undefined
    # because 0'th layer has no beginning.

    VW_list = np.zeros((num_inc_layers, 2, len(lam_vac)))
    VW_list[0, :, :] = np.nan
    VW = np.array([T, np.zeros(len(lam_vac))])
    VW_list[-1] = VW

    for i in range(num_inc_layers - 2, 0, -1):
        VW = np.matmul(L_list[i], VW_list[i+1].T[:, :, None])
        VW_list[i, :, :] = VW.transpose()


    # stackFB_list[n]=[F,B] means that F is light traveling forward towards n'th
    # stack and B is light traveling backwards towards n'th stack.
    # Reminder: inc_from_stack[i] = j means that the i'th stack comes after the
    # layer with incoherent index j.

    stackFB_list = []
    for stack_index, prev_inc_index in enumerate(inc_from_stack):
        if prev_inc_index == 0:  # stack starts right after semi-infinite layer.
            F = np.ones_like(lam_vac)
        else:
            F = VW_list[prev_inc_index, 0] * P_list[prev_inc_index]
        B = VW_list[prev_inc_index + 1, 1]
        stackFB_list.append([F, B])

    # power_entering_list[i] is the normalized Poynting vector crossing the
    # interface into the i'th incoherent layer from the previous (coherent or
    # incoherent) layer. See manual.
    power_entering_list = [np.ones_like(lam_vac)]  # "1" by convention for infinite 0th layer.
    for i in range(1, num_inc_layers):
        prev_stack_index = stack_from_inc[i]
        if np.isnan(prev_stack_index):
            # case where this layer directly follows another incoherent layer
            if i == 1:  # special case because VW_list[0] & A_list[0] are undefined
                power_entering_list.append(T_list[0, 1]
                                           - VW_list[1, 1] * T_list[1, 0])
            else:
                power_entering_list.append(
                    VW_list[i - 1][0] * P_list[i - 1] * T_list[i - 1, i]
                    - VW_list[i][1] * T_list[i, i - 1])
        else:  # case where this layer follows a coherent stack
            power_entering_list.append(
                stackFB_list[prev_stack_index][0] *
                coh_tmm_data_list[prev_stack_index]['T']
                - stackFB_list[prev_stack_index][1] *
                coh_tmm_bdata_list[prev_stack_index]['power_entering'])

    if len(stackFB_list) > 0:
        stackFB_list_ans = np.stack(stackFB_list).transpose(2, 0, 1)
    else:
        stackFB_list_ans = []
    #('VWlist', VW_list.transpose(2, 0, 1))
    ans = {'T': T, 'R': R, 'VW_list': VW_list.transpose(2, 0, 1),
           'coh_tmm_data_list': coh_tmm_data_list,
           'coh_tmm_bdata_list': coh_tmm_bdata_list,
           'stackFB_list': stackFB_list_ans,
           'power_entering_list': np.stack(power_entering_list).T}
    ans.update(group_layers_data)
    return ans


def inc_absorp_in_each_layer(inc_data):
    """
    A list saying what proportion of light is absorbed in each layer.

    Assumes all reflected light is eventually absorbed in the 0'th medium, and
    all transmitted light is eventually absorbed in the final medium.

    Returns a list [layer0absorp, layer1absorp, ...]. Entries should sum to 1.

    inc_data is output of inc_tmm()
    """
    # Reminder: inc_from_stack[i] = j means that the i'th stack comes after the
    # layer with incoherent index j.
    # Reminder: stack_from_inc[i] = j means that the layer
    # with incoherent index i comes immediately after the j'th stack (or j=nan
    # if it's not immediately following a stack).

    stack_from_inc = inc_data['stack_from_inc']
    power_entering_list = inc_data['power_entering_list'].T
    # stackFB_list[n]=[F,B] means that F is light traveling forward towards n'th
    # stack and B is light traveling backwards towards n'th stack.
    stackFB_list = inc_data['stackFB_list']
    absorp_list = []

    # loop through incoherent layers, excluding the final layer
    for i, power_entering in enumerate(power_entering_list[:-1]):
        if np.isnan(stack_from_inc[i + 1]):
            # case that incoherent layer i is right before another incoherent layer
            absorp_list.append(power_entering_list[i] - power_entering_list[i + 1])
        else:  # incoherent layer i is immediately before a coherent stack
            j = stack_from_inc[i + 1]
            coh_tmm_data = inc_data['coh_tmm_data_list'][j]
            coh_tmm_bdata = inc_data['coh_tmm_bdata_list'][j]
            # First, power in the incoherent layer...
            power_exiting = (
                    stackFB_list[:, j, 0] * coh_tmm_data['power_entering']
                    - stackFB_list[:, j, 1] * coh_tmm_bdata['T'])
            absorp_list.append(power_entering_list[i] - power_exiting)
            # Next, power in the coherent stack...
            stack_absorp = ((stackFB_list[:, j, 0] *
                             absorp_in_each_layer(coh_tmm_data))[1:-1]
                            + (stackFB_list[:, j, 1] *
                               absorp_in_each_layer(coh_tmm_bdata))[-2:0:-1])
            absorp_list.extend(stack_absorp)
    # final semi-infinite layer
    absorp_list.append(inc_data['T'])

    absorp_list = np.array(absorp_list)
    absorp_list[absorp_list < 0] = 0

    return absorp_list


def inc_find_absorp_analytic_fn(layer, inc_data):
    """
    Outputs an absorp_analytic_fn object for a coherent layer within a
    partly-incoherent stack.

    inc_data is output of inc_tmm()
    """
    j = inc_data['stack_from_all'][layer]
    if np.any(np.isnan(j)):
        raise ValueError('layer must be coherent for this function!')
    [stackindex, withinstackindex] = j
    forwardfunc = absorp_analytic_fn()
    forwardfunc.fill_in(inc_data['coh_tmm_data_list'][stackindex],
                        withinstackindex)
    forwardfunc.scale(inc_data['stackFB_list'][:, stackindex, 0])
    backfunc = absorp_analytic_fn()
    backfunc.fill_in(inc_data['coh_tmm_bdata_list'][stackindex],
                     -1 - withinstackindex)

    backfunc.scale(inc_data['stackFB_list'][:, stackindex, 1])

    backfunc.flip()

    return forwardfunc.add(backfunc)


def inc_position_resolved(layer, dist, inc_tmm_data, coherency_list, alphas, zero_threshold=1e-6):
    """
    This function is vectorized. Analogous to position_resolved, but
    for layers (incoherent or coherent) in (partly) incoherent stacks.
    This is a new function, not from Steven Byrnes' tmm package.
    Starting with output of inc_tmm(), calculate the Poynting vector
    and absorbed energy density a distance "dist" into layer number "layer"
    """

    layers = list(set(layer)) # unique layer indices
    A_per_layer = np.array(inc_absorp_in_each_layer(inc_tmm_data))
    fraction_reaching = 1 - np.cumsum(A_per_layer, axis = 0)
    A_local = np.zeros((len(alphas[0]), len(dist)))
    for i, l in enumerate(layers):
        if coherency_list[l] == 'c':

            fn = inc_find_absorp_analytic_fn(l, inc_tmm_data)
            A_layer = fn.run(dist[layer == l])

        else:
            A_layer = beer_lambert(alphas[l] * 1e9, fraction_reaching[i], dist[layer == l] * 1e-9)

        A_layer[fraction_reaching[i] < zero_threshold, :] = 0
        A_local[:, layer == l] = A_layer

    return A_local


def beer_lambert(alphas, fraction, dist):
    """
    Calculates absorption profile according to the Beer-Lambert law given alphas (in m-1)
    and a vector of distance into the layer (in m) and the fraction of incident light
    reaching the front of the layer. This is used to calculate the absorption profile in
    incoherent layers within (partly) incoherent stacks. Vectorized over wavelengths.
    """

    expn = np.exp(- alphas[:, None] * dist[None,:])

    output = fraction[:, None]*alphas[:, None]*expn

    return output/1e9

