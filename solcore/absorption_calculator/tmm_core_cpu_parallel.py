"""
This is implementation of the coh_tmm function using numba jit (CPU parallelization)
For large number of calculations (e.g. s and p polarization, 10,000 wavelengths x angles)
the speed increase over tmm_core_vec implementation is signficant, execution time @ 35% 
(tested with 24 cores)
if further set detailed = False, will skip calculations of power_entering 
and vw_list, and the function will be faster still, execution time roughly 85% x 35% = 30%
Note however that the first time the function executes, the jit will consume a lot of time,
so the speed increase is only observed after an initial run
"""

import numpy as np
import sys
from numba import jit

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

@jit(nopython=True, parallel=True)
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
    # Note: scipy.arcsin is deprecated, hence scimath replacement
    return np.arcsin(n_list[0] * np.sin(th_0) / n_list)

@jit(nopython=True, parallel=True)
def power_entering_from_r(pol, r, n_i, cos_th_i):
    """
    Calculate the power entering the first interface of the stack, starting with
    reflection amplitude r. Normally this equals 1-R, but in the unusual case
    that n_i is not real, it can be a bit different than 1-R. See manual.

    n_i is refractive index of incident medium.

    th_i is (complex) propegation angle through incident medium
    (in radians, where 0=normal). "th" stands for "theta".
    """
    if (pol == 's'):
        return ((n_i * cos_th_i * (1 + np.conj(r)) * (1 - r)).real
                / (n_i * cos_th_i).real)
    elif (pol == 'p'):
        return ((n_i * np.conj(cos_th_i) * (1 + r) * (1 - np.conj(r))).real
                / (n_i * np.conj(cos_th_i)).real)
    else:
        raise ValueError("Polarization must be 's' or 'p'")

@jit(nopython=True, parallel=True)
def construct_r_t_list(n_list, cos_th_list, pol, r_list, t_list):
    if pol == 's':
        t_list[:,:-1] = np.transpose((2 * n_list[:-1] * cos_th_list[:-1]) /
                (n_list[:-1] * cos_th_list[:-1] + n_list[1:] * cos_th_list[1:]))
        r_list[:,:-1] = np.transpose((n_list[:-1] * cos_th_list[:-1] - n_list[1:] * cos_th_list[1:]) /
                (n_list[:-1] * cos_th_list[:-1] + n_list[1:] * cos_th_list[1:]))
    else:
        t_list[:,:-1] = np.transpose((2 * n_list[:-1] * cos_th_list[:-1]) /
                (n_list[1:] * cos_th_list[:-1] + n_list[:-1] * cos_th_list[1:]))
        r_list[:,:-1] = np.transpose((n_list[1:] * cos_th_list[:-1] - n_list[:-1] * cos_th_list[1:]) /
                (n_list[1:] * cos_th_list[:-1] + n_list[:-1] * cos_th_list[1:]))

@jit(nopython=True, parallel=True)
def construct_M_list(num_layers, delta, r_list, t_list, M_list):
    for i in range(0, num_layers - 1):
        exp_ = np.exp(-1j * delta[i])
        d = 1/t_list[:,i]
        M_list[i,:,0,0] = exp_*d
        M_list[i,:,0,1] = exp_*d*r_list[:, i]
        M_list[i,:,1,0] = d/exp_*r_list[:, i]
        M_list[i,:,1,1] = d/exp_

@jit(nopython=True, parallel=True)
def construct_Mtilde(num_layers, M_list, Mtilde):
    for i in range(0, num_layers - 1):
        M_00 = Mtilde[:,0,0]*M_list[i,:,0,0] + Mtilde[:,0,1]*M_list[i,:,1,0]
        M_01 = Mtilde[:,0,0]*M_list[i,:,0,1] + Mtilde[:,0,1]*M_list[i,:,1,1]
        M_10 = Mtilde[:,1,0]*M_list[i,:,0,0] + Mtilde[:,1,1]*M_list[i,:,1,0]
        M_11 = Mtilde[:,1,0]*M_list[i,:,0,1] + Mtilde[:,1,1]*M_list[i,:,1,1]
        
        Mtilde[:,0,0] = M_00
        Mtilde[:,0,1] = M_01
        Mtilde[:,1,0] = M_10
        Mtilde[:,1,1] = M_11

@jit(nopython=True, parallel=True)
def construct_r_t_R_T(Mtilde, pol, n_i, n_f, cos_th_i, cos_th_f):
    r = Mtilde[:, 1, 0] / Mtilde[:, 0, 0]
    t = np.ones_like(Mtilde[:, 0, 0]) / Mtilde[:, 0, 0]
    R = np.abs(r) ** 2
    if (pol == 's'):
        T = np.abs(t ** 2) * (((n_f * cos_th_f).real) / (n_i * cos_th_i).real)
    else:
        T = np.abs(t ** 2) * (((n_f * np.conj(cos_th_f)).real) /
                              (n_i * np.conj(cos_th_i)).real)
    return r, t, R, T

@jit(nopython=True, parallel=True)
def construct_vw(num_layers, M_list, vw, vw_list):
    for i in range(num_layers - 2, 0, -1):
        M_00 = M_list[i,:,0,0]*vw[:,0,0] + M_list[i,:,0,1]*vw[:,1,0]
        M_01 = M_list[i,:,0,0]*vw[:,0,1] + M_list[i,:,0,1]*vw[:,1,1]
        M_10 = M_list[i,:,1,0]*vw[:,0,0] + M_list[i,:,1,1]*vw[:,1,0]
        M_11 = M_list[i,:,1,0]*vw[:,0,1] + M_list[i,:,1,1]*vw[:,1,1]
        
        vw[:,0,0] = M_00
        vw[:,0,1] = M_01
        vw[:,1,0] = M_10
        vw[:,1,1] = M_11
        vw_list[i, :, :] = vw[:, :, 1]
    vw_list[-1, :, 1] = 0
        
def coh_tmm(pol, n_list, d_list, th_0, lam_vac, detailed=True):
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
    d_list[0] = 0.0
    d_list[-1] = 0.0
    d_list = d_list.astype(complex)

    # input tests
    # if hasattr(th_0, 'size') and th_0.size > 1 and th_0.size != lam_vac.size:
    #     raise ValueError('This function is not vectorized for angles; you need to run one angle calculation at a time.')
    # if n_list.shape[0] != d_list.shape[0]:
    #     raise ValueError("Problem with n_list or d_list!")
    # if (d_list[0] != np.inf) or (d_list[-1] != np.inf):
    #     raise ValueError('d_list must start and end with inf!')
    # if any(abs((n_list[0] * np.sin(th_0)).imag) > 100 * EPSILON):
    #     raise ValueError('Error in n0 or th0!')
    if hasattr(th_0, 'size'):
        th_0 = np.array(th_0)
    num_layers = n_list.shape[0]
    num_wl = n_list.shape[1]

    # th_list is a list with, for each layer, the angle that the light travels
    # through the layer. Computed with Snell's law. Note that the "angles" may be
    # complex!
    th_list = list_snell(n_list, th_0)
    cos_th_0 = np.cos(th_0)
    cos_th_list = np.cos(th_list)

    # kz is the z-component of (complex) angular wavevector for forward-moving
    # wave. Positive imaginary part means decaying.
    kz_list = 2 * np.pi * n_list * cos_th_list / lam_vac

    # delta is the total phase accrued by traveling through a given layer.
    # ignore warning about inf multiplication
    olderr = np.seterr(invalid='ignore')
    delta = kz_list * d_list
    np.seterr(**olderr)

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
    t_list = np.zeros((num_wl, num_layers), dtype=complex)
    r_list = np.zeros((num_wl, num_layers), dtype=complex)

    construct_r_t_list(n_list, cos_th_list, pol, r_list, t_list)

    # At the interface between the (n-1)st and nth material, let v_n be the
    # amplitude of the wave on the nth side heading forwards (away from the
    # boundary), and let w_n be the amplitude on the nth side heading backwards
    # (towards the boundary). Then (v_n,w_n) = M_n (v_{n+1},w_{n+1}). M_n is
    # M_list[n]. M_0 and M_{num_layers-1} are not defined.
    # My M is a bit different than Sernelius's, but Mtilde is the same.

    M_list = np.zeros((num_layers, num_wl, 2, 2), dtype=complex)
    construct_M_list(num_layers, delta, r_list, t_list, M_list)

    Mtilde = make_2x2_array(np.ones_like(delta[0]), np.zeros_like(delta[0]), np.zeros_like(delta[0]),
                            np.ones_like(delta[0]), dtype=complex)
    construct_Mtilde(num_layers, M_list, Mtilde)

    r, t, R, T = construct_r_t_R_T(Mtilde, pol, n_list[0], n_list[-1], cos_th_0, cos_th_list[-1])

    # vw_list[n] = [v_n, w_n]. v_0 and w_0 are undefined because the 0th medium
    # has no left interface.
    vw_list = []
    power_entering = []
    if detailed:
        vw_list = np.zeros((num_layers, num_wl, 2), dtype=complex)
        vw = np.zeros((num_wl, 2, 2), dtype=complex)
        vw[:, 0, 0] = t
        vw[:, 0, 1] = t
        vw_list[-1] = vw[:, 0, :]
        construct_vw(num_layers, M_list, vw, vw_list)

        power_entering = power_entering_from_r(
            pol, r, n_list[0], cos_th_0)

    return {'r': r, 't': t, 'R': R, 'T': T, 'power_entering': power_entering,
            'vw_list': vw_list, 'kz_list': kz_list, 'th_list': th_list,
            'pol': pol, 'n_list': n_list, 'd_list': d_list, 'th_0': th_0,
            'lam_vac': lam_vac}
