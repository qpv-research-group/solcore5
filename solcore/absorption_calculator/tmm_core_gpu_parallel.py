"""
This is implementation of the coh_tmm function using numba cuda (GPU parallelization)
For large number of calculations (e.g. s and p polarization, 10,000 wavelengths x angles)
the speed increase over tmm_core_vec implementation is significant, execution time @ 12% 
(tested with NVIDIA GeForce RTX 4060)
if further set detailed = False, will skip calculations of power_entering 
and vw_list, and the function will be faster still, execution time roughly 40% x 12% = 4.8%
Note however that the first time the function executes, the compilation will consume a lot of time,
so the speed increase is only observed after an initial run
"""

import numpy as np
import sys
from numba import cuda
import math
import cmath

EPSILON = sys.float_info.epsilon  # typical floating-point calculation error


def list_snell(d_n_list, th_0, lam_vac, d_list, d_th_list, d_cos_th_list, d_delta, threadsperblock, catch_error):
    d_kz_list = cuda.device_array((d_n_list.shape[0],d_n_list.shape[1]), dtype=complex)
    succeed = False
    if catch_error:
        for _ in range(4):
            try:
                blockspergrid = math.ceil(d_n_list.shape[0]*d_n_list.shape[1] / threadsperblock)
                list_snell_kernel_function[blockspergrid, threadsperblock](d_n_list, th_0, lam_vac, d_list, d_th_list, d_cos_th_list, d_kz_list, d_delta)
                succeed = True
            except cuda.cudadrv.driver.CudaAPIError as e:
                threadsperblock //= 2
    else:
        blockspergrid = math.ceil(d_n_list.shape[0]*d_n_list.shape[1] / threadsperblock)
        list_snell_kernel_function[blockspergrid, threadsperblock](d_n_list, th_0, lam_vac, d_list, d_th_list, d_cos_th_list, d_kz_list, d_delta)
        succeed = True

    th_list = d_th_list.copy_to_host()
    kz_list = d_kz_list.copy_to_host()
    return th_list, kz_list, threadsperblock, succeed

@cuda.jit
def list_snell_kernel_function(n_list, th_0, lam_vac, d_list, th_list, cos_th_list, kz_list, delta):
    num_of_layers = n_list.shape[0]
    num_of_rows = n_list.shape[1]
    pos = cuda.grid(1)
    if pos < num_of_layers*num_of_rows:
        pos1 = int(pos/num_of_layers)
        pos0 = pos - pos1*num_of_layers
        th_list[pos0,pos1] = np.arcsin(n_list[0,pos1] * np.sin(th_0) / n_list[pos0,pos1])
        cos_th_list[pos0,pos1] = np.cos(th_list[pos0,pos1])
        kz_list[pos0,pos1] = 2 * np.pi * n_list[pos0,pos1] * cos_th_list[pos0,pos1] / lam_vac[pos1]
        delta[pos0,pos1] = kz_list[pos0,pos1] * d_list[pos0,0]
        if delta[pos0,pos1].imag > 100:
            delta[pos0,pos1] = delta[pos0,pos1].real + 100j

def interface_r(d_n_list, d_cos_th_list, d_r_list, threadsperblock, catch_error):
    succeed = False
    if catch_error:
        for _ in range(4):
            try:
                blockspergrid = math.ceil((d_n_list.shape[0]-1)*d_n_list.shape[1]*2 / threadsperblock) 
                interface_r_kernel_function[blockspergrid, threadsperblock](d_n_list, d_cos_th_list, d_r_list)
                succeed = True
            except cuda.cudadrv.driver.CudaAPIError as e:
                threadsperblock //= 2
    else:
        blockspergrid = math.ceil((d_n_list.shape[0]-1)*d_n_list.shape[1]*2 / threadsperblock) 
        interface_r_kernel_function[blockspergrid, threadsperblock](d_n_list, d_cos_th_list, d_r_list)
        succeed = True
    return threadsperblock, succeed

@cuda.jit
def interface_r_kernel_function(n_list, cos_th_list, r_list):
    num_of_layers = n_list.shape[0]
    num_of_rows = n_list.shape[1]
    product = (num_of_layers-1)*num_of_rows
    pos = cuda.grid(1)
    if pos < product*2:
        pol = int(pos/product)
        remainder = pos - pol*product
        pos1 = int(remainder/(num_of_layers-1))
        pos0 = remainder - pos1*(num_of_layers-1)
        if pol==0: #s
            r_list[pos1,pos0,pol] = ((n_list[pos0,pos1] * cos_th_list[pos0,pos1] - n_list[pos0+1,pos1] * cos_th_list[pos0+1,pos1]) /
                    (n_list[pos0,pos1] * cos_th_list[pos0,pos1] + n_list[pos0+1,pos1] * cos_th_list[pos0+1,pos1]))
        else:
            r_list[pos1,pos0,pol] = ((n_list[pos0+1,pos1] * cos_th_list[pos0,pos1] - n_list[pos0,pos1] * cos_th_list[pos0+1,pos1]) /
                    (n_list[pos0+1,pos1] * cos_th_list[pos0,pos1] + n_list[pos0,pos1] * cos_th_list[pos0+1,pos1]))

def interface_t(d_n_list, d_cos_th_list, d_t_list, threadsperblock, catch_error):
    succeed = False
    if catch_error:
        for _ in range(4):
            try:
                blockspergrid = math.ceil((d_n_list.shape[0]-1)*d_n_list.shape[1]*2 / threadsperblock) 
                interface_t_kernel_function[blockspergrid, threadsperblock](d_n_list, d_cos_th_list, d_t_list)
                succeed = True
            except cuda.cudadrv.driver.CudaAPIError as e:
                threadsperblock //= 2
    else:
        blockspergrid = math.ceil((d_n_list.shape[0]-1)*d_n_list.shape[1]*2 / threadsperblock) 
        interface_t_kernel_function[blockspergrid, threadsperblock](d_n_list, d_cos_th_list, d_t_list)
        succeed = True
    return threadsperblock, succeed

@cuda.jit
def interface_t_kernel_function(n_list, cos_th_list, t_list):
    num_of_layers = n_list.shape[0]
    num_of_rows = n_list.shape[1]
    product = (num_of_layers-1)*num_of_rows
    pos = cuda.grid(1)
    if pos < product*2:
        pol = int(pos/product)
        remainder = pos - pol*product
        pos1 = int(remainder/(num_of_layers-1))
        pos0 = remainder - pos1*(num_of_layers-1)
        if pol==0: #s
            t_list[pos1,pos0,pol] = ((2 * n_list[pos0,pos1] * cos_th_list[pos0,pos1] ) /
                    (n_list[pos0,pos1] * cos_th_list[pos0,pos1] + n_list[pos0+1,pos1] * cos_th_list[pos0+1,pos1]))
        else:
            t_list[pos1,pos0,pol] = ((2 * n_list[pos0,pos1] * cos_th_list[pos0,pos1] ) /
                    (n_list[pos0+1,pos1] * cos_th_list[pos0,pos1] + n_list[pos0,pos1] * cos_th_list[pos0+1,pos1]))

def make_M_list(d_delta, d_t_list, d_r_list, d_M_list, threadsperblock, catch_error):
    succeed = False
    if catch_error:
        for _ in range(4):
            try:
                blockspergrid = math.ceil((d_M_list.shape[0]-1)*d_M_list.shape[1]*2 / threadsperblock) 
                make_M_list_kernel_function[blockspergrid, threadsperblock](d_delta, d_t_list, d_r_list, d_M_list)
                succeed = True
            except cuda.cudadrv.driver.CudaAPIError as e:
                threadsperblock //= 2
    else:
        blockspergrid = math.ceil((d_M_list.shape[0]-1)*d_M_list.shape[1]*2 / threadsperblock) 
        make_M_list_kernel_function[blockspergrid, threadsperblock](d_delta, d_t_list, d_r_list, d_M_list)
        succeed = True
    return threadsperblock, succeed
    
@cuda.jit
def make_M_list_kernel_function(delta, t_list, r_list, M_list):
    num_of_layers = M_list.shape[0]
    num_of_rows = M_list.shape[1]
    product = (num_of_layers-1)*num_of_rows
    pos = cuda.grid(1)
    if pos < product*2:
        pol = int(pos/product)
        remainder = pos - pol*product
        pos1 = int(remainder/(num_of_layers-1))
        pos0 = remainder - pos1*(num_of_layers-1)
        exp_ = 1.0
        if pos0 > 0:
            exp_ = cmath.exp(-1j * delta[pos0,pos1])
        d = (1 / t_list[pos1,pos0,pol])
        M_list[pos0,pos1,0,0,pol] = exp_*d
        M_list[pos0,pos1,0,1,pol] = exp_*d*r_list[pos1,pos0,pol]
        M_list[pos0,pos1,1,0,pol] = d/exp_*r_list[pos1,pos0,pol]
        M_list[pos0,pos1,1,1,pol] = d/exp_

def make_Mtilde(d_M_list, d_Mtilde, threadsperblock, catch_error):
    succeed = False
    if catch_error:
        for _ in range(4):
            try:
                blockspergrid = math.ceil(d_Mtilde.shape[0]*2 / threadsperblock) 
                make_Mtilde_kernel_function[blockspergrid, threadsperblock](d_M_list, d_Mtilde)
                succeed = True
            except cuda.cudadrv.driver.CudaAPIError as e:
                threadsperblock //= 2
    else:
        blockspergrid = math.ceil(d_Mtilde.shape[0]*2 / threadsperblock) 
        make_Mtilde_kernel_function[blockspergrid, threadsperblock](d_M_list, d_Mtilde)
        succeed = True
    return threadsperblock, succeed

@cuda.jit
def make_Mtilde_kernel_function(M_list, Mtilde):
    num_of_rows = Mtilde.shape[0]
    num_layers = M_list.shape[0]
    pos = cuda.grid(1)
    if pos < num_of_rows*2:
        pol = int(pos/num_of_rows)
        pos0 = pos - pol*num_of_rows
        Mtilde[pos0,0,0,pol] = 1.0
        Mtilde[pos0,0,1,pol] = 0.0
        Mtilde[pos0,1,0,pol] = 0.0
        Mtilde[pos0,1,1,pol] = 1.0
        for i in range(0, num_layers - 1):
            M00 = Mtilde[pos0,0,0,pol]*M_list[i,pos0,0,0,pol] + Mtilde[pos0,0,1,pol]*M_list[i,pos0,1,0,pol]
            M01 = Mtilde[pos0,0,0,pol]*M_list[i,pos0,0,1,pol] + Mtilde[pos0,0,1,pol]*M_list[i,pos0,1,1,pol]
            M10 = Mtilde[pos0,1,0,pol]*M_list[i,pos0,0,0,pol] + Mtilde[pos0,1,1,pol]*M_list[i,pos0,1,0,pol]
            M11 = Mtilde[pos0,1,0,pol]*M_list[i,pos0,0,1,pol] + Mtilde[pos0,1,1,pol]*M_list[i,pos0,1,1,pol]
            Mtilde[pos0,0,0,pol] = M00
            Mtilde[pos0,0,1,pol] = M01
            Mtilde[pos0,1,0,pol] = M10
            Mtilde[pos0,1,1,pol] = M11

def get_r_t_R_T(d_Mtilde, d_n_list, d_cos_th_list, threadsperblock, catch_error):
    d_r = cuda.device_array((d_Mtilde.shape[0],2), dtype=complex) # last index is polarization
    d_t = cuda.device_array((d_Mtilde.shape[0],2), dtype=complex)
    d_R = cuda.device_array((d_Mtilde.shape[0],2), dtype=float)
    d_T = cuda.device_array((d_Mtilde.shape[0],2), dtype=float)
    succeed = False
    if catch_error:
        for _ in range(4):
            try:
                blockspergrid = math.ceil(d_Mtilde.shape[0]*2 / threadsperblock) 
                get_r_t_R_T_kernel_function[blockspergrid, threadsperblock](d_Mtilde, d_n_list, d_cos_th_list, d_r, d_t, d_R, d_T)
                succeed = True
            except cuda.cudadrv.driver.CudaAPIError as e:
                threadsperblock //= 2
    else:
        blockspergrid = math.ceil(d_Mtilde.shape[0]*2 / threadsperblock) 
        get_r_t_R_T_kernel_function[blockspergrid, threadsperblock](d_Mtilde, d_n_list, d_cos_th_list, d_r, d_t, d_R, d_T)
        succeed = True
    
    r = d_r.copy_to_host()
    t = d_t.copy_to_host()
    R = d_R.copy_to_host()
    T = d_T.copy_to_host()
    return r, t, R, T, threadsperblock, succeed

@cuda.jit
def get_r_t_R_T_kernel_function(Mtilde, n_list, cos_th_list, r, t, R, T):
    num_of_rows = Mtilde.shape[0]
    pos = cuda.grid(1)
    if pos < num_of_rows*2:
       pol = int(pos/num_of_rows)
       pos0 = pos - pol*num_of_rows
       r[pos0,pol] = Mtilde[pos0, 1, 0, pol] / Mtilde[pos0, 0, 0, pol]
       t[pos0,pol] = 1.0 / Mtilde[pos0, 0, 0, pol]
       R[pos0,pol] = abs(r[pos0,pol]) ** 2
       if pol == 0: # s
           T[pos0,pol] = abs(t[pos0,pol] ** 2) * (((n_list[-1][pos0] * cos_th_list[-1][pos0]).real) / (n_list[0][pos0] * cos_th_list[0][pos0]).real)
       else:
           T[pos0,pol] = abs(t[pos0,pol] ** 2) * (((n_list[-1][pos0] * cos_th_list[-1][pos0].conjugate()).real) /
                              (n_list[0][pos0] * cos_th_list[0][pos0].conjugate()).real)
       
def build_vw(d_M_list,vw,vw_list, threadsperblock, catch_error):
    d_vw = cuda.to_device(vw)
    d_vw_list = cuda.to_device(vw_list)
    succeed = False
    if catch_error:
        for _ in range(4):
            try:
                blockspergrid = math.ceil(d_M_list.shape[1]*2 / threadsperblock) 
                build_vw_kernel_function[blockspergrid, threadsperblock](d_M_list,d_vw,d_vw_list)
                succeed = True
            except cuda.cudadrv.driver.CudaAPIError as e:
                threadsperblock //= 2
    else:
        blockspergrid = math.ceil(d_M_list.shape[1]*2 / threadsperblock) 
        build_vw_kernel_function[blockspergrid, threadsperblock](d_M_list,d_vw,d_vw_list)
        succeed = True

    revised_vw = d_vw.copy_to_host()
    revised_vw_list = d_vw_list.copy_to_host()
    return revised_vw, revised_vw_list, threadsperblock, succeed

@cuda.jit
def build_vw_kernel_function(M_list,vw,vw_list):
    num_layers = M_list.shape[0]
    num_of_rows = M_list.shape[1]
    pos = cuda.grid(1)
    if pos < num_of_rows*2:
        pol = int(pos/num_of_rows)
        pos0 = pos - pol*num_of_rows
        for i in range(num_layers - 2, 0, -1):
            M00 = M_list[i,pos0,0,0,pol]*vw[pos0,0,0,pol] + M_list[i,pos0,0,1,pol]*vw[pos0,1,0,pol]
            M01 = M_list[i,pos0,0,0,pol]*vw[pos0,0,1,pol] + M_list[i,pos0,0,1,pol]*vw[pos0,1,1,pol]
            M10 = M_list[i,pos0,1,0,pol]*vw[pos0,0,0,pol] + M_list[i,pos0,1,1,pol]*vw[pos0,1,0,pol]
            M11 = M_list[i,pos0,1,0,pol]*vw[pos0,0,1,pol] + M_list[i,pos0,1,1,pol]*vw[pos0,1,1,pol]
            vw[pos0,0,0,pol] = M00
            vw[pos0,0,1,pol] = M01
            vw[pos0,1,0,pol] = M10
            vw[pos0,1,1,pol] = M11
            vw_list[i,pos0,0,pol] = vw[pos0,0,1,pol]
            vw_list[i,pos0,1,pol] = vw[pos0,1,1,pol]
        vw_list[-1, pos0, 1,pol] = 0

def power_entering_from_r(r, n_i, cos_th_i):
    """
    Calculate the power entering the first interface of the stack, starting with
    reflection amplitude r. Normally this equals 1-R, but in the unusual case
    that n_i is not real, it can be a bit different than 1-R. See manual.

    n_i is refractive index of incident medium.

    th_i is (complex) propegation angle through incident medium
    (in radians, where 0=normal). "th" stands for "theta".
    """
    s_part = ((n_i * cos_th_i * (1 + np.conj(r[:,0])) * (1 - r[:,0])).real
                / (n_i * cos_th_i).real)
    p_part = ((n_i * np.conj(cos_th_i) * (1 + r[:,1]) * (1 - np.conj(r[:,1]))).real
                / (n_i * np.conj(cos_th_i)).real)
    return np.column_stack((s_part,p_part))

class Coh_tmm_GPU:
    def __init__(self):
        self.GPU_enabled = False
        if cuda.is_available():
            device = cuda.get_current_device()
            self.threadsperblock = [device.MAX_THREADS_PER_BLOCK]*7
            self.GPU_enabled = True
            self.catch_error = True

    def coh_tmm(self, n_list, d_list, th_0, lam_vac, detailed=True):
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
        d_n_list = cuda.to_device(n_list)
        d_th_list = cuda.device_array((n_list.shape[0],n_list.shape[1]), dtype=complex)
        d_cos_th_list = cuda.device_array((n_list.shape[0],n_list.shape[1]), dtype=complex)
        d_delta = cuda.device_array((n_list.shape[0],n_list.shape[1]), dtype=complex)
        d_list[0] = 0.0
        d_list[-1] = 0.0
        d_list = d_list.astype(complex)
        th_list, kz_list, threadsperblock, success = list_snell(d_n_list, th_0, lam_vac, d_list, d_th_list, d_cos_th_list, d_delta, self.threadsperblock[0], self.catch_error)
        self.threadsperblock[0] = threadsperblock

        olderr = np.seterr(invalid='ignore')
        np.seterr(**olderr)

        d_t_list = cuda.device_array((num_wl, num_layers, 2), dtype=complex) # last index is for both s and p
        d_r_list = cuda.device_array((num_wl, num_layers, 2), dtype=complex) # last index is for both s and p
        threadsperblock, success = interface_r(d_n_list, d_cos_th_list, d_r_list, self.threadsperblock[1], self.catch_error)
        if not success:
            return {}, False
        self.threadsperblock[1] = threadsperblock
        threadsperblock, success = interface_t(d_n_list, d_cos_th_list, d_t_list, self.threadsperblock[2], self.catch_error)
        if not success:
            return {}, False
        self.threadsperblock[2] = threadsperblock

        d_M_list = cuda.device_array((n_list.shape[0],n_list.shape[1],2,2,2), dtype=complex) # the last index is for both polarities
        threadsperblock, success = make_M_list(d_delta, d_t_list, d_r_list, d_M_list, self.threadsperblock[3], self.catch_error)
        if not success:
            return {}, False
        self.threadsperblock[3] = threadsperblock

        d_Mtilde = cuda.device_array((n_list.shape[1],2,2,2), dtype=complex)
        threadsperblock, success = make_Mtilde(d_M_list, d_Mtilde, self.threadsperblock[4], self.catch_error)
        if not success:
            return {}, False
        self.threadsperblock[4] = threadsperblock

        # Net complex transmission and reflection amplitudes
        r, t, R, T, threadsperblock, success = get_r_t_R_T(d_Mtilde, d_n_list, d_cos_th_list, self.threadsperblock[5], self.catch_error)
        if not success:
            return {}, False
        self.threadsperblock[5] = threadsperblock

        # vw_list[n] = [v_n, w_n]. v_0 and w_0 are undefined because the 0th medium
        # has no left interface.
        vw_list = np.zeros((num_layers, num_wl, 2, 2), dtype=complex) # last index is polarization
        power_entering = np.zeros((1,2))
        if detailed:
            vw = np.zeros((num_wl, 2, 2, 2), dtype=complex)
            vw[:, 0, 0, 0] = t[:,0]
            vw[:, 0, 1, 0] = t[:,0]
            vw[:, 0, 0, 1] = t[:,1]
            vw[:, 0, 1, 1] = t[:,1]
            vw_list[-1] = vw[:, 0, :, :]
            vw, vw_list, threadsperblock, success = build_vw(d_M_list,vw,vw_list, self.threadsperblock[6], self.catch_error)
            if not success:
                return {}, False
            self.threadsperblock[6] = threadsperblock

            # Net transmitted and reflected power, as a proportion of the incoming light
            # power.
            power_entering = power_entering_from_r(
                r, n_list[0], np.cos(th_0))

        # return both s [0] and p [1] results
        return {'r': r[:,0], 't': t[:,0], 'R': R[:,0], 'T': T[:,0], 'power_entering': power_entering[:,0],
                'vw_list': vw_list[:,:,:,0], 'kz_list': kz_list, 'th_list': th_list,
                'pol': 's', 'n_list': n_list, 'd_list': d_list, 'th_0': th_0,
                'lam_vac': lam_vac}, {'r': r[:,1], 't': t[:,1], 'R': R[:,1], 'T': T[:,1], 'power_entering': power_entering[:,1],
                'vw_list': vw_list[:,:,:,1], 'kz_list': kz_list, 'th_list': th_list,
                'pol': 'p', 'n_list': n_list, 'd_list': d_list, 'th_0': th_0,
                'lam_vac': lam_vac}, True
    
