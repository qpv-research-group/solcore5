import numpy as np
from solcore import material, si


Si = material('Si')()
Air = material('Air')()

wavelength = si(np.linspace(300, 1000, 2), 'nm')

def Bruggeman_EMA(x, mat1, mat2, wl):
    # x is the fraction of material 1
    e1 = (mat1.n(wl) + 1j*mat1.k(wl))**2
    e2 = (mat2.n(wl) + 1j*mat2.k(wl))**2

    b = e1*(3*x-1) + e2*(2-3*x)

    e_eff_plus = (1/4)*(b + np.sqrt(b**2 + 8*e1*e2))
    #e_eff_minus = (1/4)*(b - np.sqrt(b**2 + 8*e1*e2))

    nk_plus = np.sqrt(e_eff_plus)
    #nk_minus = np.sqrt(e_eff_minus)

    return nk_plus#, nk_minus

nk_eff = Bruggeman_EMA(0.5, Air, Si, wavelength)