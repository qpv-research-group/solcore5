import numpy as np
import matplotlib.pyplot as plt

from solcore import material, si
from solcore.structure import Layer, Structure
from solcore.absorption_calculator.rigorous_coupled_wave import calculate_rat_rcwa

T = 300
# define materials
Air = material("Air")(T=T)
TiO2 = material("TiO2", sopra=True)(T=T)  # for the nanoparticles
GaAs = material("GaAs")(T=T)

# define a flat layer and another with circular discs with the same thickness
Flat = Layer(si('50nm'), TiO2)
NP_layer = Layer(si('50nm'), Air, geometry=[{'type': 'circle', 'mat': TiO2, 'center': (200, 200), 'radius': 50}])

flat_struct = Structure([Flat])
np_struct = Structure([NP_layer])

# And the wavelength in which the solve the problem
wl = np.linspace(300, 1000, 150)

rat_flat = calculate_rat_rcwa(flat_struct, size=((400, 0), (0, 400)), orders=10, wavelength=wl,
                              substrate=GaAs, incidence=Air)
rat_np = calculate_rat_rcwa(np_struct, size=((400, 0), (0,  400)), orders=10, wavelength=wl,
                            substrate=GaAs, incidence=Air)

plt.plot(wl, rat_flat["R"], '-k', label="Flat (R)")
plt.plot(wl, rat_np["R"], '-b', label="Nanoparticles (R)")
plt.plot(wl, rat_flat["T"], '--k', label="Flat (T)")
plt.plot(wl, rat_np["T"], '--b', label="Nanoparticles (T)")
plt.legend()
plt.show()