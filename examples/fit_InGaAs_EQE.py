import numpy as np
import matplotlib.pyplot as plt

from solcore import material, si
from solcore.solar_cell import SolarCell
from solcore.structure import Layer, Junction
from solcore.state import State
from solcore.solar_cell_solver import solar_cell_solver

InGaAs = material('InGaAs')
InP = material('InP')()

wl = np.linspace(300, 1800, 1000) * 1e-9

emitter = InGaAs(In=0.53, Nd=si('4e18cm-3'), hole_diffusion_length=si('800nm'))
base = InGaAs(In=0.53, Na=si('2e17cm-3'), electron_diffusion_length=si('5000nm'))

solar_cell = SolarCell([
    Layer(si('30nm'), InP),
    Junction([
        Layer(si('440nm'), emitter, role='emitter'),
        Layer(si('3500nm'), base, role='base')], kind='DA',
    sn=100),
    ],
    substrate=InP,
    shading=0.02,
)

options = State()
options.wavelength = wl
options.light_iv = True
options.mpp = True
options.optics_method = 'TMM'
options.no_back_reflection = False
options.voltages = np.linspace(0, 0.35, 40)


solar_cell_solver(solar_cell, 'qe', options)
solar_cell_solver(solar_cell, 'iv', options)

experimental_data = np.loadtxt('data/InGaAs_EQE.txt')

plt.figure()
plt.plot(experimental_data[:, 0], experimental_data[:, 1], label='Experimental')
plt.plot(wl*1e9, 100*solar_cell[1].layer_absorption, '--', label='Absorbed in junction')
plt.plot(wl*1e9, 100*solar_cell[0].layer_absorption + 100*solar_cell[1].layer_absorption, '--',
         label='Absorbed in junction + window')
plt.plot(wl * 1e9, 100*solar_cell[1].eqe(wl), label='EQE')
plt.plot(wl * 1e9, 100*solar_cell[1].eqe_emitter(wl), label='EQE')
plt.plot(wl * 1e9, 100*solar_cell[1].eqe_base(wl), label='EQE')
plt.legend()
plt.show()

plt.figure()
plt.plot(options.voltages, -solar_cell[1].iv(options.voltages))
plt.ylim(-400, 400)
plt.show()

# fit diffusion lengths,