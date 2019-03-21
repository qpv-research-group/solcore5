import numpy as np
import matplotlib.pyplot as plt

from solcore import material, si

from solcore.solar_cell import SolarCell, Layer, Junction
from solcore.solar_cell_solver import solar_cell_solver, prepare_solar_cell
from solcore.state import State

GaInP = material('GaInP')(In=0.5)
GaAs = material('GaAs')()
Ge = material('Ge')()

optical_struct = SolarCell([Layer(material=GaInP, width=si('5000nm')),
                            Junction([Layer(material=GaAs, width=si('200nm')),
                            Layer(material=GaAs, width=si('5um'))], kind = 'DA'),
                            Layer(material=Ge, width=si('50um'))
                            ])

wl = np.linspace(250, 1700, 100)*1e-9


options = State()
options.position = None
prepare_solar_cell(optical_struct, options)
position = np.arange(0, optical_struct.width, 1e-10)
options.position = position
options.wavelength = wl
options.optics_method = 'TMM'
options.no_back_reflexion = False
options.pol = 'p'
options.BL_correction = True
options.coherency_list = ['c', 'c', 'c', 'c']
options.theta = 30
solar_cell_solver(optical_struct, 'optics', options)

plt.figure(1)
plt.plot(wl*1e9, optical_struct[0].layer_absorption)
plt.plot(wl*1e9, optical_struct[1].layer_absorption)
plt.plot(wl*1e9, optical_struct[2].layer_absorption)
plt.plot(wl*1e9, optical_struct.reflected, '--')
plt.plot(wl*1e9, optical_struct.transmitted, '--')
plt.plot(wl*1e9, optical_struct[0].layer_absorption+optical_struct[1].layer_absorption +
         optical_struct[2].layer_absorption + optical_struct.reflected + optical_struct.transmitted)
plt.legend(['GaInP', 'GaAs', 'Ge', 'R', 'T'])
plt.show()

#plt.figure(1)
#plt.plot(np.linspace(0, 3000, 100),
#         optical_struct[1].absorbed(np.linspace(0, 3000, 100)*1e-9)[:,0])
#plt.show()