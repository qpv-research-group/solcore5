from solcore import material, si
from solcore.solar_cell import SolarCell, Junction, Layer
from solcore.state import State
from solcore.solar_cell_solver import solar_cell_solver
from solcore.light_source import LightSource
import numpy as np
import matplotlib.pyplot as plt

options = State()
options.wavelength = np.linspace(280, 700, 20)*1e-9
options.optics_method = 'TMM'
options.light_iv = True
options.light_source = LightSource(source_type="standard",
                           version="AM1.5g", x=options.wavelength, output_units="photon_flux_per_m")
options.voltages = np.linspace(0, 2, 40)

T = 293

add_args = {'relative_permittivity': 10, 'electron_minority_lifetime': 5e-6,
            'hole_minority_lifetime': 5e-6,
            'electron_auger_recombination': 1e-45,
            'hole_auger_recombination': 1e-45}

ARC = material('Si3N4')()
window = material('AlGaAs')(T=T, Na=5e24, Al=0.8, **add_args)
p_AlGaAs = material('AlGaAs')(T=T, Na=1e24, Al=0.4, **add_args)
n_AlGaAs = material('AlGaAs')(T=T, Nd=8e22, Al=0.4, **add_args)
bsf = material('AlGaAs')(T=T, Nd=2e24, Al=0.6, **add_args)

junction = Junction([Layer(width=si('30nm'), material=window, role="Window"),
                   Layer(width=si('150nm'), material=p_AlGaAs, role="Emitter"),
                   Layer(width=si('1000nm'), material=n_AlGaAs, role="Base"),
                   Layer(width=si('200nm'), material=bsf, role="BSF")], sn=1e6, sp=1e6, T=T, kind='sesame_PDD')

widths = [layer.width for layer in junction]
junction.mesh = np.linspace(0, np.sum(widths), 1000)

solar_cell = SolarCell(
    [Layer(60e-0, ARC), junction]
)

# solar_cell_solver(solar_cell, 'iv', options)
solar_cell_solver(solar_cell, 'qe', options)

# # qe_sesame(solar_cell[1], options)
#
# iv_sesame(solar_cell[1], options)

junction_old = Junction([Layer(width=si('30nm'), material=window, role="Window"),
                   Layer(width=si('150nm'), material=p_AlGaAs, role="Emitter"),
                   Layer(width=si('1000nm'), material=n_AlGaAs, role="Base"),
                   Layer(width=si('200nm'), material=bsf, role="BSF")], sn=1e6, sp=1e6, T=T,
                        kind='PDD')


solar_cell_old = SolarCell(
    [Layer(60e-0, ARC), junction_old]
)

# solar_cell_solver(solar_cell_old, 'iv', options)
#
# plt.figure()
# plt.plot(options.voltages, solar_cell[1].iv(options.voltages), 'o-', label='IV')
# plt.plot(options.voltages, solar_cell_old[1].iv(options.voltages), 'o-', label='IV')
# # plt.xlim(-2, 1.8)
# plt.ylim(-400, 400)
# plt.show()

plt.figure()
plt.plot(options.wavelength*1e9, solar_cell.absorbed)
# plt.plot(options.wavelength*1e9, solar_cell_old.absorbed)
plt.show()

plt.figure()
plt.plot(options.wavelength*1e9, solar_cell[1].eqe(options.wavelength))
plt.plot(options.wavelength*1e9, solar_cell[1].iqe(options.wavelength))
plt.show()