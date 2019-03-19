import numpy as np
import matplotlib.pyplot as plt

from solcore import si, material
from solcore.structure import Junction, Layer
from solcore.solar_cell import SolarCell
from solcore.solar_cell_solver import solar_cell_solver, default_options
from solcore.light_source import LightSource

# user options
T = 298
wl = si(np.linspace(400, 1700, 200), 'nm')
light_source = LightSource(source_type='standard', version='AM1.5g', x=wl,
                           output_units='photon_flux_per_m', concentration=1)
opts = default_options
opts.wavelength, opts.no_back_reflexion, opts.size, opts.light_source, opts.T_ambient = \
    wl, False, [400, 400], light_source, T
# The size of the unit cell for the RCWA structure is 400 x 400 nm

# Defining all the materials we need
Air = material('Air')(T=T)
p_GaAs = material('GaAs')(T=T, Na=si('4e18cm-3'))  # for the GaAs cell emitter
n_GaAs = material('GaAs')(T=T, Nd=si('2e17cm-3'))  # for the GaAs cell base
AlAs, GaAs = material('AlAs')(T=T), material('GaAs')(T=T)  # for the DBR
SiO2 = material('SiO2', sopra=True)(T=T)  # for the spacer layer
TiO2 = material('TiO2', sopra=True)(T=T)  # for the nanoparticles
GaInP = material('GaInP')(In=0.5)
Ge = material('Ge', sopra=True)()
# some parameters for the QE solver
for mat in [n_GaAs, p_GaAs]:
    mat.hole_mobility, mat.electron_mobility, mat.permittivity = 3.4e-3, 5e-2, 9
    n_GaAs.hole_diffusion_length, p_GaAs.electron_diffusion_length = si("500nm"), si("5um")

# Define the different parts of the structure we will use. For the GaAs junction, we use the depletion approximation
GaAs_junction = [Junction([Layer(width=si('100nm'), material=p_GaAs, role="emitter"),
                           Layer(width=si('400nm'), material=n_GaAs, role="base")], T=T, kind='DA')]

# this creates 10 repetitions of the AlAs and GaAs layers, to make the DBR structure
DBR = 10 * [Layer(width=si("73nm"), material=AlAs), Layer(width=si("60nm"), material=GaAs)]

# The layer with nanoparticles
NP_layer = [Layer(si('50nm'), Air, geometry=[{'type': 'circle', 'mat': TiO2, 'center': (200, 200),
                                              'radius': 50}])]

spacer = [Layer(width=si('25nm'), material=GaInP)]

bottom = [Layer(width=si('10um'), material = Ge)]

# --------------------------------------------------------------------------
# solar cell with SiO2 coating
coherency_list = ['c', 'c', 'c', 'i']

solar_cell = SolarCell(spacer + GaAs_junction + bottom)

opts.optics_method = 'TMM'
opts.no_back_reflexion = False
opts.pol = 'u'
opts.coherency_list = coherency_list
solar_cell_solver(solar_cell, 'optics', opts)

plt.figure()
plt.plot(wl*1e9, solar_cell[0].layer_absorption)
plt.plot(wl*1e9, solar_cell[1].layer_absorption)
plt.plot(wl*1e9, solar_cell[2].layer_absorption)

plt.plot(wl*1e9, solar_cell.reflected)
plt.plot(wl*1e9, solar_cell.transmitted)
plt.legend(['spacer', 'GaAs', 'Ge', 'R', 'T'])
plt.ylim(0, 1)
plt.show()

plt.figure()
plt.plot(np.linspace(0, 500, 500), solar_cell[1].absorbed(np.linspace(0, 500, 500)*1e-9)[:,4])
plt.show()