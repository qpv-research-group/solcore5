from solcore import material, si
from solcore.solar_cell import SolarCell, Junction, Layer
from solcore.state import State
from solcore.solar_cell_solver import solar_cell_solver
from solcore.light_source import LightSource
import numpy as np
import matplotlib.pyplot as plt

options = State()
# options.wavelength = np.linspace(280, 700, 20)*1e-9
# options.optics_method = 'TMM'
# options.light_iv = True
# options.light_source = LightSource(source_type="standard",
#                            version="AM1.5g", x=options.wavelength, output_units="photon_flux_per_m")
# options.voltages = np.linspace(0, 2, 40)
# # options.minimum_spacing = 1e-09
# # options.maximum_spacing = 1e-8
#
# T = 293
#
# add_args = {'relative_permittivity': 10, 'electron_minority_lifetime': 5e-6,
#             'hole_minority_lifetime': 5e-6,
#             'electron_auger_recombination': 1e-45,
#             'hole_auger_recombination': 1e-45}
#
# ARC = material('Si3N4')()
# window = material('AlGaAs')(T=T, Na=5e24, Al=0.8, **add_args)
# p_AlGaAs = material('AlGaAs')(T=T, Na=1e24, Al=0.4, **add_args)
# n_AlGaAs = material('AlGaAs')(T=T, Nd=8e22, Al=0.4, **add_args)
# bsf = material('AlGaAs')(T=T, Nd=2e24, Al=0.6, **add_args)
#
# junction = Junction([Layer(width=si('30nm'), material=window, role="Window"),
#                    Layer(width=si('1000nm'), material=p_AlGaAs, role="Emitter"),
#                    Layer(width=si('150nm'), material=n_AlGaAs, role="Base"),
#                    Layer(width=si('200nm'), material=bsf, role="BSF")], sn=1e6, sp=1e6, T=T, kind='sesame_PDD')
#
# widths = [layer.width for layer in junction]
#
# # x = np.concatenate((np.linspace(0,3e-7, 400, endpoint=False),
# #                     np.linspace(3e-7, np.sum(widths) - 100e-9, 100, endpoint=False),
# #                     np.linspace(np.sum(widths) - 100e-9, np.sum(widths), 200)))
# # junction.mesh = x
#
# solar_cell = SolarCell(
#     [Layer(60e-0, ARC), junction]
# )
#
# solar_cell_solver(solar_cell, 'iv', options)
# solar_cell_solver(solar_cell, 'qe', options)
#
# # # qe_sesame(solar_cell[1], options)
# #
# # iv_sesame(solar_cell[1], options)
#
# junction_old = Junction([Layer(width=si('30nm'), material=window, role="Window"),
#                    Layer(width=si('150nm'), material=p_AlGaAs, role="Emitter"),
#                    Layer(width=si('1000nm'), material=n_AlGaAs, role="Base"),
#                    Layer(width=si('200nm'), material=bsf, role="BSF")], sn=1e6, sp=1e6, T=T,
#                         kind='PDD')
#
#
# solar_cell_old = SolarCell(
#     [Layer(60e-0, ARC), junction_old]
# )
#
# # solar_cell_solver(solar_cell_old, 'iv', options)
# #
# plt.figure()
# plt.plot(options.voltages, solar_cell[1].iv(options.voltages), 'o-', label='IV')
# # plt.plot(options.voltages, solar_cell_old[1].iv(options.voltages), 'o-', label='IV')
# # plt.xlim(-2, 1.8)
# plt.ylim(-200, 200)
# plt.show()
#
# plt.figure()
# plt.plot(options.wavelength*1e9, solar_cell.absorbed)
# # plt.plot(options.wavelength*1e9, solar_cell_old.absorbed)
# plt.show()
#
# plt.figure()
# plt.plot(options.wavelength*1e9, solar_cell[1].eqe(options.wavelength))
# plt.plot(options.wavelength*1e9, solar_cell[1].iqe(options.wavelength))
# plt.show()
#
# solar_cell[1].absorbed(solar_cell[1].mesh)
#
# plt.figure()
# plt.plot(junction.mesh, 'x')
# plt.show()
#
#

#####

import os
from solcore.constants import q
from solcore.light_source import LightSource
from solcore.interpolate import interp1d
from solcore.absorption_calculator import search_db
from scipy.special import erfc
import matplotlib.pyplot as plt
from rayflare.textures import regular_pyramids

from rayflare.transfer_matrix_method import tmm_structure

from rayflare.textures import regular_pyramids, rough_pyramids
from solcore import material, si
from rayflare.options import default_options
from solcore.structure import Layer
from rayflare.ray_tracing import rt_structure
from rayflare.utilities import make_absorption_function
from os.path import exists

# Note Sesame assumes all quantities of length are input in units of cm. Other assumed
# input units are: time in s, energy in eV.

# LONGi paper:
# To further increase JSC, a 150-nm-thick MgF2 film was evaporated on the front TCO layer
# as a second antireflective coating. For 26.81% cell, an additional 120-nm-thick MgF2/150-nm-thick
# Ag stack was evaporated on the rear TCO layer, which means this cell is a monofacial solar cell.

options = State()
options.wavelength = np.linspace(280, 1200, 20)*1e-9
options.optics_method = 'TMM'
options.light_iv = True
options.light_source = LightSource(source_type="standard",
                           version="AM1.5g", x=options.wavelength, output_units="photon_flux_per_m")
options.voltages = np.linspace(0, 0.8, 20)
options.mpp = True

# ray-tracing for optics first
d_bulk = 130e-6

# ITO_front = material("ITO_front_transparent")()
ITO_front = material("ITO_front_transparent")()
ITO_back = material("ITO_back")()
Ag = material("Ag_Jiang")()
aSi_i = material("aSi_i")()
aSi_p = material("aSi_p")()
aSi_n = material("aSi_n")()
Air = material("Air")()
MgF2_pageid = str(search_db(os.path.join("MgF2", "Rodriguez-de Marcos"))[0][0])
MgF2 = material(MgF2_pageid, nk_db=True)()
Ag = material("Ag_Jiang")()

Si_pn = material("Si")(Nc=3e25, Nv=1e25, electron_mobility=si("1e4cm2"), hole_mobility=si("1e3cm2"),
                       electron_minority_lifetime=0.05, hole_minority_lifetime=0.05,
                       radiative_recombination=si("1.89e-15cm3"),
                       electron_auger_recombination=si("3e-31cm6"), hole_auger_recombination=si("1e-31cm6"))

front_materials = [Layer(110e-9, MgF2), Layer(40 * 1e-9, ITO_front),
                   Layer(1e-9, aSi_i), Layer(1e-9, aSi_n)]
back_materials = [Layer(2e-9, aSi_i), Layer(2e-9, aSi_p), Layer(100e-9, ITO_back),
                  Layer(120e-9, MgF2), Layer(150e-9, Ag)]

front_width = np.sum([layer.width for layer in front_materials])
back_width = np.sum([layer.width for layer in back_materials])

L = 130e-6

# rear junction

nD = si("1e20cm-3")
nA = si("1e19cm-3")
bulk_doping = si("1e15cm-3") # n type bulk
def doping_profile_func(x):

    L = d_bulk

    doping_profile = - nA * erfc(x/150e-9) # characteristic depth of 150 nm

    doping_profile_rear = nD * erfc((L - x)/200e-9) # characteristic depth of 200 nm

    return doping_profile + doping_profile_rear + bulk_doping

L_cm = d_bulk*100 # length of the system in the x-direction [cm]

# Mesh
x = np.concatenate((np.linspace(0,2e-4, 2000, endpoint=False),
                    np.linspace(2e-4, L_cm - 2e-4, 2000, endpoint=False),
                    np.linspace(L_cm - 2e-4, L_cm, 2000)))

options.position = np.concatenate((np.linspace(0, front_width, 100), front_width + x/100,
                                  np.linspace(front_width + d_bulk,
                                              front_width + d_bulk + back_width, 100)))

Si_junction = [Junction([Layer(d_bulk, Si_pn)],
                        doping_profile=doping_profile_func, kind='sesame_PDD',
                        mesh=x/100)]

Si_cell = SolarCell(front_materials +
                     Si_junction +
                     back_materials)

Si_ind = len(front_materials)

solar_cell_solver(Si_cell, 'qe', options)

plt.figure()
plt.plot(options.wavelength*1e9, Si_cell[Si_ind].layer_absorption, '--', label='Si absorption')
plt.plot(options.wavelength*1e9, Si_cell[Si_ind].eqe(options.wavelength), '-k', label='Si EQE')
plt.plot(options.wavelength*1e9, Si_cell[Si_ind].iqe(options.wavelength), '-r', label='Si IQE')
plt.legend()
plt.show()

solar_cell_solver(Si_cell, 'iv', options)

jsc = Si_cell.iv.Isc/10

plt.figure()
plt.plot(options.voltages, -Si_cell[Si_ind].iv(options.voltages)/10, 'o-', label='IV')
# plt.plot(options.voltages, solar_cell_old[1].iv(options.voltages), 'o-', label='IV')
# plt.xlim(-2, 1.8)
plt.text(0.02, 0.9*jsc, r'Jsc = {:.2f} mA/cm$^2$'.format(jsc))
plt.text(0.02, 0.8*jsc, 'Voc = {:.3f} V'.format(Si_cell.iv.Voc))
plt.text(0.02, 0.7*jsc, 'FF = {:.3f} %'.format(Si_cell.iv.FF*100))
plt.text(0.02, 0.6*jsc, 'Eff = {:.3f} %'.format(Si_cell.iv.Eta*100))
plt.text(0.02, 0.5*jsc, r'Jmpp = {:.3f} mA/cm$^2$'.format(Si_cell.iv.Impp/10))
plt.text(0.02, 0.4*jsc, 'Vmpp = {:.3f} V'.format(Si_cell.iv.Vmpp))
plt.ylim(0, 1.1*jsc)
plt.xlim(np.min(options.voltages), np.max(options.voltages))
plt.xlabel('Voltage (V)')
plt.ylabel('Current density (mA/cm$^2$)')
plt.show()

# plt.figure()
# plt.plot(x*1e4, Si_cell[Si_ind].sesame_sys.g, 'o')
# plt.xlim(0, 0.1)
# plt.show()
#
# plt.figure()
# plt.semilogy(Si_junction[0].mesh*1e6, Si_cell[Si_ind].absorbed(Si_junction[0].mesh), 'o', markerfacecolor='none')
# plt.xlim(0, 1)
# plt.ylim(1e-3,1e8)
# plt.show()


x = np.linspace(0, d_bulk, int(1e7))
int_A = np.trapz(Si_cell[Si_ind].absorbed(x), x, axis=0)
