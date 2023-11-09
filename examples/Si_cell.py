from solcore.solar_cell import SolarCell, Junction, Layer
from solcore.state import State
from solcore.solar_cell_solver import solar_cell_solver
from solcore.sesame_drift_diffusion.solve_pdd import process_structure
import numpy as np
import os
from solcore.light_source import LightSource
from solcore.absorption_calculator import search_db
from scipy.special import erfc
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from cycler import cycler
from solcore.constants import q
from scipy.interpolate import interp1d, RectBivariateSpline

from solcore import material, si

d_bulk = 130e-6

Air = material("Air")()
MgF2_pageid = str(search_db(os.path.join("MgF2", "Rodriguez-de Marcos"))[0][0])
MgF2 = material(MgF2_pageid, nk_db=True)()
Ag = material("Ag_Jiang")()

Si_pn = material("Si")(electron_mobility=si("1e4cm2"), hole_mobility=si("1e3cm2"),
                       electron_minority_lifetime=0.01, hole_minority_lifetime=0.01)

wavelengths = np.linspace(280, 1200, 100)*1e-9

TCO = material('ITO2')()


front_materials = [Layer(80e-9, MgF2), Layer(55e-9, TCO)]

back_materials = [Layer(55e-9, TCO),
                  Layer(120e-9, MgF2)]

options = State()
options.wavelength = wavelengths
options.optics_method = 'TMM'
options.light_iv = True
options.T = 298

options.light_source = LightSource(source_type="standard",
                           version="AM1.5g", x=options.wavelength, output_units="photon_flux_per_m")


options.voltages = np.linspace(0, 0.8, 30)
options.internal_voltages = options.voltages
options.mpp = True
options.recalculate_absorption = True

nD = si("1e20cm-3")
nA = si("1e19cm-3")
bulk_doping = si("5e15cm-3") # n type bulk

# rear junction (n-type)
def doping_profile_func(x):

    L = d_bulk

    doping_profile = - nA * erfc(x/150e-9) # characteristic depth of 150 nm

    doping_profile_rear = nD * erfc((L - x)/200e-9) # characteristic depth of 200 nm

    return doping_profile + doping_profile_rear + bulk_doping

L_cm = d_bulk*100 # length of the system in the x-direction [cm]

surface_recomb = dict(sn_front=si('5e3 cm s-1'),  # important
                          sp_front=si('1e4 cm s-1'),
                          sn_rear=si('1e4 cm s-1'),
                          sp_rear=si('1e4 cm s-1'))  # important

Si_junction = [Junction([Layer(d_bulk, Si_pn)],
                        doping_profile=doping_profile_func, kind='sesame_PDD',
                        **surface_recomb,
                        )
               ]


Si_cell = SolarCell(front_materials +
                     Si_junction +
                     back_materials,
                    shading=0.02,
                    substrate=Air,
                    )


solar_cell_solver(Si_cell, 'iv', options)

solar_cell_solver(Si_cell, 'qe', options)

shading = 0.02

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10, 3.5))
# ax.stackplot(wavelengths * 1e9, 100 * (1 - shading) * np.array(result_stack), linewidth=0.5)
# ax.fill_between(wavelengths * 1e9, 100 * (1 - shading), 100, color='white')
ax.plot(wavelengths * 1e9, 100 * Si_cell(0).eqe(wavelengths), '-r', linewidth=2)
ax.legend(['Si absorption', 'Front stack absorption', 'Rear stack aborption',
           'Front reflection', 'Front escape', 'Electrode shading', 'Simulated EQE'],
          loc='lower center')
ax.set_xlim(280, 1200)
ax.set_ylim(0, 100)
ax.set_xlabel("Wavelength (nm)")
ax.set_ylabel("R / A / EQE (%)")
ax.set_title('a) EQE and cell optics', loc='left')
# plt.show()

jsc = Si_cell.iv.Isc / 10

ax2.plot(options.voltages, -(1 - shading) * Si_cell(0).iv(options.voltages) / 10, '-', label='IV',
         linewidth=2, color='k')
# plt.plot(options.voltages, solar_cell_old[1].iv(options.voltages), 'o-', label='IV')
# plt.xlim(-2, 1.8)

ax2.set_ylim(0, 1.03 * jsc)
ax2.set_xlim(np.min(options.voltages), np.max(options.voltages))
ax2.set_xlabel('Voltage (V)')
ax2.set_ylabel('Current density (mA/cm$^2$)')
ax2.set_title('b) IV characteristics and power output', loc='left')

ax3 = ax2.twinx()
ax3.plot(options.voltages, Si_cell.iv['IV'][0] * Si_cell.iv['IV'][1],
         '-r', label='Power', linewidth=2)
ax3.set_ylabel('Power density (W m$^{-2}$)')
ax3.set_ylim(0, 1.03 * jsc * 10)

# ax3.spines['right'].set_color(pal[2])
# ax3.yaxis.label.set_color(pal[2])
# ax3.tick_params(axis='y', colors=pal[2])

ax2.set_axisbelow(True)
ax3.set_axisbelow(True)

ax2.text(0.02, 0.9 * jsc, r'$J_{SC}$', zorder=5)
ax2.text(0.02, 0.8 * jsc, r'$V_{OC}$')
ax2.text(0.02, 0.7 * jsc, 'FF')
ax2.text(0.02, 0.6 * jsc, r'$\eta$')
ax2.text(0.02, 0.5 * jsc, r'$J_{MPP}$')
ax2.text(0.02, 0.4 * jsc, r'$V_{MPP}$')

ax2.text(0.1, 0.9 * jsc, r'= {:.2f} mA/cm$^2$'.format(jsc))
ax2.text(0.1, 0.8 * jsc, r'= {:.3f} V'.format(Si_cell.iv.Voc))
ax2.text(0.1, 0.7 * jsc, '= {:.2f} %'.format(Si_cell.iv.FF * 100))
ax2.text(0.1, 0.6 * jsc, r'= {:.2f} %'.format((1 - shading) * Si_cell.iv.Pmpp / 10))
ax2.text(0.1, 0.5 * jsc, r'= {:.2f} mA/cm$^2$'.format((1 - shading) * Si_cell.iv.Impp / 10))
ax2.text(0.1, 0.4 * jsc, r'= {:.3f} V'.format(Si_cell.iv.Vmpp))
ax2.grid(which='major', alpha=0.35)

ax3.grid(False)
plt.tight_layout()

plt.show()
