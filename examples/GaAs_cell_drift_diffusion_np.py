import numpy as np
import matplotlib.pyplot as plt
from solcore import material, si
from solcore.structure import Junction, Layer
from solcore.solar_cell_solver import solar_cell_solver, SolarCell
from solcore.state import State
from solcore.light_source import LightSource

# set up structure similar to https://doi.org/10.1109/16.46385 (but with doping flipped)
# Compare performance of the two PDD solvers.

# define user options
options = State()
options.wavelength = np.linspace(280, 950, 70)*1e-9
options.voltages = np.linspace(-1.2, 0, 60)
options.internal_voltages = np.linspace(-1.2, 0.0, 60)
options.light_iv = True
options.mpp = True
options.optics_method = 'TMM'
options.light_source = LightSource(source_type="standard",
                           version="AM1.5g", x=options.wavelength, output_units="photon_flux_per_m")
options.no_back_reflection = False
options.recalculate_absorption = True


# define materials

MgF2 = material("MgF2")()
ZnS = material("ZnScub")()
Ag = material("Ag")()

ARC_layers = [Layer(si("100nm"), material=MgF2),
              Layer(si("50nm"), material=ZnS)]

AlGaAs_window = material('AlGaAs')(Al=0.8, Nd=si('4e18cm-3'),
                                   relative_permittivity=11,
                                   electron_minority_lifetime=1e-9,
                                   hole_minority_lifetime=1e-9,
                                   )
GaAs_emitter = material('GaAs')(Nd=si('4e18cm-3'),
                                electron_minority_lifetime=5e-9,
                                hole_minority_lifetime=5e-9)
GaAs_base = material('GaAs')(Na=si('2e17cm-3'),
                             electron_minority_lifetime=5e-9,
                             hole_minority_lifetime=5e-9)
GaAs_bsf = material('GaAs')(Na=si('2e18cm-3'),
                            electron_minority_lifetime=5e-9,
                            hole_minority_lifetime=5e-9,
                            )
GaAs_substrate = material('GaAs')()

junction_layers = [
        Layer(si('30nm'), AlGaAs_window, role='window'),
        Layer(si('100nm'), GaAs_emitter, role='emitter'),
        Layer(si('3000nm'), GaAs_base, role='base'),
        Layer(si('100nm'), GaAs_bsf, role='bsf'),
        ]

solar_cell = SolarCell(
    ARC_layers + [Junction(junction_layers, kind='PDD', sn=1e6, sp=1e6,
                           R_shunt=1,
                           ),
                  Layer(si('100um'), GaAs_substrate)
                  ],
    substrate=Ag,
    R_series=0.0001,
)

solar_cell_sesame = SolarCell(
    ARC_layers + [Junction(junction_layers, kind='sesame_PDD', sn=1e6, sp=1e6,
                           R_shunt=1,
                           ),
                  Layer(si('100um'), GaAs_substrate)
                  ],
    substrate=Ag,
    R_series=0.0001,
)

solar_cell_optics = SolarCell(
    ARC_layers + junction_layers + [Layer(si('100um'), GaAs_substrate)], substrate=Ag,
)

solar_cell_solver(solar_cell, 'iv', options)
solar_cell_solver(solar_cell, 'qe', options)

solar_cell_solver(solar_cell_sesame, 'iv', options)
solar_cell_solver(solar_cell_sesame, 'qe', options)

solar_cell_solver(solar_cell_optics, 'optics', options)

absorption_per_layer = np.array([layer.layer_absorption for layer in solar_cell_optics])

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

ax1.plot(solar_cell.iv['IV'][0], solar_cell.iv['IV'][1]/10, '--k', label='Overall IV')
ax1.plot(options.voltages, -solar_cell[2].iv(options.voltages)/10, '-r', label='junction IV')

ax1.legend()
ax1.set_title("Fortran PDD solver")

ax2.stackplot(options.wavelength*1e9, 100*absorption_per_layer[::-1], alpha=0.4)
ax2.plot(options.wavelength*1e9, 100*solar_cell[2].eqe(options.wavelength), '-k')

ax2.legend(['substrate', 'GaAs BSF', 'GaAs base', 'GaAs emitter', 'AlGaAs window', 'ZnS', 'MgF$_2$'])

ax1.set_xlabel('Voltage (V)')
ax1.set_ylabel('Current density (mA/cm$^2$)')
ax1.set_ylim(-35, 10)

ax2.set_xlabel('Wavelength (nm)')
ax2.set_ylabel('EQE / Absorption (%)')

plt.tight_layout()
plt.show()

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

ax1.plot(solar_cell_sesame.iv['IV'][0], solar_cell_sesame.iv['IV'][1]/10, '--k', label='Overall IV')
ax1.plot(options.voltages, -solar_cell_sesame[2].iv(options.voltages)/10, '-r', label='junction IV')

ax1.legend()
ax1.set_title("Sesame PDD solver")

ax2.stackplot(options.wavelength*1e9, 100*absorption_per_layer[::-1], alpha=0.4)
ax2.plot(options.wavelength*1e9, 100*solar_cell_sesame[2].eqe(options.wavelength), '-k')

ax2.legend(['substrate', 'GaAs BSF', 'GaAs base', 'GaAs emitter', 'AlGaAs window', 'ZnS', 'MgF$_2$'])

ax1.set_xlabel('Voltage (V)')
ax1.set_ylabel('Current density (mA/cm$^2$)')
ax1.set_ylim(-35, 10)

ax2.set_xlabel('Wavelength (nm)')
ax2.set_ylabel('EQE / Absorption (%)')

plt.tight_layout()
plt.show()


# show that generation rate in optical solver, integrated with the spectrum, and used/saved by
# Sesame are the same:

z_pos = np.linspace(0, 3230, 4000)*1e-9

gen_per_wl = solar_cell_sesame[2].absorbed(z_pos)

gen_total = np.trapezoid(gen_per_wl * options.light_source.spectrum(options.wavelength, output_units="photon_flux_per_m")[1][None, :], options.wavelength)

plt.figure()
plt.plot(solar_cell_sesame[2].mesh*1e6, solar_cell_sesame[2].pdd_output.G, '-k', label='Generation')
plt.plot(z_pos*1e6, gen_total, '--r')
plt.xlabel('Position (um)')
plt.ylabel('Generation rate (m$^{-3}$s$^{-1}$)')
plt.show()