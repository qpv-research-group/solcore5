import numpy as np
import matplotlib.pyplot as plt
from solcore import material, si
from solcore.structure import Junction, Layer
from solcore.solar_cell_solver import solar_cell_solver, SolarCell
from solcore.state import State
from solcore.light_source import LightSource

# set up structure similar to https://doi.org/10.1109/16.46385

# define user options
options = State()
options.wavelength = np.linspace(280, 950, 100)*1e-9
options.voltages = np.linspace(0, 1.3, 20)
options.light_iv = True
options.mpp = True
options.optics_method = 'TMM'
options.recalculate_absorption = True
options.light_source = LightSource(source_type="standard",
                           version="AM1.5g", x=options.wavelength, output_units="photon_flux_per_m")


# define materials

MgF2 = material("MgF2")()
ZnS = material("ZnScub")()

ARC_layers = [Layer(si("100nm"), material=MgF2),
              Layer(si("50nm"), material=ZnS)]

AlGaAs_window = material('AlGaAs')(Al=0.8, Na=si('4e18cm-3'),
                                   relative_permittivity=11,
                                   electron_minority_lifetime=1e-6,
                                   hole_minority_lifetime=1e-6
                                   )
GaAs_emitter = material('GaAs')(Na=si('4e18cm-3'),
                                electron_minority_lifetime=1e-6,
                                hole_minority_lifetime=1e-6)
GaAs_base = material('GaAs')(Nd=si('2e17cm-3'),
                             electron_minority_lifetime=1e-6,
                             hole_minority_lifetime=1e-6)
GaAs_bsf = material('GaAs')(Nd=si('2e18cm-3'),
                            electron_minority_lifetime=1e-6,
                            hole_minority_lifetime=1e-6
                            )
GaAs_substrate = material('GaAs')()

junction_layers = [
        Layer(si('30nm'), AlGaAs_window, role='window'),
        Layer(si('100nm'), GaAs_emitter, role='emitter'),
        Layer(si('3000nm'), GaAs_base, role='base'),
        Layer(si('100nm'), GaAs_bsf, role='bsf'),
        ]

# define solar cell
mesh = np.linspace(0, 3230, 4000)*1e-9

solar_cell = SolarCell(
    ARC_layers + [Junction(junction_layers, kind='PDD', sn=1e6, sp=1e6,)],
        substrate=GaAs_substrate)

solar_cell_sesame = SolarCell(
    ARC_layers + [Junction(junction_layers, kind='sesame_PDD', sn=1e6, sp=1e6,
                           mesh=mesh)],
        substrate=GaAs_substrate)

solar_cell_optics = SolarCell(
    ARC_layers + junction_layers, substrate=GaAs_substrate
)

solar_cell_solver(solar_cell_optics, 'optics', options)

solar_cell_solver(solar_cell, 'iv', options)
solar_cell_solver(solar_cell, 'qe', options)

solar_cell_solver(solar_cell_sesame, 'iv', options)
solar_cell_solver(solar_cell_sesame, 'qe', options)

absorption_per_layer = np.array([layer.layer_absorption for layer in solar_cell_optics])

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

ax1.plot(solar_cell.iv['IV'][0], solar_cell.iv['IV'][1]/10, label='IV')

ax2.stackplot(options.wavelength*1e9, 100*absorption_per_layer[::-1], alpha=0.4)
ax2.plot(options.wavelength*1e9, 100*solar_cell[2].eqe(options.wavelength), '-k')

ax2.legend(['GaAs BSF', 'GaAs base', 'GaAs emitter', 'AlGaAs window', 'ZnS', 'MgF$_2$'])

ax1.set_xlabel('Voltage (V)')
ax1.set_ylabel('Current density (mA/cm$^2$)')
ax1.set_ylim(-10, 35)

ax2.set_xlabel('Wavelength (nm)')
ax2.set_ylabel('EQE (%)')

plt.tight_layout()
plt.show()

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

ax1.plot(solar_cell_sesame.iv['IV'][0], solar_cell_sesame.iv['IV'][1]/10, label='IV')

ax2.stackplot(options.wavelength*1e9, 100*absorption_per_layer[::-1], alpha=0.4)
ax2.plot(options.wavelength*1e9, 100*solar_cell_sesame[2].eqe(options.wavelength), '-k')

ax2.legend(['GaAs BSF', 'GaAs base', 'GaAs emitter', 'AlGaAs window', 'ZnS', 'MgF$_2$'])

ax1.set_xlabel('Voltage (V)')
ax1.set_ylabel('Current density (mA/cm$^2$)')
ax1.set_ylim(-10, 35)

ax2.set_xlabel('Wavelength (nm)')
ax2.set_ylabel('EQE / Absorption (%)')

plt.tight_layout()
plt.show()

# from solcore import material, si
# from solcore.structure import Junction, Layer
# from solcore.solar_cell_solver import solar_cell_solver, SolarCell
# from solcore.analytic_solar_cells import iv_depletion
# from solcore.sesame_drift_diffusion import iv_sesame
# from solcore.state import State
# from solcore.optics import solve_tmm
# from solcore.light_source import LightSource
# import numpy as np
#
# options = State()
# options.wavelength = np.linspace(300, 950, 100) * 1e-9
# options.voltages = np.linspace(-1.2, 1.2, 20)
# options.internal_voltages = np.linspace(-1.2, 1.2, 20)
# options.T = 300
# options.light_iv = True
# options.light_source = LightSource(source_type='standard', version='AM1.5g', x=options.wavelength,
#                                    output_units='photon_flux_per_m')
# options.da_mode = 'green'
# options.optics_method = 'TMM'
#
# mesh = np.linspace(0, 3100, 2000)*1e-9
#
# GaAs_n = material('GaAs')(T=300, Na=si('4e18cm-3'), hole_minority_lifetime=1e-6, electron_minority_lifetime=1e-6,
#                           electron_mobility=1, hole_mobility=1)
# GaAs_p = material('GaAs')(T=300, Nd=si('2e17cm-3'), hole_minority_lifetime=1e-6, electron_minority_lifetime=1e-6,
#                           electron_mobility=1, hole_mobility=1)
#
# pn_junction = Junction([Layer(si('100nm'), GaAs_n, role='emitter'), Layer(si('3000nm'), GaAs_p, role='base')],
#                        kind='sesame_PDD',
#                        # mesh=mesh
#                        )
#
# solar_cell_solver(SolarCell([pn_junction]), 'optics', options)
#
# plt.figure()
# iv_sesame(pn_junction, options)
# plt.plot(options.voltages, pn_junction.iv(options.voltages), label='Sesame')
# plt.ylim(-250, 250)
# plt.show()

