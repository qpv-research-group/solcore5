from solcore import material
from solcore.structure import Layer, Junction, TunnelJunction
from solcore.solar_cell import SolarCell
from solcore.solar_cell_solver import solar_cell_solver
from solcore.light_source import LightSource
import solcore.poisson_drift_diffusion as PDD

import numpy as np
import matplotlib.pyplot as plt

T = 298
wl = np.linspace(350, 1050, 601) * 1e-9

# First, we create the materials of the QW
QWmat = material('InGaAs')(T=T, In=0.2, strained=True)
Bmat = material('GaAsP')(T=T, P=0.1, strained=True)
i_GaAs = material('GaAs')(T=T)

# The QW is 7 nm wide, with GaAs interlayers 2 nm thick at each side and GaAsP barriers 10 nm thick.
# The final device will have 30 of these QWs.
QW = PDD.QWunit([Layer(width=10e-9, material=Bmat, role="barrier"),
                 Layer(width=2e-9, material=i_GaAs, role="well"),
                 Layer(width=7e-9, material=QWmat, role="well"),
                 Layer(width=2e-9, material=i_GaAs, role="well"),
                 Layer(width=10e-9, material=Bmat, role="barrier")], T=T, repeat=30, substrate=i_GaAs)

# We solve the quantum properties of the QW, leaving the default values of all parameters
QW_list = QW.GetEffectiveQW(wavelengths=wl)

# Materials for the BOTTOM junction
window_bottom = material('GaInP')(T=T, Nd=5e24, In=0.49)
n_GaAs = material('GaAs')(T=T, Nd=1e24)
p_GaAs = material('GaAs')(T=T, Na=8e22)
bsf_bottom = material('GaInP')(T=T, Na=5e24, In=0.49)

# If you want to test the code without QWs, to make ti a bit faster, comment the line with QW_list
GaAs_junction = Junction([Layer(width=10e-9, material=window_bottom, role="Window"),
                          Layer(width=150e-9, material=n_GaAs, role="Emitter")] +
                         QW_list +
                         [Layer(width=2000e-9, material=p_GaAs, role="Base"),
                          Layer(width=200e-9, material=bsf_bottom, role="BSF")],
                         sn=1e6, sp=1e6, T=T, kind='PDD')

# Materials for the TOP junction
window_top = material('AlInP')(T=T, Nd=5e23, Al=0.53, electron_mobility=0.01, hole_mobility=7e-4)
n_GaInP = material('GaInP')(T=T, Nd=5e23, In=0.49)
p_GaInP = material('GaInP')(T=T, Na=8e22, In=0.49)
bsf_top = material('AlInP')(T=T, Na=5e23, Al=0.53, electron_mobility=0.01, hole_mobility=7e-4)

# We create the Top junction, leaving outside the window and bsf layers. They work well in the dark, but cause
# convergence issues under illumination.
GaInP_junction = Junction([Layer(width=120e-9, material=n_GaInP, role="Emitter"),
                           Layer(width=800e-9, material=p_GaInP, role="Base")], sn=1e3, sp=1e3, T=T, kind='PDD')

# We create the tunnel junction
tunnel = TunnelJunction([Layer(width=40e-9, material=n_GaInP, role="TJ")],
                        v_peak=0.2, j_peak=7.5e4, v_valley=1, j_valley=4e4, prefactor=5, j01=1e-23, kind='parametric',
                        pn=True)

# And the materials needed for the anti reflecting coating
MgF2 = material('MgF2')()
ZnS = material('ZnScub')()

# Finally, we put everithing together to make a solar cell
my_solar_cell = SolarCell([Layer(width=110e-9, material=MgF2, role="ARC1"),
                           Layer(width=60e-9, material=ZnS, role="ARC2"),
                           Layer(width=30e-9, material=window_top, role="window"),
                           GaInP_junction,
                           Layer(width=100e-9, material=bsf_top, role="BSF"),
                           tunnel,
                           GaAs_junction],
                          T=T, substrate=n_GaAs)

light_source = LightSource(source_type='standard', version='AM1.5g', x=wl,
                           output_units='photon_flux_per_m', concentration=1)

# The definitions are all done, so we just start solving the properties, starting with the QE.

# We calculate the QE curve under illumination
solar_cell_solver(my_solar_cell, 'qe',
                  user_options={'light_source': light_source, 'wavelength': wl, 'optics_method': 'TMM'})

# And now, the IV curves under various concentration levels.
# NOTE: Due to the presence of QWs and the fact we calculate things a 19 different concentrations, this might take a
# while (~4 hours). Remove the QWs as indicated above to test the code much faster.

num_con = 3
con = np.logspace(0, 3, num_con)
vint = np.linspace(-3.5, 4, 600)
V = np.linspace(-3.5, 0, 300)

allI = []
isc = []
voc = []
FF = []
pmpp = []

fig3, axIV = plt.subplots(1, 1, figsize=(6, 4))
for c in con:
    light_source.options['concentration'] = c

    solar_cell_solver(my_solar_cell, 'iv',
                      user_options={'light_source': light_source, 'wavelength': wl, 'optics_method': None,
                                    'light_iv': True, 'mpp': True, 'voltages': V, 'internal_voltages': vint})

    isc.append(my_solar_cell.iv['Isc'])
    voc.append(my_solar_cell.iv['Voc'])
    FF.append(my_solar_cell.iv['FF'])
    pmpp.append(my_solar_cell.iv['Pmpp'])
    allI.append(my_solar_cell.iv['IV'][1])

    # And now, everything is plotting...
    axIV.plot(-V, my_solar_cell.iv['IV'][1] / isc[-1], label=int(c))

axIV.legend(loc='lower left', frameon=False)
axIV.set_ylim(0, 1.1)
axIV.set_xlim(0, 3.5)
axIV.set_xlabel('Voltage (V)')
axIV.set_ylabel('Normalized current (-)')

plt.tight_layout()

fig2, axes = plt.subplots(2, 2, figsize=(11.25, 8))

axes[0, 0].semilogx(con, np.array(pmpp) / con / 10, 'r-o')
axes[0, 0].set_xlabel('Concentration (suns)')
axes[0, 0].set_ylabel('Efficiency (%)')

axes[0, 1].loglog(con, abs(np.array(isc)), 'b-o')
axes[0, 1].set_xlabel('Concentration (suns)')
axes[0, 1].set_ylabel('I$_{SC}$ (Am$^{-2}$)')

axes[1, 0].semilogx(con, abs(np.array(voc)), 'g-o')
axes[1, 0].set_xlabel('Concentration (suns)')
axes[1, 0].set_ylabel('V$_{OC}$ (V)')

axes[1, 1].semilogx(con, abs(np.array(FF)) * 100, 'k-o')
axes[1, 1].set_xlabel('Concentration (suns)')
axes[1, 1].set_ylabel('Fill Factor (%)')

plt.tight_layout()

# We can plot the electron and hole densities in equilibrium and at short circuit
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.25, 4))
for j in my_solar_cell.junction_indices:
    zz = my_solar_cell[j].short_circuit_data.Bandstructure['x'] + my_solar_cell[j].offset
    n = my_solar_cell[j].short_circuit_data.Bandstructure['n']
    p = my_solar_cell[j].short_circuit_data.Bandstructure['p']
    ax1.semilogy(zz * 1e9, n, 'b')
    ax1.semilogy(zz * 1e9, p, 'r')

    zz = my_solar_cell[j].equilibrium_data.Bandstructure['x'] + my_solar_cell[j].offset
    n = my_solar_cell[j].equilibrium_data.Bandstructure['n']
    p = my_solar_cell[j].equilibrium_data.Bandstructure['p']
    ax1.semilogy(zz * 1e9, n, 'b--')
    ax1.semilogy(zz * 1e9, p, 'r--')

ax1.set_xlabel('Position (nm)')
ax1.set_ylabel('Carrier density (m$^{-3}$)')
plt.tight_layout()

# And we plot the QE
labels = ['EQE GaInP', 'EQE GaAs']
colours = ['b', 'r']
for i, j in enumerate(my_solar_cell.junction_indices):
    ax2.plot(wl * 1e9, my_solar_cell[j].eqe(wl), colours[i], label=labels[i])

ax2.plot(wl * 1e9, my_solar_cell.absorbed, 'k', label='Total Absorbed')

ax2.legend(loc='upper right', frameon=False)
ax2.set_xlabel('Wavelength (nm)')
ax2.set_ylabel('EQE')
ax2.set_ylim(0, 1.1)
ax2.set_xlim(350, 1150)

plt.tight_layout()

plt.show()
