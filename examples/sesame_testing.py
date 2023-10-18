import sesame
import numpy as np
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

force_recalculate = True

step_size = 20

shading = 0.02

wls = np.arange(280, 1200, step_size) * 1e-9

wls_ls = np.linspace(250, 4500, 4000)
base_light = LightSource(source_type="standard",
                           version="AM1.5g", x=wls_ls, output_units="photon_flux_per_nm")

base_flux = base_light.spectrum()[1]

# ray-tracing for optics first
d_bulk = 130e-6

# photons at wavelength shorter than the singlet energy S1 are assumed to be absorbed,
# and turned into two photons through singlet fission. Assuming thermally neutral SF, these
# two photons will have half the energy/two times the wavelength of S1.

absorption_data = f'rt_data_{step_size}3.npy'

# ITO_front = material("ITO_front_transparent")()
ITO_front = material("ITO_front_transparent")()
ITO_back = material("ITO_back")()
Ag = material("Ag_Jiang")()
aSi_i = material("aSi_i")()
aSi_p = material("aSi_p")()
aSi_n = material("aSi_n")()
Si = material("Si")()
Air = material("Air")()
SiNx = material("SiNx_Ox")()
MgF2_pageid = str(search_db(os.path.join("MgF2", "Rodriguez-de Marcos"))[0][0])
MgF2 = material(MgF2_pageid, nk_db=True)()
Ag = material("Ag_Jiang")()


plt.figure()
plt.plot(wls*1e9, SiNx.n(wls))
plt.plot(wls*1e9, SiNx.k(wls))
plt.plot(wls*1e9, ITO_front.n(wls), '--')
plt.plot(wls*1e9, ITO_front.k(wls), '--')
plt.plot(wls*1e9, ITO_back.n(wls), '-.')
plt.plot(wls*1e9, ITO_back.k(wls), '-.')
plt.show()
# build structure

options = default_options()
options.wavelengths = wls
options.pol = "u"
options.I_thresh = 1e-3
options.randomize_surface = True
options.n_jobs = -3
options.periodic = True
options.project_name = "basic_cell_model"
options.nx = 20
options.ny = options.nx
options.n_rays = 4*options.nx ** 2
options.lookuptable_angles = 500
options.depth_spacing = 1e-9
options.depth_spacing_bulk = 1e-9
options.theta_in = 0*np.pi/180
options.phi_in = 0*np.pi/180
options.parallel = True

options.theta_in = 0*np.pi/180


front_materials = [Layer(110e-9, MgF2), Layer(40 * 1e-9, ITO_front),
                   Layer(1e-9, aSi_i), Layer(1e-9, aSi_n)]
back_materials = [Layer(2e-9, aSi_i), Layer(2e-9, aSi_p), Layer(100e-9, ITO_back),
                  Layer(120e-9, MgF2), Layer(150e-9, Ag)]

# thickness of front:
d_front = np.sum([x.width*1e9 for x in front_materials])

tmm_str = tmm_structure(front_materials + [Layer(d_bulk, Si)] + back_materials, incidence=Air, transmission=Air)

options.coherent = False
options.coherency_list = ['c']*len(front_materials) + ['i']*(len(back_materials) + 1)

tmm_result = tmm_str.calculate_profile(options)

bulk_profile = tmm_result['profile']

Si_profile = bulk_profile[:, int(d_front):(int(1e9*d_bulk)+int(d_front)+1)]


# plot integrated absorption - this is what we actually need to compare the EQE to see where losses are
# coming from.

bulk_positions_cm = np.arange(0, d_bulk * 1e2, options.depth_spacing_bulk * 1e2)

# convert bulk_profile from nm-1 to cm-1

profile_func = interp1d(bulk_positions_cm, 1e7 * Si_profile, kind='linear', bounds_error=False, fill_value=0)

L = d_bulk*100 # length of the system in the x-direction [cm]

# Mesh
x = np.concatenate((np.linspace(0,2e-4, 2000, endpoint=False),
                    np.linspace(2e-4, L - 2e-4, 2000, endpoint=False),
                    np.linspace(L - 2e-4, L, 2000)))

A_intgr = np.trapz(profile_func(x), x, axis=1)

Si = material("Si")()

nD = 1e19 # [cm^-3]
nA = 1e20 # [cm^-3]
bulk_doping = 1e15 # p-type

doping_profile = nD * erfc(x/150e-7) # characteristic depth of 150 nm

doping_profile_rear = - nA * erfc((L - x)/200e-7) # characteristic depth of 200 nm

overall_profile = doping_profile + doping_profile_rear - bulk_doping
# Create a system
sys = sesame.Builder(x)

# https://www.pveducation.org/es/fotovoltaica/materials/general-properties-of-silicon

# Dictionary with the material parameters
material_Si = {'Nc': 3e19, 'Nv': 1e19, 'Eg': 1.12, 'affinity': 4.05, 'epsilon': 11.7,
        'mu_e': 1000, 'mu_h': 10000, 'tau_e': 0.05, 'tau_h': 0.05, 'Et': 0,
            'B': 1.89e-15,
            'Cn': 3e-31, 'Cp': 1e-31
            }

# Auger: https://ui.adsabs.harvard.edu/abs/1977ApPhL..31..346D/abstract#:~:text=The%20Auger%20coefficients%20at%20300,31%20cm6%20s%2D1.


# Et: energy levels of the bulk recombination centers [eV]
# B: radiation recombination constant [cm\ :sup:`3`/s], Cn (Cp): Auger
            # recombination constant for electrons (holes) [cm\ :sup:`6`/s],
            # mass_e (mass_h): effective mass of electrons (holes). All parameters
            # can be scalars or callable functions similar to the location
            # argument.
# can all be location dependent



# Add the material to the system
sys.add_material(material_Si)

# fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
# ax1.plot(x, overall_profile, label='erfc')
# # ax1.plot(x, sys.rho*sys.scaling.density, label='flat')
# ax1.set_xlim(0, 400e-7)
# ax1.legend()
#
# ax2.plot(x, overall_profile)
# # ax2.plot(x, sys.rho*sys.scaling.density)
# ax2.set_xlim(L - 600e-7, L)
#
# plt.show()

sys.rho = overall_profile/sys.scaling.density

# Define Ohmic contacts
sys.contact_type('Ohmic', 'Ohmic')

# Define the surface recombination velocities for electrons and holes [cm/s]
Sn_left, Sp_left, Sn_right, Sp_right = 1, 1e3, 1e3, 1
# lower than 1e-6: diverges
sys.contact_S(Sn_left, Sp_left, Sn_right, Sp_right)

def gfcn_fun_rt(wl_index, flux):

    # profile_func is in m-1, but as a function of cm

    def gcfn(x, y):
        return flux*profile_func(x)[wl_index]

    return gcfn

test_rt = gfcn_fun_rt(8, 1e20)
rt_data = test_rt(x, 0)

eqe = np.zeros_like(wls)

for i1, wl in enumerate(wls):

    flux = 1e20
    print(wl)
    sys.generation(gfcn_fun_rt(i1, flux))
    print(gfcn_fun_rt(i1, flux)(0, 0))
    print(gfcn_fun_rt(i1, flux)(1e-6, 0))
    voltages = [0]

    if i1 == 0:
        guess = None

    else:
        guess = {key: result[key][0,:] for key in result.keys()}

    j, result = sesame.IVcurve(sys, voltages, guess=guess)

    # A[i1] = 1 - np.exp(-alpha*L)
    # eqe[i1] = j/(q*flux)

    eqe[i1] = j/(q*flux)


A = tmm_result['A_per_layer'][:, len(front_materials)]
iqe = eqe/A

# convert dimensionless current to dimension-ful current
eqe = eqe * sys.scaling.current
iqe = iqe * sys.scaling.current

plt.figure()
plt.plot(wls*1e9, eqe,'-b', label='EQE')
plt.plot(wls*1e9, iqe,'--r', label='IQE')
plt.plot(wls*1e9, A_intgr, '--k', label='A')
plt.xlabel('Wavelength (nm)')
plt.ylabel('EQE')
plt.legend()
plt.grid()
plt.show()

# def gen_wl(x):
#     return 1e-2*Si.alpha(wls)[None,:]*1e-2 * np.exp(-1e-2*Si.alpha(wls)[None,:]* x[:,
#                                                                               None])

gen_wl = profile_func

# gg = light_source.spectrum()[1][None,:] * gen_wl(x)
# photon flux m-2 / m, generation cm-1
gg = base_light.spectrum(wls, output_units="photon_flux_per_m")[1][:, None] * gen_wl(x)
# g_vs_z = np.trapz(gg, wls, axis=1) / 1e4 # m-2 cm-1 -> cm-3
g_vs_z = (1-shading)*np.trapz(gg, wls, axis=0) / 1e4
g_vs_z[np.isnan(g_vs_z)] = 0


# somehow missing a factor of 100 somewhere (with Beer-Lambert)...
# g_vs_z = 100*g_vs_z

# result: entries are v, efn, and efp. Store in arrays instead.

# equilibrium solution:
j0, result0 = sesame.IVcurve(sys, [0])

voltages = np.concatenate((np.linspace(0,0.55, 8, endpoint=False),
                    np.linspace(0.55, 0.7, 12, endpoint=False),
                    np.linspace(0.7, 0.8, 12)))
sys.generation(g_vs_z)
j, result = sesame.IVcurve(sys, voltages)#, verbose=False)
j = j * sys.scaling.current

jsc = j[0]*1e4

# jsc = q*np.trapz(eqe*light_source.spectrum()[1], wls) # A/m2

zero_crossing = np.where(np.diff(np.sign(j)))[0][0]
j_above = j[zero_crossing]
j_below = j[zero_crossing + 1]

voc = voltages[zero_crossing] + (voltages[zero_crossing + 1] - voltages[zero_crossing]) * j_above / (j_above - j_below)

vmpp = voltages[np.nanargmax(j*voltages)]
jmpp = j[np.nanargmax(j*voltages)]*1e4

FF = vmpp*jmpp/(jsc*voc)

plt.figure()

plt.plot(voltages, 1e4*j, '-k', label='Current')
plt.plot(voltages, 1e4*j*voltages, '--r', label='Power')

plt.xlabel('Voltage [V]')
plt.ylabel(r'Current (A/m$^2$) / Power (W)')
plt.text(0.02, 0.9*jsc, r'Jsc = {:.2f} mA/cm$^2$'.format(jsc/10))
plt.text(0.02, 0.8*jsc, 'Voc = {:.3f} V'.format(voc))
plt.text(0.02, 0.7*jsc, 'FF = {:.3f} %'.format(FF*100))
plt.text(0.02, 0.6*jsc, 'Pmax = {:.3f} W'.format(jmpp*vmpp))
plt.text(0.02, 0.5*jsc, 'Eff = {:.3f} %'.format(jmpp*vmpp/10))
plt.text(0.02, 0.4*jsc, r'Jmpp = {:.3f} mA/cm$^2$'.format(jmpp/10))
plt.text(0.02, 0.3*jsc, 'Vmpp = {:.3f} V'.format(vmpp))

plt.ylim(-100, 1.1e4*j[0])
plt.grid()     # add grid
plt.legend()
plt.show()

# print("Jsc", jsc)
# print("Voc", voc)
# print("eff", np.max(1e4*j*voltages)/1000)
# print("FF", FF)
# print("Vmpp", vmpp)
# print("Jmpp", jmpp)

line = ((0, 0), (L, 0))

# sys, result = sesame.load_sim(f'1dhomo_V_0.gzip')
analyzer = sesame.Analyzer(sys, {key: result[key][0, :] for key in result.keys()})
x, s = analyzer.line(sys, *line)

electron_density_eq = analyzer.electron_density(location=line)*sys.scaling.density
hole_density_eq = analyzer.hole_density(location=line)*sys.scaling.density

result_maxPP = {key: result[key][np.nanargmax(j*voltages), :] for key in result.keys()}
# sys, result = sesame.load_sim(f'1dhomo_V_0.gzip')
analyzer = sesame.Analyzer(sys, result_maxPP)
x, s = analyzer.line(sys, *line)

electron_density = analyzer.electron_density(location=line)
hole_density = analyzer.hole_density(location=line)

srh = analyzer.bulk_srh_rr(location=line)
auger = analyzer.auger_rr(location=line)
rad = analyzer.radiative_rr(location=line)

fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
ax1.semilogy(sys.xpts[s]*1e4, electron_density*sys.scaling.density, label='e-')
ax1.semilogy(sys.xpts[s]*1e4, hole_density*sys.scaling.density, label='h+')
ax1.set_xlim(0, 10)
ax1.set_xlabel('x (um)')
ax1.set_ylabel('density (cm-3)')
ax1.legend()

ax2.semilogy(sys.xpts[s]*1e4, electron_density*sys.scaling.density, label='e-')
ax2.semilogy(sys.xpts[s]*1e4, hole_density*sys.scaling.density, label='h+')
# ax2.plot(x, sys.rho*sys.scaling.density)
ax2.set_xlim(d_bulk*1e6 - 8, d_bulk*1e6)
ax2.set_xlabel('x (um)')
plt.show()


fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
ax1.semilogy(sys.xpts[s]*1e4, electron_density*sys.scaling.density, label='e-')
ax1.semilogy(sys.xpts[s]*1e4, hole_density*sys.scaling.density, label='h+')
ax1.set_xlim(0, 10)
ax1.set_xlabel('x (um)')
ax1.set_ylabel('density (cm-3)')
ax1.legend()

ax2.semilogy(sys.xpts[s]*1e4, electron_density*sys.scaling.density, label='e-')
ax2.semilogy(sys.xpts[s]*1e4, hole_density*sys.scaling.density, label='h+')
# ax2.plot(x, sys.rho*sys.scaling.density)
ax2.set_xlim(d_bulk*1e6 - 8, d_bulk*1e6)
ax2.set_xlabel('x (um)')
plt.show()

fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(5, 3))
ax1.semilogy(sys.xpts[s]*1e4, srh, label='SRH')
ax1.semilogy(sys.xpts[s]*1e4, auger, label='Auger')
ax1.semilogy(sys.xpts[s]*1e4, rad, label='Radiative')
ax1.set_xlim(0, 2)
ax1.set_xlabel('x (um)')
ax1.legend()

ax2.semilogy(sys.xpts[s]*1e4, srh, label='SRH')
ax2.semilogy(sys.xpts[s]*1e4, auger, label='Auger')
ax2.semilogy(sys.xpts[s]*1e4, rad, label='Radiative')
# ax2.plot(x, sys.rho*sys.scaling.density)
ax2.set_xlim(d_bulk*1e6 - 2, d_bulk*1e6)
ax2.set_xlabel('x (um)')
plt.tight_layout()
plt.show()

# analyzer.full_current()*sys.scaling.current * sys.scaling.length

# jsc = analyzer.full_current()*sys.scaling.current

# bulk SRH:
bulk_SRH = analyzer.integrated_bulk_srh_recombination()*sys.scaling.generation * sys.scaling.length
bulk_aug = analyzer.integrated_auger_recombination()*sys.scaling.generation * sys.scaling.length
bulk_rad = analyzer.integrated_radiative_recombination()*sys.scaling.generation * sys.scaling.length

# current_loss = q*(bulk_SRH + bulk_aug + bulk_rad)/jsc

# integrated current:

jsc_max = 1e9*q*np.trapz(A_intgr*base_light.spectrum(wls*1e9)[1], wls)

# current_loss_2 = 1 - jsc*1e4/jsc_max
#
# print(current_loss, current_loss_2)

# account for all the losses

auger_loss = 1e4*q*bulk_aug
srh_loss = 1e4*q*bulk_SRH
rad_loss = 1e4*q*bulk_rad

# current available above 1200 nm:

total_J = 1e9*q*np.trapz(base_light.spectrum(wls*1e9)[1], wls)

R_loss = 1e9*q*np.trapz(tmm_result['R']*base_light.spectrum(wls*1e9)[1], wls)

escape_loss = tmm_result['T']
escape_loss = 1e9*q*np.trapz(escape_loss*base_light.spectrum(wls*1e9)[1], wls)

parasitic_abs = np.sum(tmm_result['A_per_layer'], 1) - A

parasitic_loss = 1e9*q*np.trapz(parasitic_abs*base_light.spectrum(wls*1e9)[1], wls)
# units of SRV: cm/s
# units of carrier density: cm-3
# SRV*carrier density = cm-2/s
# q*SRV*carrier density = q*cm-2/s = A cm-2

shading_loss = 1e9*q*np.trapz(shading*base_light.spectrum(wls*1e9)[1], wls)

# think this should be excess carrier density
n_front_eq = electron_density_eq[0]
p_rear_eq = hole_density_eq[-1]

n_rear_eq = electron_density_eq[-1]
p_front_eq = hole_density_eq[0]

e_sr_left = Sn_left*1e4*q*(electron_density[0]*sys.scaling.density - n_front_eq)
h_sr_left = Sp_left*1e4*q*(hole_density[0]*sys.scaling.density - p_front_eq)
e_sr_right = Sn_right*1e4*q*(electron_density[-1]*sys.scaling.density - n_rear_eq)
h_sr_right = Sp_right*1e4*q*(hole_density[-1]*sys.scaling.density - p_rear_eq)

# only minority carrier recombination matters
sr_left = e_sr_left if hole_density[0]*sys.scaling.density > electron_density[0]*sys.scaling.density else h_sr_left
sr_right = e_sr_right if hole_density[-1]*sys.scaling.density > electron_density[-1]*sys.scaling.density else h_sr_right

print(shading_loss + R_loss + parasitic_loss + escape_loss + jmpp + auger_loss + srh_loss + rad_loss + sr_left + sr_right)
print(total_J)

J_list = [jmpp, shading_loss, R_loss, parasitic_loss, escape_loss, auger_loss, srh_loss, rad_loss, sr_left, sr_right]
fig, ax = plt.subplots(figsize=(4, 6))

labels = [r'$J_{MPP}$', 'Shading', r'R$_{front}$', r'A$_{parasitic}$', 'Escape', 'Auger', 'SRH', 'Radiative', 'Front surface', 'Rear surface']

for j1, item in enumerate(J_list):

    ax.bar(['J'], item, bottom=np.sum(J_list[:j1]), label=labels[j1])

ax.set_ylim(300, 1.01*total_J)
ax.axhline(total_J, color='k', linestyle='--', label='Available current')
plt.legend()
plt.show()

print(total_J - np.sum(J_list))

# LONGi:
# Jsc: 41.16
# Voc: 0.751
# FF: 86.5

# 33.91
# 0.777
# 22.766