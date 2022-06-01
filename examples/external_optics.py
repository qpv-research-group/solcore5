from solcore import si
from solcore import material
from solcore.solar_cell import SolarCell, Layer, Junction
from solcore.solar_cell_solver import solar_cell_solver, default_options
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

T = 300
wavelengths_optics = np.linspace(300, 1200, 800)*1e-9

Si = material("Si")
SiO2 = material("SiO2")()

n_material = Si(T=T, Nd=si(1e21, "cm-3"), hole_diffusion_length=si("50um"),
                electron_mobility=50e-4)
p_material = Si(T=T, Na=si(1e16, "cm-3"), electron_diffusion_length=si("150um"),
                hole_mobility=400e-4)

ARC_width = si("100nm")
n_material_width = si("500nm")
p_material_width = si("50um")

solar_cell = SolarCell(
    [
        Layer(width=ARC_width, material=SiO2),
        Junction([Layer(width=n_material_width, material=n_material, role='emitter'),
                  Layer(width=p_material_width, material=p_material, role='base'),
		 ], sn=1, sp=1, kind='DA'),
    ])

total_width = ARC_width + n_material_width + p_material_width

options = default_options
options.optics_method = "TMM"
options.wavelength = wavelengths_optics
# options.position = np.linspace(0, total_width, 100000)
options.light_iv = True
V = np.linspace(0, 1.2, 200)

solar_cell_solver(solar_cell, 'iv', options)

reflected = solar_cell.reflected
absorbed_in_Si = solar_cell[1].layer_absorption

interp_ref = interp1d(options.wavelength, reflected)
interp_totalA = interp1d(options.wavelength, solar_cell[1].layer_absorption)

wavelengths_external = np.linspace(301, 1199, 800)*1e-9

alpha = n_material.alpha(wavelengths_external)
A_layer = interp_totalA(wavelengths_external)

junction_width = n_material_width + p_material_width

def make_absorb_fn(alpha, A_layer, junction_width):
    norm = A_layer * alpha / (1 - np.exp(-alpha * junction_width))
    def profile(z):
        xy = norm[None, :] * np.exp(-alpha * z[:, None])
        return xy.T
    return profile

diff_absorb_fn = make_absorb_fn(alpha, A_layer, junction_width)

solar_cell_external = SolarCell(
    [
        Junction([Layer(width=n_material_width, material=n_material, role='emitter'),
                  Layer(width=p_material_width, material=p_material, role='base'),
		 ], sn=1, sp=1, kind='DA'),
    ], external_reflected=interp_ref(wavelengths_external), external_absorbed=diff_absorb_fn)

options.optics_method = "external"
options.wavelength = wavelengths_external

solar_cell_solver(solar_cell_external, 'iv', options)

plt.figure()
plt.plot(wavelengths_external*1e9, solar_cell_external.reflected, label='Reflected - external')
plt.plot(wavelengths_external*1e9, solar_cell_external.absorbed, label='Absorbed - external')
plt.plot(wavelengths_optics*1e9, reflected, 'k--', label='Reflected - TMM')
plt.plot(wavelengths_optics*1e9, absorbed_in_Si, '--', label='Absorbed - TMM')
plt.legend()
plt.xlabel('Wavelength (nm)')
plt.ylabel('R/A')
plt.show()

plt.figure(1)
plt.plot(V, -solar_cell[1].iv(V), 'b', label='TMM calculation')
plt.plot(V, -solar_cell_external[0].iv(V), 'k--', label='External optics')
plt.legend()
plt.ylim(-20, 350)
plt.xlim(0, 1)
plt.ylabel('Current (A/m$^2$)')
plt.xlabel('Voltage (V)') #The expected values of Isc and Voc are 372 A/m^2 and 0.63 V respectively

plt.show()

position_plot = np.linspace(0, 200, 100)*1e-9

absorption_profile_TMM = solar_cell[0].diff_absorption(position_plot)
absorption_profile_constructed = solar_cell_external[0].diff_absorption(position_plot)

plt.figure()
plt.plot(position_plot, solar_cell[1].absorbed(position_plot))
plt.show()

plt.figure(figsize=(10, 4))
plt.subplot(121)
plt.imshow(absorption_profile_TMM, aspect='auto', extent=[0, 400, 1200, 300])
plt.colorbar()
plt.xlabel('Depth (nm)')
plt.ylabel('Wavelength (nm)')
plt.title(r'Differential absorption (m$^{-1}$) in front surface of cell')

plt.subplot(122)
plt.imshow(absorption_profile_constructed, aspect='auto', extent=[0, 400, 1200, 300])
plt.xlabel('Depth (nm)')
plt.ylabel('Wavelength (nm)')
plt.show()