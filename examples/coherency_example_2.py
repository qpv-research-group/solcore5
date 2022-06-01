import numpy as np
import matplotlib.pyplot as plt

from solcore import siUnits, material, si
from solcore.solar_cell import SolarCell
from solcore.structure import Junction, Layer
from solcore.solar_cell_solver import solar_cell_solver, default_options
from solcore.light_source import LightSource
from solcore.state import State


incidence_angle = 45 # should be in degrees
wl = np.linspace(290, 1900, 400) * 1e-9
concX=566 # the light concentration
light_source = LightSource(source_type='standard', version='AM1.5d', x=default_options.wavelength,
                           output_units='photon_flux_per_m', concentration=concX) # define the input light source as AM1.5G

all_materials = []

Al2O3 = material('Al2O3')
TiO2 = material('TiO2')
AlInP = material("AlInP")
GaInP = material("GaInP")
GaAs = material('GaAs')
Ge = material("Ge")
Al02Ga08As = material('AlGaAs')
Al08Ga02As = material('AlGaAs')

# TOP CELL - GaInP
ARC1= Al2O3()
ARC2 = TiO2()

top_window_material = AlInP(Al=0.5)
top_cell_n_material = GaInP(In=0.51,Nd=siUnits(2e18, "cm-3"), hole_diffusion_length=si("300nm"))
top_cell_p_material = GaInP(In=0.51,Na=siUnits(1.5e17, "cm-3"), electron_diffusion_length=si("2um"))
top_cell_TJ_material = Al08Ga02As(Al=0.8)

for mat in [top_cell_n_material, top_cell_p_material]:
    mat.band_gap = material('GaInP')(In=0.51).band_gap
    mat.eff_mass_hh_z = material('GaInP')(In=0.51).eff_mass_hh_z
    mat.eff_mass_electron = material('GaInP')(In=0.51).eff_mass_electron

all_materials.append(ARC1)
all_materials.append(ARC2)
all_materials.append(top_window_material)
all_materials.append(top_cell_n_material)
all_materials.append(top_cell_p_material)
all_materials.append(top_cell_TJ_material)

# MID CELL  - InGaAs
mid_window_material = GaInP(In=0.51)
mid_cell_n_material = GaAs(Nd=siUnits(2e18, "cm-3"), hole_diffusion_length=si("300nm"))
mid_cell_p_material = GaAs(Na=siUnits(1.5e17, "cm-3"), electron_diffusion_length=si("3um"))
mid_BSF_material = GaInP(In=0.51)
mid_cell_TJ_material = Al08Ga02As(Al=0.8)

for mat in [mid_cell_n_material, mid_cell_p_material]:
    mat.band_gap = material('GaAs')(In=0.01).band_gap
    mat.eff_mass_hh_z = material('GaAs')(In=0.01).eff_mass_hh_z
    mat.eff_mass_electron = material('GaAs')(In=0.01).eff_mass_electron

all_materials.append(mid_window_material)
all_materials.append(mid_cell_n_material)
all_materials.append(mid_cell_p_material)
all_materials.append(mid_BSF_material)
all_materials.append(mid_cell_TJ_material)

DBR1 = Al02Ga08As(Al=0.2)
DBR2 = Al08Ga02As(Al=0.8)

all_materials.append(DBR1)
all_materials.append(DBR2)

# BOTTOM CELL - Ge

bot_buffer_material = GaAs()
bot_nucleation_material = GaInP(In=0.51)
bot_cell_n_material = Ge(Nd=siUnits(2e18, "cm-3"), hole_diffusion_length=si("800nm"))
bot_cell_p_material = Ge(Na=siUnits(1e17, "cm-3"), electron_diffusion_length=si("50um"))

for mat in [bot_cell_n_material, bot_cell_p_material]:
    mat.band_gap = material('Ge')().band_gap
    mat.eff_mass_hh_z = material('Ge')().eff_mass_hh_z
    mat.eff_mass_electron = material('Ge')().eff_mass_electron

all_materials.append(bot_buffer_material)
all_materials.append(bot_nucleation_material)
all_materials.append(bot_cell_n_material)
all_materials.append(bot_cell_p_material)


# We add some other properties to the materials, assumed the same in all cases, for simplicity.
# If different, we should have added them above in the definition of the materials.

for mat in all_materials:

    mat.hole_mobility = 3.4e-3
    mat.electron_mobility = 5e-2


ARC = [Layer(si('80nm'), material = ARC1), Layer(si('33nm'), material = ARC2)]

top_junction = [Junction([Layer(si("18nm"), material=top_window_material, role='window'),
                  Layer(si("100nm"), material=top_cell_n_material, role='emitter'),
                  Layer(si("891.248nm"), material=top_cell_p_material, role='base'),
                  Layer(si("111.445nm"), material = top_cell_TJ_material, role = 'TJ')
                 ], sn=1, sp=1, kind='DA')]
middle_junction = [Junction([Layer(si("18nm"), material=mid_window_material, role='window'),
                  Layer(si("100nm"), material=mid_cell_n_material, role='emitter'),
                  Layer(si("1632.091nm"), material=mid_cell_p_material, role='base'),
                  Layer(si("10nm"), material = mid_BSF_material, role = 'BSF'),
                  Layer(si("91.084nm"), material=mid_cell_TJ_material, role='TJ')
                 ], sn=1, sp=1, kind='DA')]
DBRa = 16 * [Layer(width=si("62.638nm"), material=DBR1), Layer(width=si("71.980nm"), material=DBR2)]
DBRb = 16 * [Layer(width=si("68.919nm"), material=DBR1), Layer(width=si("78.725nm"), material=DBR2)]
DBRc = 16 * [Layer(width=si("75.838nm"), material=DBR1), Layer(width=si("86.805nm"), material=DBR2)]
        # the 4* here makes the two layers given repeat 4 times (so 8 layers total)
bottom_junction = [Junction([Layer(si("405.048nm"), material=bot_buffer_material, role='window'),
                  Layer(si("14.369nm"), material=bot_nucleation_material, role='window'),
                  Layer(si("200nm"), material=bot_cell_n_material, role='emitter'),
                  Layer(si("29800nm"), material = bot_cell_p_material, role = 'base')
                 ], sn=1, sp=1, kind='DA')]
# And, finally, we put everything together, adding also the surface recombination velocities sn and sp.
# setting kind = 'DA' in the Junction definition tells the electrical solver later to use the depletion approximation
optical_struct = SolarCell(ARC + top_junction + middle_junction + DBRa + DBRb + DBRc + bottom_junction,
                           shading = 0.05)


wl = np.linspace(250, 1700, 400)*1e-9


options = State()
options.wavelength = wl
options.optics_method = 'TMM'
options.no_back_reflection = False
options.pol = 'p'
options.BL_correction = True
options.coherency_list = 111*['c']
options.theta = 30
solar_cell_solver(optical_struct, 'qe', options)

plt.figure()
plt.plot(wl*1e9, optical_struct[0].layer_absorption+optical_struct[1].layer_absorption)
plt.plot(wl*1e9, optical_struct[2].layer_absorption)
plt.plot(wl*1e9, optical_struct[3].layer_absorption)
plt.plot(wl*1e9, optical_struct[100].layer_absorption)
plt.plot(wl*1e9, optical_struct.absorbed, '--')
plt.plot(wl*1e9, optical_struct.transmitted, '--')
plt.plot(wl*1e9, optical_struct.reflected, '--')
plt.legend(['ARC', 'top', 'middle', 'bottom', 'A', 'T', 'R'])
plt.ylim(0,1)
plt.ylabel('Absorption/Transmission/Reflection')
plt.xlabel('Wavelength (nm)')
plt.show()

plt.figure()
plt.plot(wl*1e9, 100*optical_struct[2].eqe(wl))
plt.plot(wl*1e9, 100*optical_struct[3].eqe(wl))
plt.plot(wl*1e9, 100*optical_struct[100].eqe(wl))
plt.plot(wl*1e9, 100*optical_struct.absorbed, '--')
plt.legend(['top', 'middle', 'bottom', 'A'])
plt.ylim(0,100)
plt.ylabel('EQE (%)')
plt.xlabel('Wavelength (nm)')
plt.show()
