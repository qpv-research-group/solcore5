import numpy as np
import os
from solcore import siUnits, material, eVnm, si
from solcore.constants import h, c, vacuum_permittivity
from solcore.structure import Structure, Junction, Layer
from solcore.state import State
from scipy.interpolate import interp1d

# ********************************
# Triple junction cell QE calculation

all_materials = []

# We need to build the solar cell layer by layer.
# We start from the AR coating. In this case, we load it from an an external file
def this_dir_file(f):
    return os.path.join(os.path.split(__file__)[0], f)

refl_nm = np.loadtxt(this_dir_file("MgF-ZnS_AR.csv"), unpack=True, delimiter=",")
refl_J = np.array((h * c / siUnits(refl_nm[0], "nm")[::-1], refl_nm[1][::-1]))
ref = interp1d(siUnits(refl_nm[0], "nm"), refl_nm[1], bounds_error=False, fill_value=(0, 0))


# Next is the window layer, made of AlInP. We load the absorption coefficent of AlInP from an external file
AlInP = material("AlInP")
window_material = AlInP(Al=0.52)
EeV, alpha = np.loadtxt(this_dir_file("alinp.csv"), unpack=True, delimiter=",")
WLm = siUnits(eVnm(EeV), 'nm')
window_material.alpha = interp1d(x=WLm[::-1], y=alpha[::-1], bounds_error=False, fill_value=0)

all_materials.append(window_material)


# TOP CELL - GaInP
# Now we build the top cell, which requires the n and p sides of GaInP.
# We also load the absorption coefficient from an external file. We also add some extra parameters needed for the
# calculation such as the minority carriers diffusion lengths
InGaP = material("GaInP")
top_cell_n_material = InGaP(In=0.49, Nd=siUnits(2e18, "cm-3"), role="n")
top_cell_p_material = InGaP(In=0.49, Na=siUnits(1e17, "cm-3"), role="p")
EeV, alpha = np.loadtxt(this_dir_file("in048ga052p.csv"), unpack=True, delimiter=",")
WLm = siUnits(eVnm(EeV), 'nm')
top_cell_n_material.alpha = interp1d(x=WLm[::-1], y=alpha[::-1], bounds_error=False, fill_value=0)
top_cell_p_material.alpha = top_cell_n_material.alpha
top_cell_n_material.hole_diffusion_length = si("200nm")
top_cell_p_material.electron_diffusion_length = si("1um")

all_materials.append(top_cell_n_material)
all_materials.append(top_cell_p_material)


# MID CELL  - InGaAs
# The same thing.
InGaAs = material("GaAs")
mid_cell_n_material = InGaAs(In=0.01, Nd=siUnits(3e18, "cm-3"), role="n")
mid_cell_p_material = InGaAs(In=0.01, Na=siUnits(1e17, "cm-3"), role="p")
EeV, alpha = np.loadtxt(this_dir_file("in01gaas.csv"), unpack=True, delimiter=",")
WLm = siUnits(eVnm(EeV), 'nm')
mid_cell_n_material.alpha = interp1d(x=WLm[::-1], y=alpha[::-1], bounds_error=False, fill_value=0)
mid_cell_p_material.alpha = mid_cell_n_material.alpha
mid_cell_n_material.hole_diffusion_length = si("500nm")
mid_cell_p_material.electron_diffusion_length = si("5um")

all_materials.append(mid_cell_n_material)
all_materials.append(mid_cell_p_material)

# BOTTOM CELL - Levinshtein
# Idem
Ge = material("Ge")
bot_cell_n_material = Ge(Nd=siUnits(2e18, "cm-3"), role="n")
bot_cell_p_material = Ge(Na=siUnits(1e17, "cm-3"), role="p")
EeV, alpha = np.loadtxt(this_dir_file("Ge-Palik.csv"), unpack=True, delimiter=",")
WLm = siUnits(eVnm(EeV), 'nm')
bot_cell_n_material.alpha = interp1d(x=WLm[::-1], y=alpha[::-1], bounds_error=False, fill_value=0)
bot_cell_p_material.alpha = bot_cell_n_material.alpha
bot_cell_n_material.hole_diffusion_length = si("800nm")
bot_cell_p_material.electron_diffusion_length = si("50um")

all_materials.append(bot_cell_n_material)
all_materials.append(bot_cell_p_material)


# We add some other properties to the materials, assumed the same in all cases.
# If different, we should have added them above.
for mat in all_materials:
    mat.hole_mobility = 3.4e-3
    mat.electron_mobility = 5e-2
    mat.permittivity = 9


# And, finally, we put everithing together, adding also the surface recombination velocities.
triplejunction = Structure(
    [
        Layer(material=window_material, width=si("25nm")),
        Junction(
            (
                Layer(si("100nm"), material=top_cell_n_material),
                Layer(si("600nm"), material=top_cell_p_material),
            ),
            sn=1,
            sp=1,
        ),
        Junction(
            (
                Layer(si("100nm"), material=mid_cell_n_material),
                Layer(si("3.5um"), material=mid_cell_p_material),
            ),
            sn=1,
            sp=1,
        ),
        Junction(
            (
                Layer(si("400nm"), material=bot_cell_n_material),
                Layer(si("100um"), material=bot_cell_p_material),
            ),
            sn=1,
            sp=1,
        ),
    ]
)

# Now the 'solar_cell' object below contains all the information regarding the solar cell. So, we use it
solar_cell = State()
solar_cell.structure = triplejunction
solar_cell.reflectivity = refl_J
