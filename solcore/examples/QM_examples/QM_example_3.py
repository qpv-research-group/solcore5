"""
We use the same two options for solving QWs to calculate the absorption coefficient. In the case of mode='kp4x4',
the correct 1D absorption calculator is not implemented yet, so it uses the parabolic aproximation to calculate the
effective mass at k=0... In other words, it is wrong and, for now, inconsistent, so use carefully.
"""
import numpy as np
import matplotlib.pyplot as plt
from solcore import material, si
from solcore.structure import Layer
from solcore.constants import vacuum_permittivity, q
import solcore.quantum_mechanics as QM
from solcore import si, asUnit

GaAs = material("GaAs")
InGaAs = material("InGaAs")
GaAsP = material("GaAsP")
GaP = material("GaP")

InGaAs.strained = True
GaAsP.strained = True
GaP.strained = True
GaAs.strained = False

bulkMaterial = GaAs(T=300)
wellMaterial = InGaAs(In=0.245, T=300)
barrierMaterial = GaAsP(P=0.1, T=300)
interlayer_material = InGaAs(In=0.14, T=300)

my_structure = QM.assemble_qw_structure(
    repeats=1,
    well=Layer(si("7.2nm"), wellMaterial),
    bulk_l_top=Layer(si("1nm"), barrierMaterial),
    bulk_l_bottom=Layer(si("1nm"), barrierMaterial),
    barrier=Layer(si("28nm"), barrierMaterial),
    well_interlayer=Layer(si("3nm"), interlayer_material)
)
my_structure.substrate = bulkMaterial

# We have to define a few extra parameters allowing the solver to calculate the absorption.
wl = np.linspace(800, 1300, 1000)
E = 1240 / wl * q

alpha_params = {
    "well_width": my_structure.width(),
    "theta": 0,
    "eps": 12.9 * vacuum_permittivity,
    "espace": E,
    "hwhm": si("4meV"),
    "dimensionality": 0.2,
}

band_edge_kp8x8, b = QM.schrodinger(my_structure, mode='kp8x8_bulk', calculate_absorption=True,
                                    alpha_params=alpha_params)
band_edge_kp4x4, b = QM.schrodinger(my_structure, mode='kp4x4', calculate_absorption=True, alpha_params=alpha_params)

# You will notice the big difference between both (a factor of 2, more or less) related with the uncorrect treatement
# of the kp4x4 case. The number of steps in the curve is related with the number of energy levels that each of the
# method produces, and therefor, the number of possible trnsitions.
plt.plot(asUnit(band_edge_kp8x8['alpha'][0], 'eV'), asUnit(band_edge_kp8x8['alpha'][1], 'cm-1'), 'b', label='kp8x8_bulk')
plt.plot(asUnit(band_edge_kp4x4['alpha'][0], 'eV'), asUnit(band_edge_kp4x4['alpha'][1], 'cm-1'), 'r', label='kp4x4')
plt.legend()
plt.show()