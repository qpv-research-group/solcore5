"""Reading in material parameters and using the modified critical point parabolic band model (CPPB) to approximate the complex dielectric
function of GaAs. This model is based on:

- Adachi, S., Physical Properties of III-V Semiconductor Compounds, John Wiley & Sons (1992)
- C. C. Kim et al, 'Modelling the optical dielectric function of semiconductors: Extension of the critical-point parabolic band approximation', Physical Review B 45(20) 11749, 1992

"""

import matplotlib.pyplot as plt
import numpy as np

from solcore.graphing.Custom_Colours import colours
from solcore.absorption_calculator.cppm import Custom_CPPB as cppb
from solcore.absorption_calculator.dielectric_constant_models import Oscillator
from solcore.structure import Structure

# First, read in experimental data for GaAs dielectric function (from Palik)...
Palik_Eps1 = np.loadtxt("data/Palik_GaAs_Eps1.csv", delimiter=',', unpack=False)
Palik_Eps2 = np.loadtxt("data/Palik_GaAs_Eps2.csv", delimiter=',', unpack=False)

# Generate a list of energies over which to calculate the model dielectric function.
E = np.linspace(0.2, 5, 1000)

# Class object is created, CPPB_Model
CPPB_Model = cppb()

# The Material_Params method loads in the desired material parameters as a dictionary variable...
MatParams = CPPB_Model.Material_Params("GaAs")

# Parameters can be customised by assigning to the correct dictionary key...
MatParams["B1"] = 5.8
MatParams["B1s"] = 1.0
MatParams["Gamma_Eg_ID"] = 0.3
MatParams["Alpha_Eg_ID"] = 0.0
MatParams["E1"] = 2.8
MatParams["E1_d1"] = 2.9
MatParams["Gamma_E1"] = 0.1
MatParams["E2"] = 4.72
MatParams["C"] = 3.0
MatParams["Alpha_E2"] = 0.04
MatParams["Gamma_E2"] = 0.19

# Must define a structure object containing the required oscillator functions. The oscillator type and material
# parameters are both passed to individual 'Oscillators' in the structure...
Adachi_GaAs = Structure([
    Oscillator(oscillator_type="E0andE0_d0", material_parameters=MatParams),
    Oscillator(oscillator_type="E1andE1_d1", material_parameters=MatParams),
    Oscillator(oscillator_type="E_ID", material_parameters=MatParams),
    Oscillator(oscillator_type="E2", material_parameters=MatParams)
])

Output = CPPB_Model.eps_calc(Adachi_GaAs, E)

# PLOT OUTPUT...
fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(7, 4.5))

# Subplot I :: Real part of the dielectric function.
ax1.set_yscale("linear")
ax1.set_xlim(0, 5.3)
ax1.set_ylim(-14, 27)

ax1.plot(Palik_Eps1[:, 0], Palik_Eps1[:, 1], label="Exp. Data (Palik)",
         marker='o', ls='none', markerfacecolor='none', markeredgecolor=colours("Red"))

ax1.plot(E, Output["eps"].real, color=colours("Navy"), label="Total")
ax1.plot(E, Output["components"][0].real, color=colours("Orange Red"), ls='--', label="$E_0$ and $E_0+\Delta_0$")
ax1.plot(E, Output["components"][1].real, color=colours("Dodger Blue"), ls='--', label="$E_1$ and $E_1+\Delta_1$")
ax1.plot(E, Output["components"][2].real, color=colours("lime green"), ls='--', label="$E_{ID}$ (Indirect)")
ax1.plot(E, Output["components"][3].real, color=colours("gold"), ls='--', label="$E_2$")

ax1.set_xlabel("Energy (eV)")
ax1.set_ylabel("$\epsilon_1 (\omega)$")
ax1.text(0.05, 0.05, '(a)', transform=ax1.transAxes, fontsize=12)

# Subplot II :: Imaginary part of the dielectric function.
ax2.set_yscale("linear")

ax2.plot(Palik_Eps2[:, 0], Palik_Eps2[:, 1], label="Exp. Data (Palik)",
         marker='o', ls='none', markerfacecolor='none', markeredgecolor=colours("Red"))

ax2.plot(E, Output["eps"].imag, color=colours("Navy"), label="Total")
ax2.plot(E, Output["components"][0].imag, color=colours("Orange Red"), ls='--', label="$E_0$ and $E_0+\Delta_0$")
ax2.plot(E, Output["components"][1].imag, color=colours("Dodger Blue"), ls='--', label="$E_1$ and $E_1+\Delta_1$")
ax2.plot(E, Output["components"][2].imag, color=colours("lime green"), ls='--', label="$E_{ID}$ (Indirect)")
ax2.plot(E, Output["components"][3].imag, color=colours("gold"), ls='--', label="$E_2$")
ax2.set_xlim(0, 5.3)
ax2.set_ylim(0, 27)

ax2.set_xlabel("Energy (eV)")
ax2.set_ylabel("$\epsilon_2 (\omega)$")
ax2.text(0.05, 0.05, '(b)', transform=ax2.transAxes, fontsize=12)
ax2.legend(loc="upper left", frameon=False)

plt.tight_layout()
plt.show()
