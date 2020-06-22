from scipy.interpolate import InterpolatedUnivariateSpline
import matplotlib.pyplot as plt
import numpy as np

from solcore.absorption_calculator import calculate_ellipsometry
from solcore.structure import Structure
from solcore.data_analysis_tools.ellipsometry_analysis import EllipsometryData
from solcore.graphing.Custom_Colours import colours
from solcore.absorption_calculator.cppm import Custom_CPPB as cppb
from solcore.absorption_calculator.dielectric_constant_models import Oscillator

E_eV = np.linspace(0.7, 4.2, 1000)

# Load in ellipsomery data from file...
Exp_Data = EllipsometryData("data/ge_ellipsometry_data.dat")
Exp_Angles = Exp_Data.angles

# Load in some experimental Ge n-k to compare fit with this...
Ge_nk_Exp = np.loadtxt("data/Ge_nk.csv", delimiter=",", unpack=False)

# Smooth the data with spline fitting...
n_spline = InterpolatedUnivariateSpline(x=Ge_nk_Exp[::5, 0], y=Ge_nk_Exp[::5, 1], k=3)(E_eV)
k_spline = InterpolatedUnivariateSpline(x=Ge_nk_Exp[::5, 2], y=Ge_nk_Exp[::5, 3], k=3)(E_eV)

## Step 1 :: n and k modelling...
# First model the Ge02 layer with the Sellmeier model

# Define Oscillator Structure
GeO2 = Structure([
    Oscillator(oscillator_type="Sellmeier", material_parameters=None,
               A1=0.80686642, L1=0.68972606E-1,
               A2=0.71815848, L2=0.15396605,
               A3=0.85416831, L3=0.11841931E2)
])

GeO2_nk = cppb().nk_calc(oscillator_structure=GeO2, energy_array=E_eV)

# Step 2 :: use this modelled n and k to calculate the ellipsometry data...
# Define a structure for the optical stack...
stack = Structure([
    [4.4, 1240 / E_eV, GeO2_nk["n"], GeO2_nk["k"]],  # Layer 1 :: GeO2 native oxide layer
    [350000, 1240 / E_eV, n_spline, k_spline]  # Layer 2/ Substrate :: Bulk Ge
])

# Calculate Ellipsometry data...
Out = calculate_ellipsometry(stack, 1240 / E_eV, angle=Exp_Angles)

# Define functions for the quick conversion of data
# We show this with the angle = 79º, which is the third one (i = 2)
i = 2

rho = lambda psi, delta: np.tan(psi) * np.exp(1j * delta)
eps = lambda r, theta: np.sin(theta) ** 2 * (1 + np.tan(theta) ** 2 * ((1 - r) / (1 + r)) ** 2)

# Experimental data...
Exp_rho = rho(np.radians(Exp_Data.data[Exp_Angles[i]][1]), np.radians((Exp_Data.data[Exp_Angles[i]][3])))
Exp_eps = eps(Exp_rho, np.radians(Exp_Angles[i]))

# Modelled data...
Mod_rho = rho(np.radians(Out["psi"][:, i]), np.radians(Out["Delta"][:, i]))
Mod_eps = eps(Mod_rho, np.radians(Exp_Angles[i]))

## Step 3 :: Data Plotting...
fig, ax1 = plt.subplots(1, 1)
ax1b = ax1.twinx()
ax1.set_xlim([400, 1500])

# Experimental data...
ax1.plot(Exp_Data.data[Exp_Angles[i]][0] * 1000, Exp_eps.real, lw=2, marker="o", ls='none', color=colours("Orange Red"),
         label="$\epsilon_1 (\lambda)$ :: $ %3.1f^{\circ}$" % Exp_Angles[i])
ax1b.plot(Exp_Data.data[Exp_Angles[i]][0] * 1000, abs(Exp_eps.imag), lw=2, marker="s", ls='none',
          color=colours("Dodger Blue"),
          label="$\epsilon_2 (\lambda)$ :: $ %3.1f^{\circ}$" % Exp_Angles[i])

# Modelled Data...
ax1.plot(1240 / E_eV, Mod_eps.real, label="Model $\epsilon_1 (\lambda)$ :: $ %3.1f^{\circ}$" % Exp_Angles[i],
         color=colours("Maroon"))
ax1b.plot(1240 / E_eV, abs(Mod_eps.imag), label="Model $\epsilon_2 (\lambda)$ :: $ %3.1f^{\circ}$" % Exp_Angles[i],
          color=colours("Navy"))

ax1.set_xlabel("Wavelength (nm)")
ax1.set_ylabel('$\epsilon_1 (\lambda)$')
ax1b.set_ylabel('$\epsilon_2 (\lambda)$')
ax1.text(0.05, 0.9, '(b)', transform=ax1.transAxes, fontsize=12)

ax1.legend(loc="lower left")
ax1b.legend(loc="upper right")

plt.show()
