from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import numpy as np

from solcore.absorption_calculator import calculate_ellipsometry
from solcore.structure import Structure
from solcore.data_analysis_tools.ellipsometry_analysis import EllipsometryData
from solcore.graphing.Custom_Colours import colours
from solcore.absorption_calculator.cppm import Custom_CPPB as cppb
from solcore.absorption_calculator.dielectric_constant_models import Oscillator
from solcore.absorption_calculator import calculate_rat


fig, (ax2, ax1) = plt.subplots(1, 2, figsize=(11.25,4))

# Figure 10 (a)

E_eV = np.linspace(0.65, 4, 1000)

def nk_convert(fname, energy):
    """ Designed to handle nk data files"""

    # Import data...
    E_n, n, E_k, k = np.loadtxt(fname, delimiter=",", unpack=True)
    print("File :: " + fname + " :: Imported Successfully!")

    # Interpolate data...
    n_interp = interp1d(E_n, n, bounds_error=False, fill_value=(n[0], n[-1]))(energy)
    k_interp = interp1d(E_k, k, bounds_error=False, fill_value=(k[0], k[-1]))(energy)

    return (energy, n_interp, k_interp)


# Load nk data using nk_convert function...
alinp_nk = nk_convert("AlInP_nk.csv", energy=E_eV)
gainp_nk = nk_convert("GaInP_nk.csv", energy=E_eV)
mgf_nk = nk_convert("MgF_nk.csv", energy=E_eV)
sic_nk = nk_convert("SiCx_nk.csv", energy=E_eV)
zns_nk = nk_convert("ZnS_nk.csv", energy=E_eV)

# Build the optical stack...
OptiStack = Structure([
    [117, 1240/E_eV, mgf_nk[1], mgf_nk[2]],
    [80, 1240/E_eV, sic_nk[1], sic_nk[2]],
    [61, 1240/E_eV, zns_nk[1], zns_nk[2]],
    [25, 1240/E_eV, alinp_nk[1], alinp_nk[2]],
    [350000, 1240/E_eV, gainp_nk[1], gainp_nk[2]]
])

angles = np.linspace(0, 80, 10)

RAT_angles = []

print("Calculate RAT ::")
for theta in angles:

    rat_data = []
    print("Calculating at angle :: %4.1f deg" % theta)
    # Calculate RAT data...
    rat_data = calculate_rat(OptiStack, angle=theta, wavelength=1240 / E_eV)

    RAT_angles.append((theta, rat_data["R"], rat_data["A"]))

colors = plt.cm.jet(np.linspace(1,0,len(RAT_angles)))

ax2.set_ylim([0, 100])
ax2.set_xlim([300, 1800])

for i, RAT in enumerate(RAT_angles):

    ax2.plot(1240/E_eV, RAT[1]*100, ls="-", color=colors[i], label="%4.1f$^\circ$" % RAT[0])
    ax2.plot(1240/E_eV, RAT[2]*100, ls="--", color=colors[i])
    # plt.plot(1240/E_eV, RAT[1] + RAT[2], lw=3, ls="--", color=cmap_col[i])

ax2.set_xlabel("Wavelength (nm)")
ax2.set_ylabel("Reflection and Transmission (%)")
ax2.legend(loc=5)
ax2.text(0.05, 0.45, '(a)', transform=ax2.transAxes, fontsize=12)

plt.tight_layout(w_pad=4)

###########
# Figure 10 (b)

E_eV = np.linspace(0.7, 4.2, 1000)

# Load in ellipsomery data from file...
Exp_Data = EllipsometryData("ge_sylarus_wafer_scan_300_1700nm_75_77_79_alignment_out_ar.dat")
Exp_Angles = Exp_Data.angles

# Load in some experimental Ge n-k to compare fit with this...
Ge_nk_Exp = np.loadtxt("Ge_nk.csv", delimiter=",", unpack=False)

# Smooth the data with spline fitting...
n_spline = InterpolatedUnivariateSpline(x=Ge_nk_Exp[::5,0], y=Ge_nk_Exp[::5,1], k=3)(E_eV)
k_spline = InterpolatedUnivariateSpline(x=Ge_nk_Exp[::5,2], y=Ge_nk_Exp[::5,3], k=3)(E_eV)

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

## Step 2 :: use this modelled n and k to calculate the ellipsometry data...
# Define a structure for the optical stack...
OptiStack = Structure([
    [4.4, 1240/E_eV, GeO2_nk["n"], GeO2_nk["k"]], # Layer 1 :: GeO2 native oxide layer
    [350000, 1240/E_eV, n_spline, k_spline]      # Layer 2/ Substrate :: Bulk Ge
])

# Calculate Ellipsometry data...
Out = calculate_ellipsometry(OptiStack, 1240/E_eV, angle=Exp_Angles)

# Convert Ellipsometry data to complex dielectric function...
# Define functions for the quick conversion of data...
i = 2

rho = lambda psi, delta: np.tan(psi) * np.exp(1j * delta)
eps = lambda r, theta: np.sin(theta)**2 * (1 + np.tan(theta)**2 * ((1 - r)/(1 + r))**2)

# Experimental data...
Exp_rho = rho(np.radians(Exp_Data.data[Exp_Angles[i]][1]), np.radians((Exp_Data.data[Exp_Angles[i]][3])))
Exp_eps = eps(Exp_rho, np.radians(Exp_Angles[i]))

# Modelled data...
Mod_rho = rho(np.radians(Out["psi"][:,i]), np.radians(Out["Delta"][:,i]))
Mod_eps = eps(Mod_rho, np.radians(Exp_Angles[i]))

## Step 3 :: Data Plotting...
ax1b = ax1.twinx()
ax1.set_xlim([400, 1500])

# Experimental data...
ax1.plot(Exp_Data.data[Exp_Angles[i]][0]*1000, Exp_eps.real, lw=2, marker="o", ls='none', color=colours("Orange Red"),
         label="$\epsilon_1 (\lambda)$ :: $ %3.1f^{\circ}$" % Exp_Angles[i])
ax1b.plot(Exp_Data.data[Exp_Angles[i]][0]*1000, abs(Exp_eps.imag), lw=2, marker="s", ls='none', color=colours("Dodger Blue"),
         label="$\epsilon_2 (\lambda)$ :: $ %3.1f^{\circ}$" % Exp_Angles[i])

# Modelled Data...
ax1.plot(1240/E_eV, Mod_eps.real, label="Model $\epsilon_1 (\lambda)$ :: $ %3.1f^{\circ}$" % Exp_Angles[i], color=colours("Maroon"))
ax1b.plot(1240/E_eV, abs(Mod_eps.imag), label="Model $\epsilon_2 (\lambda)$ :: $ %3.1f^{\circ}$" % Exp_Angles[i], color=colours("Navy"))

ax1.set_xlabel("Wavelength (nm)")
ax1.set_ylabel('$\epsilon_1 (\lambda)$')
ax1b.set_ylabel('$\epsilon_2 (\lambda)$')
ax1.text(0.05, 0.9, '(b)', transform=ax1.transAxes, fontsize=12)

ax1.legend(loc="lower left")
ax1b.legend(loc="upper right")

plt.show()