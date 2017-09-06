"""
Example Script :: Illustrates how to import data from the SOPRA material constant database.
"""
import numpy as np
import matplotlib.pyplot as plt
from solcore.absorption_calculator import sopra_database

# Import material constant data for Gallium Arsenide :: Do this by placing the material name as the sole argument...
SOPRA_Material = sopra_database("GaAs")

# The materials_list() method opens the complete list of available materials in the SOPRA databse...
# SOPRA_Material.material_list()

# Example 1 :: Load n and k.txt data...
GaAs_n = SOPRA_Material.load_n()
GaAs_k = SOPRA_Material.load_k()

# Can also load alpha data...
GaAs_alpha = SOPRA_Material.load_alpha()

print(" GaAs Absorption coefficient")
print(GaAs_alpha)

# Example 2 :: Load extinction coefficient data for varying AlGaAs compositions...
SOPRA_Material = sopra_database("AlGaAs")

AlGaAs_k = []
# Define a fine step wavelength range to interpolatet data...
Wav = np.linspace(250, 800, 1000)
Al_frac = np.linspace(10, 100, 10)

# For a range of alloy compositions between 10 and 100%...
for comp in Al_frac:

    AlGaAs_k.append((SOPRA_Material.load_composition(Lambda=Wav, Al=comp)))

# plot example n and k.txt data for Gallium Arsenide...
fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(12, 10), dpi=100, facecolor='w')
fig.subplots_adjust(wspace=0.3)

# Subplot I... GaAs n and k.txt
AX1 = axs[0]
AX1_2 = AX1.twinx()

AX1.plot(GaAs_n[0], GaAs_n[1], linewidth=2)
AX1_2.plot(GaAs_k[0], GaAs_k[1], linewidth=2, color='red')

AX1.set_xlabel("Wavelength (nm)", fontweight="bold", fontstyle="italic")
AX1.set_ylabel("Refractive Index, n", fontweight="bold", fontstyle="italic", color='blue')
AX1.set_title("Example 1 :: n and k.txt data for GaAs", fontweight="bold", fontstyle="italic")

AX1_2.set_ylabel("Extinction Coefficient, k.txt", fontweight="bold", fontstyle="italic", color='red')

for tl in AX1.get_yticklabels():
    tl.set_color('blue')

for tl in AX1_2.get_yticklabels():
    tl.set_color('red')

# Subplot II... Varying AlGaAs compositions with increasing Al content...
AX2 = axs[1]
i = 0

for Wav, n, k in AlGaAs_k:

    AX2.plot(Wav, k, linewidth=2, label=("Al = %4.1f %%" % Al_frac[i]))
    i += 1

AX2.set_title("Eaxmple 2 :: AlGaAs k.txt data for varying Al fraction", fontweight="bold", fontstyle="italic")
AX2.set_xlabel("Wavelength (nm)", fontweight="bold", fontstyle="italic")
AX2.set_ylabel("Extinction Coefficient, k.txt", fontweight="bold", fontstyle="italic")
AX2.legend(loc="upper right")
plt.show()


