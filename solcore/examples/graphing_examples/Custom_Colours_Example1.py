# Custom_Colours Example 1 :: Example script showing how to make use of the Custom_Colours module within Solcore for
# adding some additional colour to your plots.

import numpy as np
import matplotlib.pyplot as plt
from solcore.graphing.Custom_Colours import *

# For this example, we will generate a number of Gaussian curves to illustrate the use of the 'colours' and
# 'colour_cycle' functions within Custom_Colours...

# First define a range of energies over which to calculate the Gaussian...
E = np.linspace(0, 3.5, 1000)
# E0_range defines a sequence of changing E0 to geneate the multiple data plots...
E0_range = np.linspace(0.75, 2.75, 23)

# Define the Gaussian function using the inline lambda statement...
Gauss = lambda E, A, E0, Sigma: A*np.exp( -1*(E - E0)**2 / (2*Sigma**2))

# Set up the plot...
fig = plt.figure(figsize=(10, 8), dpi=100, facecolor="w", edgecolor="w")

# To set a custom colour cycle, use the gca() method to assign axes properties to a variable.
AX = plt.gca()

# Set the custom colour cycle, by calling the 'colour_cycle' function with one of the defined sequences.
AX.set_prop_cycle(colour_cycle("origin_system_color_list"))

# Here, we plot to individual Gaussian instances using the 'colours' function to define their colours.
# They are the only curves plotted with a dashed line...
plt.plot(E, Gauss(E, 1, 0.5, 0.1), color=colours("Dark Green", type="HEX"), lw=3, ls="--")
plt.plot(E, Gauss(E, 1, 3, 0.1), color=colours("Dark Blue", type="RGB"), lw=3, ls="--")

# Plotting a succession of curves using the values in E0_range...
for x in E0_range:

    plt.plot(E, Gauss(E, 1, x, 0.1), lw=2)

plt.title("Custom_Colours Example", fontweight="bold", fontstyle="italic", fontsize=18)
plt.xlabel("Energy (eV)", fontweight="bold", fontstyle="italic", fontsize=18)
plt.ylabel("Amplitude (a.u.)", fontweight="bold", fontstyle="italic", fontsize=18)

plt.show()