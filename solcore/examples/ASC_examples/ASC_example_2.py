import matplotlib.pyplot as plt
import numpy as np
import solcore.analytic_solar_cells as ASC
from solcore import si, asUnit, eVnm

# We import the solar cell structure we created in the ASC_example_1
from ASC_example_1 import solar_cell

# Now we define the energies in which we want to calculate the QE.
energy_bounds = (0.6, 3.4, 1000)
E = si(np.linspace(*energy_bounds), "eV")

# We add this energy to the solar cell structure.
solar_cell.energies = E

# We run the simulation
qe_result = ASC.spectral_response_all_junctions(solar_cell, verbose=True)

# And here we have the QE of all of our junctions. Check the documentation to find out all the information provided as
# output, which is a lot.
qe_top = qe_result["junctions"][0]["qe_tot"]
qe_mid = qe_result["junctions"][1]["qe_tot"]
qe_bot = qe_result["junctions"][2]["qe_tot"]

wl = eVnm(asUnit(E, 'eV'))[::-1]
plt.plot(wl, qe_top[::-1], label='GaInP')
plt.plot(wl, qe_mid[::-1], label='InGaAs')
plt.plot(wl, qe_bot[::-1], label='Ge')
plt.xlim((300, 2000))
plt.xlabel('Wavelength (nm)')
plt.ylabel('Quantum Efficiency (-)')
plt.legend()
plt.show()
