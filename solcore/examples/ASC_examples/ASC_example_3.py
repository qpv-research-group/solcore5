import numpy as np
import solcore.analytic_solar_cells as ASC
from solcore import si
from solcore.light_source import calculate_spectrum_spectral2
from solcore.interpolate import interp1d

# We import the solar cell structure we created in the ASC_example_1
from ASC_example_1 import solar_cell

# Configure cell geometry and concentration
# Cell area is stated to be 7mm x 7mm
cell_area = 0.7 * 0.7 / 1e4  # expressed in m-2
# Geometrical concentration is 820, but the optical efficiency is approximately 0.85 and cell metalisation shading 8%
concentration = 1000
concentration_factor = concentration * 0.85 * (1 - 0.08)

# We create a solar spectrum using SPECTRAL2 and the default configuration
# - See documentation and example of that package for more info
spectrum = calculate_spectrum_spectral2()

# We use the spectrum is in SI units: Watts m-2 joule-1.
incident_x_J_y_per_J = spectrum["incident spectrum energy si"]
incident_function = interp1d(y=incident_x_J_y_per_J[1] * concentration_factor * cell_area, x=incident_x_J_y_per_J[0])
power_density = spectrum["incident power density"]
print('Power density = {:.2f} W m-2\n'.format(power_density))

# Now we define the energies in which we want to calculate the QE.
energy_bounds = (0.6, 3.4, 1000)
E = si(np.linspace(*energy_bounds), "eV")       # In SI units
solar_cell.incident_light = E, incident_function(E) / E
qe_result = ASC.spectral_response_all_junctions(solar_cell, verbose=False)

print("Subcell photocurrent density [mA cm-2]:")
print('\tJsc(Top) = {:.2f}'.format(1e3 * qe_result["junctions"][0]["J"] / (cell_area * 1e4 * concentration_factor)))
print('\tJsc(Mid) = {:.2f}'.format(1e3 * qe_result["junctions"][1]["J"] / (cell_area * 1e4 * concentration_factor)))
print('\tJsc(Bot) = {:.2f}'.format(1e3 * qe_result["junctions"][2]["J"] / (cell_area * 1e4 * concentration_factor)))