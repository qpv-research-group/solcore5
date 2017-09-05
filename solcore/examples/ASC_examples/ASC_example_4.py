import numpy as np
from solcore.structure import Structure, Junction
import solcore.analytic_solar_cells as ASC
from solcore.state import State
from solcore.graphing import *

# We import the QE results obtained in ASC_example_3 which contain the short circuit currents of each junction
# when the solar cell is illuminated with the default SPECTAL2 spectrum
from ASC_example_3 import qe_result, power_density, cell_area, concentration_factor

# Ref temperature ºC. Let's say we have the reverse saturation currents J01 and J02 of each junction at a
# reference temperature and we want the IV curve at a different one. We have to make a correction to the J01 and J02
cell_temp = 60
ref_temp = 25

Tcell = 273 + cell_temp
Tref = 273 + ref_temp

# The IV data will be stored in a State object. We create it, including the cell and reference temperatures.
IV_calculation_state = State(T=Tcell, Tref=Tref)

# From the QE object we get the short circuit currents
Isc_array = [qe_result["junctions"][0]["J"], qe_result["junctions"][1]["J"], qe_result["junctions"][2]["J"]]

# And we create a list with the reverse saturation currents. In this case, we don't calculate them but just assume we
# have them from somewhere.
I01_array = [4.93e-24, 1.0e-21, 4.93e-6]
I02_array = [3.28e-15, 2.7e-10, 1.0e-5]

# This is the structure to calculate.
IV_calculation_state.structure = Structure(
    [
        Junction(Eg=1.9, j01=I01_array[0], j02=I02_array[0], R_shunt=3e6,
                 R_series=0.0236, n1=1.00, n2=2.0, photocurrent=Isc_array[0]),
        Junction(Eg=1.4, j01=I01_array[1], j02=I02_array[1], R_shunt=1.5e6,
                 R_series=0.0012, n1=1.00, n2=2.0, photocurrent=Isc_array[1]),
        Junction(Eg=0.66, j01=I01_array[2], j02=I02_array[2], R_shunt=115,
                 R_series=8e-4, n1=1.00, n2=2.0, photocurrent=Isc_array[2]),
    ]
)

# We solve it, including explicitely the range of voltages we are interested
IV_result = ASC.multijunctionIV(IV_calculation_state, V=np.linspace(0, 4, 1000))

# We use the tools of the graphing package to get a nice plot of the IV curves.
junction_colors = ["blue", "green", "red"]
graph_lines = [GraphData(iv, label="Junction {}".format(i + 1), color=junction_colors[i])
               for i, iv in enumerate(IV_result["junction IV"])]
graph_lines.append(GraphData(IV_result["IV"], linewidth=2, color="black", label="Multijunction"))
g = Graph(graph_lines, ylim=(0, 7), xlabel="Bias (V)", ylabel="Current (A)", legend="best").draw()

eta = IV_result["Pmpp"] / (power_density * cell_area * concentration_factor)

print('\nThe solar cell properties are: ')
print('\tIsc = {:.2f} mA cm-2'.format(1e3 * IV_result["Isc"] / (cell_area * 1e4 * concentration_factor) ))
print('\tVoc = {:.2f} V'.format(IV_result["Voc"]))
print('\tFF = {:.2f} % '.format(IV_result["FF"] * 100))
print('\tEta = {:.2f} %'.format(eta * 100))
