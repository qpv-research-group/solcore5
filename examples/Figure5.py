from solcore import si, material
from solcore.structure import Layer, Structure
import matplotlib.pyplot as plt
import numpy as np
import solcore.quantum_mechanics as QM
from solcore.constants import vacuum_permittivity, q


# First we create the materials we need
bulk = material("GaAs")(T=293)
barrier = material("GaAsP")(T=293, P=0.1)

bulk.strained = False
barrier.strained = True

# As well as some of the layers
top_layer = Layer(width=si("15nm"), material=barrier)
inter = Layer(width=si("3nm"), material=bulk)
barrier_layer = Layer(width=si("5nm"), material=barrier)
bottom_layer = top_layer

E = np.linspace(1.15, 1.5, 300) * q
alfas = np.zeros((len(E), 6))

alfas[:, 0] = E / q

# We define some parameters need to calculate the shape of the excitonic absorption
alpha_params = {
    "well_width": si("10nm"),
    "theta": 0,
    "eps": 12.9 * vacuum_permittivity,
    "espace": E,
    "hwhm": si("6meV"),
    "dimensionality": 0.16,
    "line_shape": "Gauss"
}

# We create the QW material at the given composition
QW = material("InGaAs")(T=293, In=0.15)
QW.strained = True

# And the layer
well_layer = Layer(width=si("7.2nm"), material=QW)

# The following lines create the QW structure, with different number of QWs and interlayers
# test_structure = Structure([top_layer, barrier_layer, inter] + 10 * [well_layer, inter, barrier_layer, inter] +
#                            [bottom_layer])

test_structure = Structure([top_layer, inter] + 1 * [well_layer, inter] +
                            [bottom_layer])

#test_structure = Structure([top_layer, barrier_layer] + 10 * [well_layer, barrier_layer] +
#                           [bottom_layer])

test_structure.substrate = bulk

# Finally, the quantum properties are claculated here
output = QM.schrodinger(test_structure, quasiconfined=0.1, graphtype='potentials', periodic=False,
                     num_eigenvalues=200, alpha_params=alpha_params, calculate_absorption=False, Efield=0)

