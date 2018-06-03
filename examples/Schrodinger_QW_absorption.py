import matplotlib.pyplot as plt
import numpy as np

from solcore import si, material
from solcore.structure import Layer, Structure
import solcore.quantum_mechanics as QM
from solcore.constants import vacuum_permittivity, q

# First we create the materials we need
bulk = material("GaAs")(T=293, strained=False)
barrier = material("GaAsP")(T=293, P=0.1, strained=True)

# As well as some of the layers
top_layer = Layer(width=si("30nm"), material=bulk)
inter = Layer(width=si("3nm"), material=bulk)
barrier_layer = Layer(width=si("15nm"), material=barrier)
bottom_layer = top_layer

# We are going to calculate the absorption coefficient of InGaAs QWs of fixed thickness but different compositions
num_comp = 5
comp = np.linspace(0.05, 0.25, num_comp)
colors = plt.cm.jet(np.linspace(0, 1, num_comp))

# The absorption coefficients will be calculated at these energies and stored in alfas
num_energy = 300
E = np.linspace(1.15, 1.5, num_energy) * q

# We define some parameters need to calculate the shape of the excitonic absorption
alpha_params = {
    "well_width": si("7.2nm"),
    "theta": 0,
    "eps": 12.9 * vacuum_permittivity,
    "espace": E,
    "hwhm": si("6meV"),
    "dimensionality": 0.16,
    "line_shape": "Gauss"
}

# plt.figure(figsize=(4, 4.5))
for j, i in enumerate(comp):
    # We create the QW material at the given composition
    QW = material("InGaAs")(T=293, In=i, strained=True)

    # And the layer
    well_layer = Layer(width=si("7.2nm"), material=QW)

    # The following lines create the QW structure, with different number of QWs and interlayers
    test_structure = Structure([barrier_layer, inter, well_layer, inter, barrier_layer], substrate=bulk)

    # Finally, the quantum properties are claculated here
    output = QM.schrodinger(test_structure, quasiconfined=0, num_eigenvalues=20, alpha_params=alpha_params,
                            calculate_absorption=True)

    alfa = output[0]['alphaE'](E)
    plt.plot(1240 / (E / q), alfa / 100, label='{}%'.format(int(i * 100)))

plt.xlim(826, 1100)
plt.ylim(0, 23000)
plt.xlabel('Wavelength (nm)')
plt.ylabel('$\\alpha$ cm$^{-1}$')
plt.legend(loc='upper right', frameon=False)
plt.tight_layout()

plt.show()
