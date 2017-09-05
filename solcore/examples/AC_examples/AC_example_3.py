import matplotlib.pyplot as plt
import numpy as np
from solcore.absorption_calculator import calculate_ellipsometry
from solcore import material, si
from solcore.structure import Structure, Layer

# First we defined a couple of materials, for example, GaAs and AlGAs
GaAs = material('GaAs')(T=300)
AlGaAs = material('AlGaAs')(T=300, Al=0.3)

# Now, let's build the structure. We don't add a substrate and assume that there is an infinitely thick absorbing
# material in the back.
my_structure = Structure([
    Layer(si(40, 'nm'), material=AlGaAs),
    Layer(si(3000, 'nm'), material=GaAs),
])

# We want to calculate the ellipsometry of this structure as a function of the wavelength for a few angles
wavelength = np.linspace(400, 1000, 200)
angles = [60, 65, 70]

out = calculate_ellipsometry(my_structure, wavelength, angle=angles)

# This results must be taken with care. As the ellipsometry routine can only deal with coherent light, there might be
# strange oscillations related with intereference in thick layers, something that should not happen.
# Setting no_back_reflexion=True (the default) should take care of most of this effects, but might add other unexpected
# ones.
plt.plot(wavelength, out['psi'][:, 0], 'b', label=r'$\Psi$')
plt.plot(wavelength, out['Delta'][:, 0], 'r', label=r'$\Delta$')
for i in range(1, len(angles)):
    plt.plot(wavelength, out['psi'][:, i], 'b')
    plt.plot(wavelength, out['Delta'][:, i], 'r')

plt.xlabel('Wavelength (µm)')
plt.ylabel(r'$\Psi$ and $\Delta$ (º)')
plt.legend()
plt.show()
