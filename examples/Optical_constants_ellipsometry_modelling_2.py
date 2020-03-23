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
angles = [65, 70, 75]

out = calculate_ellipsometry(my_structure, wavelength, angle=angles)

# This results must be taken with care. As the ellipsometry routine can only deal with coherent light, there might be
# strange oscillations related with intereference in thick layers, something that should not happen.
# Setting no_back_reflection=True (the default) should take care of most of this effects, but might add other unexpected
# ones.
fig, ax1 = plt.subplots(1, 1)
ax2 = ax1.twinx()

ax1.plot(wavelength, out['psi'][:, 0], 'b', label=r'$\Psi$')
ax2.plot(wavelength, out['Delta'][:, 0], 'r', label=r'$\Delta$')
for i in range(1, len(angles)):
    ax1.plot(wavelength, out['psi'][:, i], 'b')
    ax2.plot(wavelength, out['Delta'][:, i], 'r')

ax1.set_xlabel("Wavelength (nm)")
ax1.set_ylabel(r'$\Psi$ (º)')
ax2.set_ylabel(r'$\Delta$ (º)')

ax1.legend(loc="upper left")
ax2.legend(loc="upper right")

plt.show()
