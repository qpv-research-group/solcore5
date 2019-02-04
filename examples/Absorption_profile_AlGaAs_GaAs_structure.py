import matplotlib.pyplot as plt
import numpy as np
from solcore.absorption_calculator import calculate_absorption_profile
from solcore import material, si
from solcore.structure import Structure, Layer

# First we defined a couple of materials, for example, GaAs and AlGAs
GaAs = material('GaAs')(T=300)
AlGaAs = material('AlGaAs')(T=300, Al=0.3)

# Now, let's build the structure
my_structure = Structure([
    Layer(si(30, 'nm'), material=AlGaAs),
    Layer(si(3000, 'nm'), material=GaAs),
    Layer(si(300, 'um'), material=GaAs),
])

# We want to calculate the absorption profile of this structure as a function of the position and wavelength
wl = np.linspace(400, 1000, 200)
out = calculate_absorption_profile(my_structure, wl, steps_size=1, z_limit=3000)

# Finally, we plot the absorption profile. Note that absorption at short wavelengths take place near the surface of the
# structure, in the AlGaAs layer and top of the GaAs layer, while longer wavelengths penetrate more. Wavelengths beyond
# the GaAs band edge are not absorbed.
plt.figure(1)
ax = plt.contourf(out['position'], wl, out['absorption'], 200)
plt.xlabel('Position (nm)')
plt.ylabel('Wavelength (nm)')
cbar = plt.colorbar()
cbar.set_label('Absorption (1/nm)')

# We can also check what is the overall light absorption in the AlGaAs and GaAs epilayers as that represents the total
# light absorption in our solar cell and therefore the maximum EQE (light absorbed in the thick substrate lost).
# The absorption is mostly limited by the reflexion in the front surface. Clearly, this solar cell needs an AR coating.
A = np.zeros_like(wl)

for i, absorption in enumerate(out['absorption'][:]):
    A[i] = np.trapz(absorption, out['position'])

plt.figure(2)
plt.plot(wl, A)
plt.xlabel('Wavelength (nm)')
plt.ylabel('Absorption')

plt.show()
