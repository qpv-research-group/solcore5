import matplotlib.pyplot as plt
import numpy as np
from solcore import material

from solcore.solar_cell import SolarCell
from solcore.structure import Layer
from solcore.absorption_calculator import calculate_rat, OptiStack

wavelengths = np.linspace(310, 1900, 400)
AlInP = material("AlInP")(Al=0.5)
GaInP = material("GaInP")(In=0.5)
MgF = material("MgF2")()
SiC = material("Si3N4")()
ZnS = material("ZnScub")()

plt.figure()
plt.plot(wavelengths, MgF.n(wavelengths*1e-9))
plt.plot(wavelengths, SiC.n(wavelengths*1e-9))
plt.plot(wavelengths, ZnS.n(wavelengths*1e-9))
plt.plot(wavelengths, AlInP.n(wavelengths*1e-9))
plt.plot(wavelengths, GaInP.n(wavelengths*1e-9))
plt.legend(['MgF', "SiC", "ZnS", "AlInP", "Ga"])
plt.show()

plt.figure()
plt.plot(wavelengths, MgF.k(wavelengths*1e-9))
plt.plot(wavelengths, SiC.k(wavelengths*1e-9))
plt.plot(wavelengths, ZnS.k(wavelengths*1e-9))
plt.plot(wavelengths, AlInP.k(wavelengths*1e-9))
plt.plot(wavelengths, GaInP.k(wavelengths*1e-9))
plt.legend(['MgF', "SiC", "ZnS", "AlInP", "Ga"])
plt.show()


# Build the optical stack...
stack = SolarCell([Layer(117e-9, MgF),
                   Layer(80e-9, SiC),
                   Layer(61e-9, ZnS),
                   Layer(25e-9, AlInP)], substrate=GaInP)

angles = np.linspace(0, 80, 10)
RAT_angles = []

print("Calculate RAT ::")
for theta in angles:
    print("Calculating at angle :: %4.1f deg" % theta)
    # Calculate RAT data...
    rat_data = calculate_rat(stack, angle=theta, wavelength=wavelengths, no_back_reflection=False)

    RAT_angles.append((theta, rat_data["R"], rat_data["A"], rat_data["T"]))

colors = plt.cm.jet(np.linspace(1, 0, len(RAT_angles)))

fig, ax2 = plt.subplots(1, 1)

for i, RAT in enumerate(RAT_angles):
    ax2.plot(wavelengths, RAT[1] * 100, ls="-", color=colors[i], label="%4.1f$^\circ$" % RAT[0])
    ax2.plot(wavelengths, (RAT[2] + RAT[3]) * 100, ls="--", color=colors[i])

ax2.set_ylim([0, 100])
ax2.set_xlim([300, 1800])
ax2.set_xlabel("Wavelength (nm)")
ax2.set_ylabel("Reflection and transmission into cell (%)")
ax2.legend(loc=5)
ax2.text(0.05, 0.45, '(a)', transform=ax2.transAxes, fontsize=12)

plt.tight_layout(w_pad=4)
plt.show()


fig, ax2 = plt.subplots(1, 1)

for i, RAT in enumerate(RAT_angles):
    ax2.plot(wavelengths, RAT[2] * 100, ls="-", color=colors[i], label="%4.1f$^\circ$" % RAT[0])
    ax2.plot(wavelengths, (RAT[3]) * 100, ls="--", color=colors[i])

ax2.set_ylim([0, 100])
ax2.set_xlim([300, 1800])
ax2.set_xlabel("Wavelength (nm)")
ax2.set_ylabel("Reflection and transmission into cell (%)")
ax2.legend(loc=5)
ax2.text(0.05, 0.45, '(a)', transform=ax2.transAxes, fontsize=12)

plt.tight_layout(w_pad=4)
plt.show()