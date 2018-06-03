""" Quick example code to generate some data from the SOPRA database"""

from solcore import material

import numpy as np
import matplotlib.pyplot as plt

# Initiate some variables...
wl = np.linspace(300, 1800, 1000)
Al_fraction = np.linspace(10, 100, 10)

f, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.25,4))
ax1b = ax1.twinx()

# Load GaAs n and k data from the Sopra database...
GaAs = material("GaAs", sopra=True)(T=300)
GaAs_n = GaAs.n(wl*1e-9)
GaAs_k = GaAs.k(wl*1e-9)

# Load Ge n and k data from Sopra database...
Ge = material("Ge", sopra=True)(T=300)
Ge_n = Ge.n(wl*1e-9)
Ge_k = Ge.k(wl*1e-9)

# Load AlGaAs k data for a range of compositions...
AlGaAs = material("AlGaAs", sopra=True)

AlGaAs_k = [GaAs_k]

for comp in Al_fraction:
    AlGaAs_k.append(AlGaAs(T=300, Al=comp/100).k(wl*1e-9))

# Plot the data...
colors = plt.cm.jet(np.linspace(0,1,len(Al_fraction)+1))

lns1 = ax1.plot(wl, GaAs_n, 'b', label='n, GaAs')
lns2 = ax1b.plot(wl, GaAs_k, 'r', label='k, GaAs')

lns3 = ax1.plot(wl, Ge_n, ls="--", color='blue', label='n, Ge')
lns4 = ax1b.plot(wl, Ge_k,ls="--", color='red', label='k, Ge')

ax1.set_xlim([300,1800])
ax1b.set_xlim([300,1800])
ax1b.set_ylim([0, 3.8])

# added these three lines
lns = lns1+lns2+lns3+lns4
labs = [l.get_label() for l in lns]
ax1.legend(lns, labs, loc="upper right", frameon=False)
ax1.text(0.05, 0.9, '(a)', transform=ax1.transAxes, fontsize=12)

ax1.set_xlabel("Wavelength (nm)")
ax1.set_ylabel("Refractive Index, n")
ax1b.set_ylabel("Extinction Coefficient, k")

for i, k in enumerate(Al_fraction):
    ax2.plot(wl, AlGaAs_k[i], color=colors[i+1], label='{}%'.format(int(Al_fraction[i])))

ax2.set_xlim([300, 900])
ax2.set_ylim([0, 2.8])

ax2.set_xlabel("Wavelength (nm)")
ax2.set_ylabel("Extinction Coefficient, k")
ax2.legend(loc="upper right", frameon=False)
ax2.text(0.05, 0.9, '(b)', transform=ax2.transAxes, fontsize=12)
plt.tight_layout(w_pad=4)

plt.show()

