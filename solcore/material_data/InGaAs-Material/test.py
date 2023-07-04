import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

cols = sns.color_palette("husl", 7)
fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(8, 4))
for i1, comp in enumerate(["000", "010", "024", "035", "053", "075", "100"]):

    nk = np.loadtxt(comp + "_InGaAs.txt")
    ax1.plot(nk[:,0], nk[:,1], color=cols[i1], label=comp)
    ax2.plot(nk[:,0], nk[:,2], color=cols[i1])

ax1.legend()
ax1.set_xlabel("E (eV)")
ax2.set_xlabel("E (eV)")
ax1.set_ylabel("e1")
ax2.set_ylabel("e2")
plt.show()