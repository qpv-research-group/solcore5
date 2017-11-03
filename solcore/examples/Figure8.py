import numpy as np
import matplotlib.pyplot as plt

from solcore.light_source import LightSource

plt.figure(figsize=(6, 4.5))

# The wavelength range of the spectra
wl = np.linspace(300, 3000, 200)

gauss = LightSource(source_type='laser', x=wl, center=800, linewidth=50, power=200)
bb = LightSource(source_type='black body', x=wl, T=5800, entendue='Sun')
am15g = LightSource(source_type='standard', x=wl, version='AM1.5g')
smarts = LightSource(source_type='SMARTS', x=wl)
spectral = LightSource(source_type='SPECTRAL2', x=wl)

plt.plot(*gauss.spectrum(), label='Gauss')
plt.plot(*bb.spectrum(), label='Black body')
plt.plot(*am15g.spectrum(), label='AM1.5G')
plt.plot(*smarts.spectrum(), label='SMARTS')
plt.plot(*spectral.spectrum(), label='SPECTRAL2')

plt.xlim(300, 3000)
plt.xlabel('Wavelength (nm)')
plt.ylabel('Power density (Wm$^{-2}$nm$^{-1}$)')
plt.tight_layout()
plt.legend(frameon=False)

plt.show()
