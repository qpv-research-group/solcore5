import numpy as np
import matplotlib.pyplot as plt

from solcore.light_source import LightSource

# TODO If SMARTS is not installed, it keeps trying and producing errors.
#  LightSource is poorly designed...
# The wavelength range of the spectra
wl = np.linspace(300, 3000, 200)

# Now different types of light sources can be defined
gauss = LightSource(source_type='laser', x=wl, center=800, linewidth=50, power=200)
bb = LightSource(source_type='black body', x=wl, T=5800, entendue='Sun')
am15g = LightSource(source_type='standard', x=wl, version='AM1.5g')
spectral = LightSource(source_type='SPECTRAL2', x=wl)

# Plot comparing the different light sources
plt.figure(1)
plt.plot(*gauss.spectrum(), label='Gauss')
plt.plot(*bb.spectrum(), label='Black body')
plt.plot(*am15g.spectrum(), label='AM1.5G')
plt.plot(*spectral.spectrum(), label='SPECTRAL2')

try:
    smarts = LightSource(source_type='SMARTS', x=wl)
    plt.plot(*smarts.spectrum(), label='SMARTS')
except TypeError:
    pass

plt.xlim(300, 3000)
plt.xlabel('Wavelength (nm)')
plt.ylabel('Power density (Wm$^{-2}$nm$^{-1}$)')
plt.tight_layout()
plt.legend()

try:
    # Plot comparing the spectra calculated with SMARTS at different hours of the day
    for h in range(8, 20):
        plt.plot(*smarts.spectrum(HOUR=h), label='{} h'.format(h))

    plt.xlim(300, 3000)
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Power density (Wm$^{-2}$nm$^{-1}$)')
    plt.tight_layout()
    plt.legend()
except TypeError:
    pass

plt.show()