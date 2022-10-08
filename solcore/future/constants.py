import numpy as np

from pint import Quantity


pi = np.pi
pi_ = Quantity(np.pi, "dimensionless")
"""Ehmmm... the Pi number"""

q = 1.60217646e-19
q_ = Quantity(q, "C")
"""Electron charge"""

kb = 1.3806503e-23
kb_ = Quantity(kb, "J/K")
"""Boltzmann's constant"""

h = 6.626068e-34
h_ = Quantity(h, "J*s")
"""Plank's constant"""

hbar = h / (2 * pi)
hbar_ = h_ / (2 * pi_)
"""Reduced Plank's constant"""

electron_mass = 9.10938188e-31
electron_mass_ = Quantity(electron_mass, "kg")
"""Electron rest mass"""

vacuum_permittivity = 8.854187817e-12
vacuum_permittivity_ = Quantity(vacuum_permittivity, "F/m")
"""Permitivity of the vacuum"""

c = 299792458.0
c_ = Quantity(c, "m/s")
"""Speed of light"""

fs = 6.8e-5
fs_ = Quantity(fs, "steradian")
"""Solid angle of the Sun - or entendu"""

Ts = 5762.0
Ts_ = Quantity(Ts, "k")
"""Temperature of the Sun when considered as a black body"""

solar_constant = 1361.0
solar_constant_ = Quantity(solar_constant, "kW/m**2")
"""Solar constant - density of solar irradiance out of the atmosphere"""
