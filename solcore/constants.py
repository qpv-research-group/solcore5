from pint import Quantity
import numpy as np

pi = Quantity(np.pi, "dimensionless")
"""Ehmmm... the Pi number"""

q = Quantity(1.60217646e-19, "C")
"""Electron charge"""

kb = Quantity(1.3806503e-23, "J/K")
"""Boltzmann's constant"""

h = Quantity(6.626068e-34, "J*s")
"""Plank's constant"""

hbar = h / (2 * pi)
"""Reduced Plank's constant"""

electron_mass = Quantity(9.10938188e-31, "kg")
"""Electron rest mass"""

vacuum_permittivity = Quantity(8.854187817e-12, "F/m")
"""Permitivity of the vacuum"""

c = Quantity(299792458.0, "m/s")
"""Speed of light"""

fs = Quantity(6.8e-5, "steradian")
"""Solid angle of the Sun - or entendu"""

Ts = Quantity(5762.0, "k")
"""Temperature of the Sun when considered as a black body"""

solar_constant = Quantity(1361.0, "kW/m**2")
"""Solar constant - density of solar irradiance out of the atmosphere"""
