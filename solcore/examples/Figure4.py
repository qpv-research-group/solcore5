from solcore import material
from solcore.quantum_mechanics import kp_bands
import numpy as np

# Material parameters
GaAs = material("GaAs")(T=300)
InGaAs = material("InGaAs")

comp = np.linspace(0.01, 0.25, 5)
for i in comp:
    InGaAs2 = InGaAs(In=i, T=300)

    result = kp_bands(GaAs, InGaAs2, graph=True, fit_effective_mass=True, effective_mass_direction="L")
