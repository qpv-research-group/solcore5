# Package file for the Poisson - Drift-Diffusion solver

# We import all the contents of the modules of this package so they can be used directly
from .DeviceStructure import *
from .Illumination import *
from .QWunit import *
from .driftdiffusion_compiler import check_ddModel_library_ok

try:
    result = check_ddModel_library_ok()
    if result:
        from .DriftDiffusionUtilities import *
except Exception as err:
    print('\nERROR: The Drift Diffusion Utilities could not be loaded because of an issue with the solver library:')
    print(err)

