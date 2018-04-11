# Package file for the Poisson - Drift-Diffusion solver

# We import all the contents of the modules of this package so they can be used directly
from .DeviceStructure import *
from .QWunit import *

try:
    #result = check_ddModel_library_ok()
    #if result:
    from .DriftDiffusionUtilities import *
except Exception as err:
    print('WARNING: The Poisson - Drift-Diffusion solver will not be available because the ddModel fortran library could not be imported.')

