import os

SOLCORE_ROOT = os.path.split(__file__)[0]
default_config = os.path.join(SOLCORE_ROOT, "solcore_config.txt")
user_config = os.path.join(os.path.expanduser("~"), ".solcore_config.txt")

from .config_tools import SolcoreConfig

config = SolcoreConfig(default_config, user_config)
verbose = config.verbose_loading()

if config.welcome_message():
    print(
        f"""\n\tWelcome to Solcore - version {config.version()}
\tCopyright (c) 2019, Imperial College London. All rights reserved.
\tSoftware released under the GNU Lesser General Public License.\n"""
    )

from solcore.units_system import UnitsSystem

# First we populate the Units system:
UnitsSystem(config["Units"])

# And now we load some functions form it.
si = UnitsSystem().si
asUnit = UnitsSystem().asUnit
siUnits = UnitsSystem().siUnits
sensibleUnits = UnitsSystem().sensibleUnits
siUnitFromString = UnitsSystem().siUnitFromString
convert = UnitsSystem().convert
guess_dimension = UnitsSystem().guess_dimension
nmJ = UnitsSystem().nmJ
mJ = UnitsSystem().mJ
eVnm = UnitsSystem().eVnm
nmHz = UnitsSystem().nmHz
spectral_conversion_nm_ev = UnitsSystem().spectral_conversion_nm_ev
spectral_conversion_nm_hz = UnitsSystem().spectral_conversion_nm_hz
eV = UnitsSystem().eV

# And the same with the Parameter system
from solcore.parameter_system import ParameterSystem

ParameterSystem(config["Parameters"])
get_parameter = ParameterSystem().get_parameter


# And the same with the Materials system
from solcore.material_system import MaterialSystem

MaterialSystem(config["Materials"])
material = MaterialSystem().material
