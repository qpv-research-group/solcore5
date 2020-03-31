import os

SOLCORE_ROOT = os.path.split(__file__)[0]
default_config = os.path.join(SOLCORE_ROOT, "solcore_config.txt")
user_path = os.environ.get("SOLCORE_USER_DATA")
if user_path is None:
    user_path = os.path.join(os.path.expanduser("~"), ".solcore")
if not os.path.exists(user_path):
    os.mkdir(user_path)
user_config = os.path.join(user_path, "solcore_config.txt")

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
units_system = UnitsSystem(config.units)
config.register_observer("Units", units_system.read)

# And now we load some functions form it.
si = units_system.si
asUnit = units_system.asUnit
siUnits = units_system.siUnits
sensibleUnits = units_system.sensibleUnits
siUnitFromString = units_system.siUnitFromString
convert = units_system.convert
guess_dimension = units_system.guess_dimension
nmJ = units_system.nmJ
mJ = units_system.mJ
eVnm = units_system.eVnm
nmHz = units_system.nmHz
spectral_conversion_nm_ev = units_system.spectral_conversion_nm_ev
spectral_conversion_nm_hz = units_system.spectral_conversion_nm_hz
eV = units_system.eV

# And the same with the Parameter system
from solcore.parameter_system import ParameterSystem

parameters_system = ParameterSystem(config.parameters)
config.register_observer("Parameters", parameters_system.read)
get_parameter = parameters_system.get_parameter

# And the same with the Materials system
from solcore.material_system import MaterialSystem

materials_system = MaterialSystem(config.materials)
config.register_observer("Materials", materials_system.read)
material = materials_system.material
