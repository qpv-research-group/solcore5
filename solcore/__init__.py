import os
from datetime import datetime

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
__version__ = config.version()

if config.welcome_message():
    print(
        f"""\n\tWelcome to Solcore - version {config.version()}
\tCopyright (c) {datetime.today().year}, Imperial College London. All rights reserved.
\tSoftware released under the GNU Lesser General Public License.\n"""
    )

from solcore.units_system import UnitsSystem

# First we populate the Units system:
us = UnitsSystem(config.units)
config.register_observer("Units", us.read)

# And now we load some functions form it.
si = us.si
asUnit = us.asUnit
siUnits = us.siUnits
sensibleUnits = us.sensibleUnits
siUnitFromString = us.siUnitFromString
convert = us.convert
guess_dimension = us.guess_dimension
nmJ = us.nmJ
mJ = us.mJ
eVnm = us.eVnm
nmHz = us.nmHz
spectral_conversion_nm_ev = us.spectral_conversion_nm_ev
spectral_conversion_nm_hz = us.spectral_conversion_nm_hz
eV = us.eV

# And the same with the Parameter system
from solcore.parameter_system import ParameterSystem

ps = ParameterSystem(config.parameters)
config.register_observer("Parameters", ps.read)
get_parameter = ps.get_parameter

# And the same with the Materials system
from solcore.material_system import MaterialSystem

ms = MaterialSystem(config.materials)
config.register_observer("Materials", ms.read)
material = ms.material
