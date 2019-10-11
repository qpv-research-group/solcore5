from solcore import material
from solcore import si
from solcore.material_system import create_new_material
from solcore.absorption_calculator import create_nk_txt, download_db, search_db
from solcore.config_tools import add_source
import matplotlib.pyplot as plt
import numpy as np
import os

# TODO Add detailed instrucctions on how to setup the environment to run this example
#  as it is now, it fails and is unclear how to make it work.
# first we need to tell Solcore some things about where to put custom materials. for this,
# we use the add_source function from config_tools, although we could also manually edit
# the solcore configuration file (which should be in your home directory).
# You need to add two things to the config file: where to put the n and k data for new
# materials added to the database, and where to put the other parameters (these can all
# go in the same file).

home_folder = os.path.expanduser('~')
custom_nk_path = os.path.join(home_folder, 'Solcore/custommats')
nk_db_path = os.path.join(home_folder, 'Solcore/NK.db')
param_path = os.path.join(home_folder, 'Solcore/custom_params.txt')

add_source('Others', 'custom_mats', custom_nk_path)
add_source('Others', 'nk', nk_db_path)
add_source('Parameters', 'custom', param_path)

# EXAMPLE 1

# need to have n and k data, and a parameter file in the correct format -
# see examples of parameter files in e.g. material_data/Adachi/binaries.txt

# create a new material, silicon-germanium-tin, from input files. Here,
# the parameters in SiGeSn_params.txt have been copied directly from Ge.
create_new_material('SiGeSn', 'SiGeSn_n.txt', 'SiGeSn_k.txt', 'SiGeSn_params.txt')

# can now create an instance of it like any Solcore material
SiGeSn = material('SiGeSn')()

plt.figure()
plt.plot(si(np.arange(300, 1700, 5), 'nm')*1e9, SiGeSn.n(si(np.arange(300, 1700, 5), 'nm')))
plt.plot(si(np.arange(300, 1700, 5), 'nm')*1e9, SiGeSn.k(si(np.arange(300, 1700, 5), 'nm')))

plt.xlabel('Wavelength (nm)')
plt.ylabel('SiGeSn n / k')

plt.show()

# EXAMPLE 2
# Can also create a Solcore material from a material in the refractiveindex.info database:

# if necessary, download database:
download_db()

# search what options are available for diamond, then use the first result's pageid to
# create data files for the n & k of diamond:

results = search_db('Diamond')
create_nk_txt(pageid=results[0][0], file='C_Diamond')
create_new_material(mat_name = 'Diamond', n_source='C_Diamond_n.txt', k_source='C_Diamond_k.txt')

Diamond = material('Diamond')()

plt.figure()
plt.plot(si(np.arange(100, 800, 5), 'nm')*1e9, Diamond.n(si(np.arange(100, 800, 5), 'nm')))
plt.plot(si(np.arange(100, 800, 5), 'nm')*1e9, Diamond.k(si(np.arange(100, 800, 5), 'nm')))

plt.xlabel('Wavelength (nm)')
plt.ylabel('Diamond n / k')

plt.show()
