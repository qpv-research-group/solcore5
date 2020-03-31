# to add a new material to the database, so that it can be created like a built-in
# Solcore material, its path has to be added to the config file and MANIFEST.in
# Have to copy the relevant n, k and parameter files into that folder.

import os

from shutil import copyfile, move
from re import sub
from solcore import config, SOLCORE_ROOT
from solcore.parameter_system import ParameterSystem
from configparser import ConfigParser
from solcore.material_system import MaterialSystem


def create_new_material(mat_name, n_source, k_source, parameter_source = None):
    """
    This function adds a new material to Solcore's material_data folder, so that it can be called like a
    built-in material. It needs a name for the new material, and source files for the n and k data and other
    parameters which will be copied into the material_data/Custom folder.

    :param mat_name: the name of the new material
    :param n_source: path of the n values (txt file, first column wavelength in m, second column n)
    :param k_source: path of the n values (txt file, first column wavelength in m, second column k)
    :param: parameter_source: file with list of parameters for the new material
    """

    PARAMETER_PATH = os.path.join(config.user_folder, "custom_parameters.txt")
    if "custom" not in config.parameters():
        if not os.path.isfile(PARAMETER_PATH):
            open(PARAMETER_PATH, 'a').close()
        config["Parameters", "custom"] = PARAMETER_PATH

    CUSTOM_PATH = os.path.join(config.user_folder, "custom_materials")

    # check if there is already a material with this name
    if mat_name in sorted(ParameterSystem().database.sections()) or mat_name in config.materials():
        answer = input(f"A material named {mat_name} already exists in the database."
                       f"Do you want to overwrite it [y/n]?")
        if answer.lower() != "y":
            return

    # create a folder in the custom materials folders
    folder = os.path.join(CUSTOM_PATH, mat_name + '-Material')
    if not os.path.exists(folder) and folder != "":
        os.makedirs(folder)

    # copy n and k data files to the material's folder
    copyfile(n_source, os.path.join(folder, 'n.txt'))
    copyfile(k_source, os.path.join(folder, 'k.txt'))

    config["Materials", mat_name] = folder

    # append the parameters for the new material
    params = ConfigParser()
    params.optionxform = str
    if parameter_source is not None:
        params.read([PARAMETER_PATH, parameter_source])
        with open(PARAMETER_PATH, "w") as fp:
            params.write(fp)
    else:
        params.read([PARAMETER_PATH])
        params[mat_name] = {}
        with open(PARAMETER_PATH, "w") as fp:
            params.write(fp)
        print('Material created with optical constants n and k only.')

    ParameterSystem().read()
