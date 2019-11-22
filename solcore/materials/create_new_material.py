# to add a new material to the database, so that it can be created like a built-in
# Solcore material, its path has to be added to the config file and MANIFEST.in
# Have to copy the relevant n, k and parameter files into that folder.

import os
from re import sub
from shutil import copyfile

from solcore import SOLCORE_ROOT
from .parameter_system import ParameterSystem


def create_new_material(mat_name, n_source, k_source, parameter_source=None):
    """This function adds a new material to Solcore's material_data folder, so that it
    can be called like a built-in material. It needs a name for the new material, and
    source files for the n and k data and other parameters which will be copied into the
    material_data/Custom folder.

    :param mat_name: the name of the new material
    :param n_source: path of the n values
    (txt file, first column wavelength in m, second column n)
    :param k_source: path of
    the n values (txt file, first column wavelength in m, second column k)

    :return:
    parameter_source: file with list of materials for the new material
    """
    from solcore.config_tools import add_source, user_config_data as config

    CUSTOM_PATH = os.path.abspath(
        config["Others"]["custom_mats"].replace("SOLCORE_ROOT", SOLCORE_ROOT)
    )
    PARAMETER_PATH = os.path.abspath(
        config["Parameters"]["custom"].replace("SOLCORE_ROOT", SOLCORE_ROOT)
    )

    # check if there is already a material with this name
    if mat_name not in sorted(ParameterSystem().database.sections()):

        # create a folder in the custom materials folders
        folder = os.path.join(CUSTOM_PATH, mat_name + "-Material")
        if not os.path.exists(folder) and folder != "":
            os.makedirs(folder)
        else:
            print("This material already exists (or at least, a folder for it).")

        # copy n and k data files to the material's folder
        copyfile(n_source, os.path.join(folder, "n.txt"))
        copyfile(k_source, os.path.join(folder, "k.txt"))

        # create the parameter file if it doesn't already exist
        if not os.path.isfile(PARAMETER_PATH):
            open(PARAMETER_PATH, "a").close()

        # append the parameters for the new material
        fout = open(PARAMETER_PATH, "r")
        existing_parameters = fout.read()
        fout.close()

        if not "[" + mat_name + "]" in existing_parameters:
            # make sure the names match
            if parameter_source is not None:
                fin = open(parameter_source, "r")
                parameters = fin.read() + "\n\n"
                parameters = sub(
                    "\[[^]]*\]",
                    lambda x: x.group(0).replace(x.group(0), "[" + mat_name + "]"),
                    parameters,
                )
                fin.close()
            else:
                parameters = "[" + mat_name + "]\n\n"
                print(
                    "Material created with optical constants n and k only, no other "
                    "parameters provided."
                )

            fout = open(PARAMETER_PATH, "a")
            fout.write(parameters)
            fout.close()
        else:
            print(
                "There are already parameters for this material in the custom "
                "parameter file at " + PARAMETER_PATH
            )

        # modify the user's config file (in their home folder) to include the
        # relevant paths
        new_entry = (
            mat_name
            + " = "
            + config["Others"]["custom_mats"]
            + "/"
            + mat_name
            + "-Material\n"
        )
        if not new_entry in config:

            add_source(
                "Materials",
                mat_name,
                config["Others"]["custom_mats"] + "/" + mat_name + "-Material",
            )

        else:
            print(
                "A path for this material was already added to the Solcore config file "
            )

    else:
        print("There is already a material with this name - choose a different one.")
