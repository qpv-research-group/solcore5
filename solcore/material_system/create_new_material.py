# to add a new material to the database, so that it can be created like a built-in
# Solcore material, its path has to be added to the config file and MANIFEST.in
# Have to copy the relevant n, k and parameter files into that folder.

# Presumably there will be an issue with permissions?

import os

from shutil import copyfile, move
from re import sub
from solcore import config, SOLCORE_ROOT

CUSTOM_PATH = os.path.abspath(config['Others']['custom_mats'].replace('SOLCORE_ROOT', SOLCORE_ROOT))
PARAMETER_PATH = os.path.abspath(config['Parameters']['custom'].replace('SOLCORE_ROOT', SOLCORE_ROOT))


def create_new_material(mat_name, n_source, k_source, parameter_source):
    # create a folder in the custom materials folders
    folder = os.path.join(CUSTOM_PATH, mat_name + '-Material')
    if not os.path.exists(folder) and folder != "":
        os.makedirs(folder)
    else:
        print('This material already exists (or at least, a folder for it).')

    # copy n and k data files to the material's folder
    copyfile(n_source, os.path.join(folder, 'n.txt'))
    copyfile(k_source, os.path.join(folder, 'k.txt'))

    # create the parameter file if it doesn't already exist
    if not os.path.isfile(PARAMETER_PATH):
        open(PARAMETER_PATH, 'a').close()

    # append the parameters for the new material
    fout = open(PARAMETER_PATH, "r")
    existing_parameters = fout.read()
    fout.close()

    if not '[' + mat_name + ']' in existing_parameters:
        # make sure the names match
        fin = open(parameter_source, "r")
        parameters = fin.read() + '\n\n'
        parameters = sub("\[[^]]*\]", lambda x:x.group(0).replace(x.group(0), '['+mat_name+']'), parameters)
        fin.close()
        fout = open(PARAMETER_PATH, "a")
        fout.write(parameters)
        fout.close()
    else:
        print('There are already parameters for this material in the custom parameter file at ' + PARAMETER_PATH)

    # modify the user's config file (in their home folder) to include the relevant paths
    new_entry = mat_name + ' = ' + config['Others']['custom_mats'] + '/' + mat_name + '-Material\n'
    home_folder = os.path.expanduser('~')
    user_config = os.path.join(home_folder, '.solcore_config.txt')
    existing_config = open(user_config, 'r').read()
    if not new_entry in existing_config:

        temp = open('temp', 'w')
        with open(user_config, 'r') as f:
            for line in f:
                if line.startswith('[Materials]'):
                    line = line.strip() + '\n' + new_entry
                temp.write(line)
        temp.close()
        move('temp', user_config)

    else:
        print('A path for this material was already added to the Solcore config file in the home directory.')

    # # Finally add the relevant paths to the MANIFEST.in file.
    # # Don't need to add the full path, just relative to where the MANIFEST file is.
    # Don't think this is necessary if it's installed as a package?
    #
    # path_toadd = folder.replace(SOLCORE_ROOT, 'solcore')
    # new_entry = '\ninclude ' + os.path.join(path_toadd, 'n.txt').replace("\\","/") +' \ninclude ' + os.path.join(path_toadd, 'k.txt').replace("\\","/")
    #
    # MANIFEST_PATH = os.path.join(os.path.dirname(SOLCORE_ROOT), 'MANIFEST.in')
    #
    # existing_manifest = open(MANIFEST_PATH, 'r').read()
    #
    # # check if it's already been added
    # if not new_entry in existing_manifest:
    #     fout = open(MANIFEST_PATH, "a")
    #     fout.write(new_entry)
    #     fout.close()
    # else:
    #     print('A path for this material was already added to the MANIFEST.in file in the package directory')
