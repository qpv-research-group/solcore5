""" This module contains tools for configuring the solcore installation, adding/removing sources with material data,
etc.
"""
from configparser import ConfigParser
import os
import shutil
import solcore
import glob

home_folder = os.path.expanduser('~')
user_config = os.path.join(home_folder, '.solcore_config.txt')
user_config_data = ConfigParser()
user_config_data.read(user_config)

default_config_data = ConfigParser()
default_config_data.read(solcore.default_config)


def reset_defaults(confirm=False):
    """ Resets the default Solcore configuration in the user home folder.

    :return: None
    """
    global user_config_data
    if not confirm:
        response = input('This action will delete any customisation of the Solcore configuration. '
                         'Are you sure you want to continue (Y/n)?')

        if response in 'Yy':
            confirm = True

    if confirm:
        shutil.copy2(solcore.default_config, user_config)

        user_config_data.read(user_config)
        user_config_data.remove_option(section='Configuration', option='version')
        save_user_config()

        print('Default Solcore configuration has been restored!')


def save_user_config():
    """ Saves the current user configuration

    :return: None
    """
    with open(user_config, 'w') as fp:
        user_config_data.write(fp)


def remove_source(section, name):
    """ General function to remove sources from the configuration files. It checks if the source exists and, if so, if
    it is a default Solcore source. In the later case, it disable the source by setting it to an empty string rather
    than removing it.

    :param section: The section to remove a source from.
    :param name: The name of the source.
    :return: None
    """
    if name not in user_config_data.options(section):
        print('Source {} does not exists. Valid options are: {}'.format(name, user_config_data.options(section)))
        return

    if name in default_config_data.options(section):
        print('Source {} is a default Solcore source. It can not be removed and will be disabled instead.'.format(name))
        add_units_source(name, '')
        return

    user_config_data.remove_option(section, name)
    save_user_config()

    print('{}:{} source removed!'.format(section, name))


def add_source(section, name, location):
    """ General function to add sources to the configuration files. If the source already exists, its value will be
    replaced by the new one.

    :param section: The section to add a source to.
    :param name: The name of the source.
    :return: None
    """
    user_config_data.set(section, name, value=location)
    save_user_config()

    print('{}:{} source added!'.format(section, name))


def restore_default_source(section, name):
    """ Restores the default value of a source, assuming the source has a default value.

    :param section: The section of the source.
    :param name: The name of the source.
    :return: None
    """
    if name not in default_config_data.options(section):
        print('There is no default value for {} source: {}. '.format(section, name))
        print('Available default {} sources are: {}'.format(section, default_config_data.options(section)))
        return

    user_config_data[section][name] = default_config_data[section][name]
    save_user_config()

    print('Default Solcore value for {} source {} has been restored!'.format(section, name))


def add_units_source(name, location):
    """ Adds a Units source to Solcore.

    :param name: The name of the source.
    :param location: The full path to the location of the source. The source must be a ConfigParser formatted file.
    :return: None
    """
    add_source('Units', name, location=location)


def add_parameters_source(name, location):
    """ Adds a Parameters source to Solcore.

    :param name: The name of the source.
    :param location: The full path to the location of the source. The source must be a ConfigParser formatted file.
    :return: None
    """
    add_source('Parameters', name, location=location)


def add_materials_source(name, location):
    """ Adds a Materials source to Solcore.

    :param name: The name of the source.
    :param location: The full path to the location of the source. The source must be a folder containing the n and k.txt
    data of the material. See the the Material System in the manual for more information.
    :return: None
    """
    add_source('Materials', name, location=location)


def remove_units_source(name):
    """ Removes a Units source from Solcore.

    :param name: The name of the source.
    :return: None
    """
    remove_source('Units', name)


def remove_parameters_source(name):
    """ Removes a Parameters source from Solcore.

    :param name: The name of the source.
    :return: None
    """
    remove_source('Parameters', name)


def remove_materials_source(name):
    """ Removes a Materials source from Solcore.

    :param name: The name of the source.
    :return: None
    """
    remove_source('Materials', name)


def restore_default_units_source(name):
    """ Restores the default value of a Units source from Solcore.

    :param name: The name of the source.
    :return: None
    """
    restore_default_source('Units', name)


def restore_default_parameters_source(name):
    """ Restores the default value of a Parameters source from Solcore.

    :param name: The name of the source.
    :return: None
    """
    restore_default_source('Parameters', name)


def restore_default_materials_source(name):
    """ Restores the default value of a Materials source from Solcore.

    :param name: The name of the source.
    :return: None
    """
    restore_default_source('Materials', name)


def welcome_message(show):
    """ Sets if the welcome message must be shown or not

    :param show: True/False for showing/hiding the welcome message
    :return: None
    """
    user_config_data['Configuration']['welcome_message'] = int(show)


def verbose_loading(show):
    """ Sets if the loading messages (besides the welcome message) must be shown or not

    :param show: True/False for showing/hiding the loading messages
    :return: None
    """
    user_config_data['Configuration']['verbose_loading'] = int(show)

#
# def get_solcore_examples(destination=home_folder):
#     """ Copies all Solcore examples to a user designated folder.
#
#     :param destination:
#     :return:
#     """
#     examples_folder = os.path.join(solcore.SOLCORE_ROOT, 'examples')
#     destination = os.path.join(destination, 'solcore', 'examples')
#
#     if os.path.isdir(destination):
#         response = input('The destination folder for the examples already exists: \n\t{}\nIt will be deleted. '
#                          'Do you want to continue (Y/n)?'.format(destination))
#         if response not in 'Yy':
#             return
#
#         shutil.rmtree(destination)
#
#     shutil.copytree(examples_folder, destination)
#
#     # We remove all the files/folders that are unnecesary
#     pushdir = os.getcwd()
#     os.chdir(destination)
#     init_files = glob.glob(os.path.join('**', '__init__.py'), recursive=True)
#     for f in init_files:
#         os.remove(f)
#     init_files = glob.glob(os.path.join('**', '__pycache__'), recursive=True)
#     for f in init_files:
#         shutil.rmtree(f)
#     os.chdir(pushdir)


def set_location_of_spice(location):
    """ Sets the location of the spice executable. It does not test if it works.

    :param location: The location of the spice executable.
    :return: None
    """
    user_config_data['External programs']['spice'] = location
    save_user_config()


def set_location_of_smarts(location):
    """ Sets the location of the SMARTS distribution (the root folder). It does not test if it works.

    :param location: The location of the SMARTS distribution.
    :return: None
    """
    user_config_data['External programs']['smarts'] = location
    save_user_config()


def set_fortran_compiler(location):
    """ Sets the fotran compiler. It does not test if it works.

    :param location: The fortran compiler or the path to the fortran compiler, if it is not in the system path.
    :return: None
    """
    user_config_data['External programs']['fortran'] = location
    save_user_config()


def compile_poisson_drift_diffusion():
    """ Compiles the fortran-based poisson-drift-diffusion solver, wrapping the resulting library using F2Py to be accessible from Python.

    :return: none
    """

    from .poisson_drift_diffusion.driftdiffusion_compiler import check_ddModel_library_ok

    check_ddModel_library_ok(force=True)


def get_current_config():
    """ Prints the current Solcore configuration

    :return: None
    """

    for section in user_config_data.sections():
        print('[{}]'.format(section))
        for option in user_config_data.options(section):
            print('{} = {}'.format(option, user_config_data.get(section, option)))

        print()
#
#
# def open_documentation():
#     """ Opens Solcore documentation in a new tab in your web browser.
#
#     :return: None
#     """
#     import webbrowser
#
#     url = 'file:' + os.path.join(solcore.SOLCORE_ROOT, 'documentation', 'Solcore.html')
#     webbrowser.open(url)
#

if len(user_config_data.sections()) == 0:
    response = input('No user configuration was detected. Do you want to create one (Y/n)?')

    if response in 'Yy':
        reset_defaults(True)
        user_config_data = ConfigParser()
        user_config_data.read(user_config)
