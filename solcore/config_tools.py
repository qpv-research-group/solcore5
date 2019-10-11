""" This module contains tools for configuring the solcore installation, adding/removing sources with material data,
etc.
"""
from configparser import ConfigParser
import os
import shutil
import solcore
import glob

home_folder: str = os.path.expanduser('~')
user_config: str = os.path.join(home_folder, '.solcore_config.txt')
user_config_data = ConfigParser()
user_config_data.read(user_config)

default_config_data = ConfigParser()
default_config_data.read(solcore.default_config)


def reset_defaults(confirm: bool = False) -> None:
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


def save_user_config() -> None:
    """ Saves the current user configuration

    :return: None
    """
    with open(user_config, 'w') as fp:
        user_config_data.write(fp)


def remove_source(section: str, name: str) -> None:
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


def add_source(section: str, name: str, location: str) -> None:
    """ General function to add sources to the configuration files. If the source already exists, its value will be
    replaced by the new one.

    :param section: The section to add a source to.
    :param name: The name of the source.
    :return: None
    """
    user_config_data.set(section, name, value=location)
    save_user_config()

    print('{}:{} source added!'.format(section, name))


def restore_default_source(section: str, name: str) -> None:
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


def add_units_source(name: str, location: str) -> None:
    """ Adds a Units source to Solcore.

    :param name: The name of the source.
    :param location: The full path to the location of the source. The source must be a ConfigParser formatted file.
    :return: None
    """
    add_source('Units', name, location=location)


def add_parameters_source(name: str, location: str) -> None:
    """ Adds a Parameters source to Solcore.

    :param name: The name of the source.
    :param location: The full path to the location of the source. The source must be a ConfigParser formatted file.
    :return: None
    """
    add_source('Parameters', name, location=location)


def add_materials_source(name: str, location: str) -> None:
    """ Adds a Materials source to Solcore.

    :param name: The name of the source.
    :param location: The full path to the location of the source. The source must be a folder containing the n and k
    data of the material. See the the Material System in the manual for more information.
    :return: None
    """
    add_source('Materials', name, location=location)


def remove_units_source(name: str) -> None:
    """ Removes a Units source from Solcore.

    :param name: The name of the source.
    :return: None
    """
    remove_source('Units', name)


def remove_parameters_source(name: str) -> None:
    """ Removes a Parameters source from Solcore.

    :param name: The name of the source.
    :return: None
    """
    remove_source('Parameters', name)


def remove_materials_source(name: str) -> None:
    """ Removes a Materials source from Solcore.

    :param name: The name of the source.
    :return: None
    """
    remove_source('Materials', name)


def restore_default_units_source(name: str) -> None:
    """ Restores the default value of a Units source from Solcore.

    :param name: The name of the source.
    :return: None
    """
    restore_default_source('Units', name)


def restore_default_parameters_source(name: str) -> None:
    """ Restores the default value of a Parameters source from Solcore.

    :param name: The name of the source.
    :return: None
    """
    restore_default_source('Parameters', name)


def restore_default_materials_source(name: str) -> None:
    """ Restores the default value of a Materials source from Solcore.

    :param name: The name of the source.
    :return: None
    """
    restore_default_source('Materials', name)


def welcome_message(show: bool) -> None:
    """ Sets if the welcome message must be shown or not

    :param show: True/False for showing/hiding the welcome message
    :return: None
    """
    user_config_data['Configuration']['welcome_message'] = int(show)


def verbose_loading(show: bool) -> None:
    """ Sets if the loading messages (besides the welcome message) must be shown or not

    :param show: True/False for showing/hiding the loading messages
    :return: None
    """
    user_config_data['Configuration']['verbose_loading'] = int(show)


def set_location_of_spice(location: str) -> None:
    """ Sets the location of the spice executable. It does not test if it works.

    :param location: The location of the spice executable.
    :return: None
    """
    user_config_data['External programs']['spice'] = location
    save_user_config()


def set_location_of_smarts(location: str) -> None:
    """ Sets the location of the SMARTS distribution (the root folder). It does not test if it works.

    :param location: The location of the SMARTS distribution.
    :return: None
    """
    user_config_data['External programs']['smarts'] = location
    save_user_config()


def get_current_config() -> None:
    """ Prints the current Solcore configuration

    :return: None
    """

    for section in user_config_data.sections():
        print('[{}]'.format(section))
        for option in user_config_data.options(section):
            print('{} = {}'.format(option, user_config_data.get(section, option)))

        print()


def check_user_config() -> None:
    """ Checks if there's a user configuration file, asking if it needs to be created.

    :return: None
    """
    if len(user_config_data.sections()) == 0:
        response = input('No user configuration was detected. Do you want to create one (Y/n)?')

        if response in 'Yy':
            reset_defaults(True)
            user_config_data = ConfigParser()
            user_config_data.read(user_config)
