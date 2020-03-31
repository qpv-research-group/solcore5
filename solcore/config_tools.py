""" This module contains tools for configuring the solcore installation, adding/removing sources with material data,
etc.
"""
from configparser import ConfigParser
import os
import shutil
from typing import Union, Optional, Callable
from pathlib import Path
from warnings import warn
from collections import defaultdict


SOLCORE_ROOT = str(Path(__file__).parent)


class SolcoreConfig:
    def __init__(self, default_config: Union[Path, str], user_config: Union[Path, str]):
        self.default_config = default_config
        self.user_config = user_config
        self.user_folder = str(Path(user_config).parent)
        self._observers = defaultdict(list)
        self._default_data = ConfigParser()
        self._user_data = ConfigParser()
        self._default_data.optionxform = str
        self._user_data.optionxform = str

        self._default_data.read(self.default_config)
        if not os.path.exists(self.user_config):
            self.reset_defaults(confirm=True)
        else:
            self._user_data.read(self.user_config)

        if (
            "version" not in self.sources("Configuration")
            or self.version() != self._default_data["Configuration"]["version"]
        ):
            self.restore_defaults()

    def reset_defaults(self, confirm: bool = False) -> None:
        """ Resets the default Solcore configuration in the user home folder.

        :return: None
        """
        if not confirm:
            response = input(
                "This action will delete any custom Solcore configuration. "
                "Are you sure you want to continue (Y/n)?"
            )

            if response in "Yy":
                confirm = True

        if confirm:
            shutil.copy2(self.default_config, self.user_config)
            self._user_data.read(self.user_config)

    def _save_user_config(self) -> None:
        """ Saves the current user configuration

        :return: None
        """
        with open(self.user_config, "w") as fp:
            self._user_data.write(fp)

    def register_observer(
        self, section: str, callback: Callable[[str], None]
    ) -> None:
        """ Registers and observer to be called when a section changes.

        :param section: Section to observe.
        :param callback: Function to execute when there are changes. It takes the name
        and value of the updated/new source as inputs.
        :return: None
        """
        if section not in self.sections:
            raise KeyError(f"Unknown section: {section}.")
        self._observers[section].append(callback)

    def _notify_observers(self, section: str, source: str, value: str) -> None:
        """ Notifies the observers of a section of a change in a source.

        :param section: Section changed.
        :param name: Source added/updated.
        :param name: New value
        :return: None
        """
        for c in self._observers[section]:
            c(source=source, value=value)

    def __getitem__(self, item):
        """Gets an item from Solcore's configuration."""
        if isinstance(item, tuple) and len(item) == 2:
            return self._user_data[item[0]][item[1]].replace(
                "SOLCORE_ROOT", SOLCORE_ROOT
            )
        elif isinstance(item, str):
            return self._user_data[item]
        else:
            raise KeyError(f"Invalid config option {item}.")

    def __setitem__(self, key, value):
        """Sets an item in Solcore's configuration."""
        if isinstance(key, tuple) and len(key) == 2 and key[0] in self._user_data:
            self._user_data[key[0]][key[1]] = value
            self._save_user_config()
            self._notify_observers(section=key[0], source=key[1], value=value)
        else:
            raise KeyError(f"Invalid config option {key}.")

    @property
    def sections(self):
        """Sections in the configuration file."""
        return self._user_data.sections()

    def sources(self, section):
        """Sources in each the requested section."""
        return self._user_data.options(section)

    def remove_source(self, section: str, name: str) -> None:
        """ General function to remove sources from the configuration files. It checks if
        the source exists and, if so, if it is a default Solcore source. In the later case,
        it disable the source by setting it to an empty string rather than removing it.

        :param section: The section to remove a source from.
        :param name: The name of the source.
        :return: None
        """
        if name in self._default_data.options(section):
            print(
                f"Source {name} is a default Solcore source. "
                "It can not be removed and will be disabled instead."
            )
            self[section, name] = ""
        else:
            self._user_data.remove_option(section, name)
            self._save_user_config()

        print(f"{section}:{name} source removed!")

    def restore_default_source(self, section: str, name: str) -> None:
        """ Restores the default value of a source, assuming the source has a
        default value.

        :param section: The section of the source.
        :param name: The name of the source.
        :return: None
        """
        if name not in self._default_data.options(section):
            options = self._default_data.options(section)
            print(
                f"There is no default value for {section} source: {name}. Available "
                f"default {section} sources are: {options}"
            )
            return

        self[section, name] = self._default_data[section][name]
        print(f"Default Solcore value for {section} source {name} has been restored!")

    def restore_defaults(self) -> None:
        """Restores all the default values without touching user additions.

        :return: None
        """
        for s in self._default_data.sections():
            for o in self._default_data.options(s):
                self._user_data[s][o] = self._default_data[s][o]
        self._save_user_config()

    def version(self) -> str:
        """ Provides the Solcore version

        :return: The version number
        """
        return self["Configuration", "version"]

    def welcome_message(self, show: Optional[bool] = None) -> bool:
        """ Sets if the welcome message must be shown or not

        :param show: True/False for showing/hiding the welcome message
        :return: None
        """
        if show is not None:
            self["Configuration", "welcome_message"] = str(show)
        return self["Configuration", "welcome_message"] == "True"

    def verbose_loading(self, show: Optional[bool] = None) -> bool:
        """ Sets if the loading messages (besides the welcome message) must be shown
        or not

        :param show: True/False for showing/hiding the loading messages
        :return: None
        """
        if show is not None:
            self["Configuration", "verbose_loading"] = str(show)
        return self["Configuration", "verbose_loading"] == "True"

    def spice(self, location: Optional[Union[Path, str]] = None) -> str:
        """ Sets or returns the location of the SPICE executable.

        :param location: The location of the spice executable.
        :return: The currently configured location for the executable
        """
        if location is not None:
            self["External programs", "spice"] = str(location)
        return self["External programs", "spice"]

    def smarts(self, location: Optional[Union[Path, str]] = None) -> str:
        """ Sets or returns the location of the SMARTS executable.

        :param location: The location of the SMARTS distribution.
        :return: The currently configured location for the executable
        """
        if location is not None:
            self["External programs", "smarts"] = str(location)
        return self["External programs", "smarts"]

    def units(
        self, name: Optional[str] = None, location: Optional[Union[Path, str]] = None
    ) -> Union[list, str]:
        """ Adds a Units source to Solcore or returns the value of an existing one.

        If called without arguments, it returns the list of available Unit sources.

        :param name: The name of the source.
        :param location: The full path to the location of the source. The source must
        be a ConfigParser formatted file.
        :return: The path to the source or a list of available sources
        """
        if name is None:
            return self.sources("Units")
        if location is not None:
            self["Units", name] = str(location)
        return self["Units", name]

    def parameters(
        self, name: Optional[str] = None, location: Optional[Union[Path, str]] = None
    ) -> Union[list, str]:
        """ Adds a Parameters source to Solcore or returns the value of an existing one.

        If called without arguments, it returns the list of available Parameter sources.

        :param name: The name of the source.
        :param location: The full path to the location of the source. The source must
        be a ConfigParser formatted file.
        :return: The path to the source or a list of available sources
        """
        if name is None:
            return self.sources("Parameters")
        if location is not None:
            self["Parameters", name] = str(location)
        return self["Parameters", name]

    def materials(
        self, name: Optional[str] = None, location: Optional[Union[Path, str]] = None
    ) -> Union[list, str]:
        """ Adds a Materials source to Solcore or returns the value of an existing one.

        If called without arguments, it returns the list of available Material sources.

        :param name: The name of the source.
        :param location: The full path to the location of the source. The source must
        be a folder containing the n and k
        data of the material. See the the Material System in the manual for more
        information.
        :return: The path to the source or a list of available sources
        """
        if name is None:
            return self.sources("Materials")
        if location is not None:
            self["Materials", name] = str(location)
        return self["Materials", name]

    def __repr__(self) -> str:
        """ Returns the current Solcore configuration

        :return: String with the configuration
        """
        result = ""
        for section in self._user_data.sections():
            result += f"[{section}]\n"
            for option in self._user_data.options(section):
                result += f"{option} = {self._user_data[section][option]}\n"
            result += "\n"
        return result


def reset_defaults(confirm: bool = False) -> None:
    """ Resets the default Solcore configuration in the user home folder.

    :return: None
    """
    warn(
        "Deprecated: this function is deprecated and will be removed in the future."
        "Use solcore.config.reset_defaults instead.",
        FutureWarning,
    )
    from . import config

    config.reset_defaults(confirm=confirm)


def save_user_config() -> None:
    """ Saves the current user configuration

    :return: None
    """
    warn(
        "Deprecated: this function is deprecated and will be removed in the future."
        "Use solcore.config.save_user_config instead.",
        FutureWarning,
    )
    from . import config

    config._save_user_config()


def remove_source(section: str, name: str) -> None:
    """ General function to remove sources from the configuration files. It checks if
    the source exists and, if so, if it is a default Solcore source. In the later case,
    it disable the source by setting it to an empty string rather than removing it.

    :param section: The section to remove a source from.
    :param name: The name of the source.
    :return: None
    """
    warn(
        "Deprecated: this function is deprecated and will be removed in the future."
        "Use solcore.config.remove_source instead.",
        FutureWarning,
    )
    from . import config

    config.remove_source(section=section, name=name)


def add_source(section: str, name: str, location: Union[Path, str]) -> None:
    """ General function to add sources to the configuration files.
    If the source already exists, its value will be replaced by the new one.

    :param section: The section to add a source to.
    :param name: The name of the source.
    :param location: Path to the source (file or folder)
    :return: None
    """
    warn(
        "Deprecated: this function is deprecated and will be removed in the future."
        "Use 'solcore.config[section, option]=value' instead.",
        FutureWarning,
    )
    from . import config

    config.add_source(section=section, name=name, location=location)


def restore_default_source(section: str, name: str) -> None:
    """ Restores the default value of a source, assuming the source has a
    default value.

    :param section: The section of the source.
    :param name: The name of the source.
    :return: None
    """
    warn(
        "Deprecated: this function is deprecated and will be removed in the future."
        "Use solcore.config.restore_default_source instead.",
        FutureWarning,
    )
    from . import config

    config.restore_default_source(section=section, name=name)


def add_units_source(name: str, location: str) -> None:
    """ Adds a Units source to Solcore.

    :param name: The name of the source.
    :param location: The full path to the location of the source.
    The source must be a ConfigParser formatted file.
    :return: None
    """
    warn(
        "Deprecated: this function is deprecated and will be removed in the future."
        "Use solcore.config.units instead.",
        FutureWarning,
    )
    from . import config

    config.units(name, location=location)


def add_parameters_source(name: str, location: str) -> None:
    """ Adds a Parameters source to Solcore.

    :param name: The name of the source.
    :param location: The full path to the location of the source.
    The source must be a ConfigParser formatted file.
    :return: None
    """
    warn(
        "Deprecated: this function is deprecated and will be removed in the future."
        "Use solcore.config.parameters instead.",
        FutureWarning,
    )
    from . import config

    config.parameters(name, location=location)


def add_materials_source(name: str, location: str) -> None:
    """ Adds a Materials source to Solcore.

    :param name: The name of the source.
    :param location: The full path to the location of the source. The source must be a folder containing the n and k
    data of the material. See the the Material System in the manual for more information.
    :return: None
    """
    warn(
        "Deprecated: this function is deprecated and will be removed in the future."
        "Use solcore.config.materials instead.",
        FutureWarning,
    )
    from . import config

    config.materials(name, location=location)


def remove_units_source(name: str) -> None:
    """ Removes a Units source from Solcore.

    :param name: The name of the source.
    :return: None
    """
    warn(
        "Deprecated: this function is deprecated and will be removed in the future."
        "Use solcore.config.remove_source instead.",
        FutureWarning,
    )
    from . import config

    config.remove_source("Units", name)


def remove_parameters_source(name: str) -> None:
    """ Removes a Parameters source from Solcore.

    :param name: The name of the source.
    :return: None
    """
    warn(
        "Deprecated: this function is deprecated and will be removed in the future."
        "Use solcore.config.remove_source instead.",
        FutureWarning,
    )
    from . import config

    config.remove_source("Parameters", name)


def remove_materials_source(name: str) -> None:
    """ Removes a Materials source from Solcore.

    :param name: The name of the source.
    :return: None
    """
    warn(
        "Deprecated: this function is deprecated and will be removed in the future."
        "Use solcore.config.remove_source instead.",
        FutureWarning,
    )
    from . import config

    config.remove_source("Materials", name)


def restore_default_units_source(name: str) -> None:
    """ Restores the default value of a Units source from Solcore.

    :param name: The name of the source.
    :return: None
    """
    warn(
        "Deprecated: this function is deprecated and will be removed in the future."
        "Use solcore.config.restore_default_source instead.",
        FutureWarning,
    )
    from . import config

    config.restore_default_source("Units", name)


def restore_default_parameters_source(name: str) -> None:
    """ Restores the default value of a Parameters source from Solcore.

    :param name: The name of the source.
    :return: None
    """
    warn(
        "Deprecated: this function is deprecated and will be removed in the future."
        "Use solcore.config.restore_default_source instead.",
        FutureWarning,
    )
    from . import config

    config.restore_default_source("Parameters", name)


def restore_default_materials_source(name: str) -> None:
    """ Restores the default value of a Materials source from Solcore.

    :param name: The name of the source.
    :return: None
    """
    warn(
        "Deprecated: this function is deprecated and will be removed in the future."
        "Use solcore.config.restore_default_source instead.",
        FutureWarning,
    )
    from . import config

    config.restore_default_source("Materials", name)


def welcome_message(show: bool) -> None:
    """ Sets if the welcome message must be shown or not

    :param show: True/False for showing/hiding the welcome message
    :return: None
    """
    warn(
        "Deprecated: this function is deprecated and will be removed in the future."
        "Use solcore.config.welcome_message instead.",
        FutureWarning,
    )
    from . import config

    config.welcome_message(show)


def verbose_loading(show: bool) -> None:
    """ Sets if the loading messages (besides the welcome message) must be shown or not

    :param show: True/False for showing/hiding the loading messages
    :return: None
    """
    warn(
        "Deprecated: this function is deprecated and will be removed in the future."
        "Use solcore.config.verbose_loading instead.",
        FutureWarning,
    )
    from . import config

    config.verbose_loading(show)


def set_location_of_spice(location: Union[Path, str]) -> None:
    """ Sets the location of the spice executable. It does not test if it works.

    :param location: The location of the spice executable.
    :return: None
    """
    warn(
        "Deprecated: this function is deprecated and will be removed in the future."
        "Use solcore.config.set_location_of_spice instead.",
        FutureWarning,
    )
    from . import config

    config.set_location_of_spice(location)


def set_location_of_smarts(location: Union[Path, str]) -> None:
    """ Sets the location of the SMARTS distribution (the root folder). It does not test if it works.

    :param location: The location of the SMARTS distribution.
    :return: None
    """
    warn(
        "Deprecated: this function is deprecated and will be removed in the future."
        "Use solcore.config.set_location_of_smarts instead.",
        FutureWarning,
    )
    from . import config

    config.set_location_of_smarts(location)


def get_current_config() -> None:
    """ Prints the current Solcore configuration

    :return: None
    """
    warn(
        "Deprecated: this function is deprecated and will be removed in the future."
        "Use print(solcore.config) instead.",
        FutureWarning,
    )
    from . import config

    print(config)


def check_user_config() -> None:
    """ Checks if there's a user configuration file, asking if it needs to be created.

    :return: None
    """
    warn(
        "Deprecated: this function is deprecated and will be removed in the future."
        "Use solcore.config.check_user_config instead.",
        FutureWarning,
    )
    from . import config

    config.check_user_config()
