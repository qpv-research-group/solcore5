"""
    Another short python 3.4 travesty by Markus.
    Updated by Diego

"""
from configparser import ConfigParser
from solcore.singleton import Singleton
from typing import Optional


class SourceManagedClassDeprecated(metaclass=Singleton):
    """ Base class that manages the sources of Solcore, from the materials to the unit conversion system.
    It must be defined as a derived class from Singleton so only one instance of this class containing all sources
    exists at any time.
    """

    def __init__(self):
        """ Constructor of the class. """
        self.sources = {}
        self.initialise_cache()

    def add_source(self, name, path):
        """ Read a materials database file in ConfigParser format with case-sensitive options. """
        self.sources[name] = path
        self.read(name)

    def remove_source(self, name, reread=True):
        """ Removes a source from the database. """
        del self.sources[name]
        if reread:
            self.read()

    def read(self, name=None):
        """ Reads the sources and add its contents to the database. """
        if name is None:
            self.initialise_cache()
            for sourcename in self.sources.keys():
                self.database.read(self.sources[sourcename], encoding='utf-8')
        else:
            self.database.read(self.sources[name], encoding='utf-8')

    def initialise_cache(self):
        """ Initialises th database, kept in a ConfigParser format. Quite convenient if you think about it. """
        self.database = ConfigParser()

        # Normally the config parser converts options to lower case; this preserves case.
        self.database.optionxform = str


class SourceManagedClass(metaclass=Singleton):
    """ Base class that manages the sources of Solcore. """

    def __init__(self, sources: dict):
        """ Constructor of the class. """
        self.sources = sources
        self.database = ConfigParser()
        self.database.optionxform = str
        self.database.read(self.sources.values(), encoding='utf-8')

    def read(self, source: Optional[str] = None, value: Optional[str] = None):
        """ Updates the information in the database. """
        if source is not None and value is not None:
            self.sources[source] = value
            self.database.read(value, encoding='utf-8')
        elif source in self.sources and value is None:
            del self.sources[source]
            self.database.read(self.sources.values(), encoding='utf-8')
        else:
            self.database.read(self.sources.values(), encoding='utf-8')
