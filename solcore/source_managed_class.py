"""
    Another short python 3.4 travesty by Markus.

"""
from configparser import ConfigParser  # the chosen file format for the materials data
from solcore.singleton import Singleton


class wrapped:
    """ Class that wraps the method that has been defined as 'breakout' with its own class, so it can be called as a
    standalone function. """

    def __init__(self, cls, function):
        self.singleton_class = cls
        self.function = function

    def __call__(self, *args, **kwargs):
        return self.function(self.singleton_class(), *args, **kwargs)


def breakoutClass(cls):
    """ Decorator definition that indicates the a class contains functions that should be available in the global scope.
    """
    breakout_candidates = list(cls.__dict__.keys())

    # We scan all the methods of the class
    for funcname in (breakout_candidates):
        func = getattr(cls, funcname)

        # And choose those tagged with the 'breakout' decorator
        if hasattr(func, "break_me_out") and func.break_me_out:
            # We remove the decorator-ralated atribute...
            del func.break_me_out
            # ... and make the function 'func' available in the global scope
            globals()[funcname] = wrapped(cls, func)
    return cls


def breakout(plainFunction):
    """ Decorator definition that indicates that a given method of a class should be available in the global scope
    """
    plainFunction.break_me_out = True
    return plainFunction


class SourceManagedClass(metaclass=Singleton):
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


if __name__ == "__main__":
    @breakoutClass
    class M(metaclass=Singleton):
        @breakout
        def moo(self):
            print('it worked!')

        @breakout
        def bah(self, something_else):
            print(something_else)


    # This should work even without creating any instance of the class
    moo()
    bah(3)
