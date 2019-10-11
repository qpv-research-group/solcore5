from solcore import verbose
from typing import Any, Dict, Callable


class Singleton(type):  # META!
    """ A metaclass that restricts children to one instance.
        Attempting to create a second instance simply returns the first instance.
    """
    # All the Singleton instances live here...
    _dufftown: Dict['Singleton', 'Singleton'] = {}

    def __call__(cls, *args: Any, **kwargs: Any) -> 'Singleton':
        if cls not in cls._dufftown:  # ... or if they don't yet, it's about time!
            # print("Instantiating cls", cls)
            cls._dufftown[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._dufftown[cls]

    def breakout(method: 'Singleton') -> 'Singleton':
        """ Decorator for singleton methods that should be available in global scope
        """
        method.break_me_out = None  # or, you know, anything.
        return method

    @staticmethod
    class Redirect_to_singleton_method:
        """ Class that wraps the method that has been defined as 'breakout' with its own class, so it can be called as a
            standalone function.
        """

        def __init__(self, cls: 'Singleton', function: Callable) -> None:
            self.singleton_class = cls  # LABEL COME_BACK_HERE; here instead.;RETURN;
            self.function = function

        def __call__(self, *args: Any, **kwargs: Any) -> str:
            """  Call the singleton method and return the result """
            return self.function(self.singleton_class(), *args, **kwargs)

    class breakoutFunctions:
        """ Decorator for singletons which should make some of their methods available in a certain namespace.
            This is identical to the 'breakoutClass' function defined in 'source_managed_class', but has the option of
            making the function available not in the global scope but somewhere else. """

        def __init__(self, namespace: str) -> None:
            self.namespace_to_breakout_into = namespace

        def __call__(self, cls: 'Singleton') -> 'Singleton':
            """ And this is really the decorator for singletons which should make some of their methods available
            somewhere else
            """
            if verbose:
                print("Making functions from '{}' available "
                      "in the '{}' scope: ".format(cls.__name__, self.namespace_to_breakout_into.__name__), end="")

            # We scan all the methods of the class
            for funcname in (cls.__dict__):
                func = getattr(cls, funcname)

                # And choose those tagged with the 'breakout' decorator
                if hasattr(func, "break_me_out"):
                    # We remove the decorator-ralated atribute...
                    del func.break_me_out
                    # ... and make the function 'func' available in the given scope
                    setattr(self.namespace_to_breakout_into, funcname,
                            Singleton.Redirect_to_singleton_method(cls, func))

                    if verbose:
                        print(funcname, end=" ")

            if verbose:
                print()

            return cls


if __name__ == "__main__":
    class test:
        __name__ = 'ASC_examples'
        pass


    moo = test()


    @Singleton.breakoutFunctions(moo)
    class ForeverAlone(metaclass=Singleton):
        @Singleton.breakout
        def wait_you_havent_instantiated_anything(self) -> str:
            return "YES I DID YOU JUST MISSED IT: " + str(self)

        @Singleton.breakout
        def a_global_instance_method_surely_not(self) -> str:
            return "YES SURELY."

        @Singleton.breakout
        def rofl_this_will_never_work(self, end_symbol: str,
                                      end_number: int) -> str:
            return "WHO'S LAUGHING NOW" + end_symbol * end_number


    a = moo.wait_you_havent_instantiated_anything()
    b = moo.a_global_instance_method_surely_not()
    c = moo.rofl_this_will_never_work('!', end_number=3)

    print("moo.wait_you_havent_instantiated_anything() -> ", a)
    print("moo.a_global_instance_method_surely_not() -> ", b)
    print("moo.rofl_this_will_never_work('!', end_number=3) -> ", c)

    # should output:
    # wait_you_havent_instantiated_anything() ->  YES I DID YOU JUST MISSED IT: <__main__.ForeverAlone object at 0x10c377b90>
    # a_global_instance_method_surely_not() ->  YES SURELY.
    # rofl_this_will_never_work('!', end_number=3) ->  WHO'S LAUGHING NOW!!!
