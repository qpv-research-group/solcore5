"""Registries are dictionaries that put together Solcore's functionality that share the
same purpose, eg. calculate solar cell properties, solve the electrical properties of
junctions or solve their optical properties. They can also perform validation and be
used to extend Solcore's functionality without the need of modifying the source code.

The main advantage of registries is that they allow you to define your own functionality
outside of the solcore code base, simplifying the process of trying new things:

>>> from solcore import registries
>>> @registries.register_action("pre-process")
... def pre_process_cell(*args, **kwargs):
...     pass
>>> "pre-process" in registries.ACTIONS_REGISTRY
True

After this, `pre-process` will be one of the `task` that you can run when executing the
`solar_cell_solver` funciton.

You can also overwrite existing functionality in case you want to implement your own
alternative to an existing approach:

>>> from solcore import registries
>>> @registries.register_action("optics", overwrite=True)
... def custom_solve_optics(*args, **kwargs):
...     pass
>>> registries.ACTIONS_REGISTRY["optics"] == custom_solve_optics
True

"""
from typing import Callable, Dict

from .solar_cell import SolarCell
from .state import State

ACTIONS_SIGNATURE = Callable[[SolarCell, State], None]
ACTIONS_REGISTRY: Dict[str, ACTIONS_SIGNATURE] = {}


def register_action(name: str, overwrite: bool = False) -> Callable:
    if name in ACTIONS_REGISTRY and not overwrite:
        raise ValueError(
            f"Action '{name}' already exist in the registry."
            "Give it another name or set `overwrite = True`."
        )

    def wrap(func: ACTIONS_SIGNATURE) -> ACTIONS_SIGNATURE:
        ACTIONS_REGISTRY[name] = func
        return func

    return wrap


OPTICS_METHOD_SIGNATURE = Callable[[SolarCell, State], None]
OPTICS_METHOD_REGISTRY: Dict[str, OPTICS_METHOD_SIGNATURE] = {}


def register_optics(name: str, overwrite: bool = False) -> Callable:
    if name in OPTICS_METHOD_REGISTRY and not overwrite:
        raise ValueError(
            f"Optics method '{name}' already exist in the registry."
            "Give it another name or set `overwrite = True`."
        )

    def wrap(func: OPTICS_METHOD_SIGNATURE) -> OPTICS_METHOD_SIGNATURE:
        OPTICS_METHOD_REGISTRY[name] = func
        return func

    return wrap
