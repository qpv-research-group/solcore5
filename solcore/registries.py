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
from typing import Callable, Dict, Optional
from warnings import warn

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


def register_optics(
    name: str, overwrite: bool = False, reason_to_exclude: Optional[str] = None
) -> Callable:
    """Registers a function that calculates the optics of a solar cell.

    The method must accept as first argument a SolarCell object and can also have as
    input a variable number of parameters needed to perform the calculation, as well as
    a generic **kwargs.

    After running the function, the input SolarCell object will be updated with the
    optical properties, such as the reflection, the transmision, the absorption, etc.

    Args:
        name (str): Name of the method.
        overwrite (bool, optional): If the method should overwrite an existing one with
            the same name. Defaults to False.
        reason_to_exclude (Optional[str], optional): If there is any reason to exclude
            this method from the registry. If not None, the method will be excluded.
            Defaults to None.

    Raises:
        ValueError: If the name of the method exist already in the registry and
            overwrite is False.

    Returns:
        Callable: The inner decorator that will actually register the function.
    """
    if name in OPTICS_METHOD_REGISTRY and not overwrite:
        raise ValueError(
            f"Optics method '{name}' already exist in the registry. "
            "Give it another name or set `overwrite = True`."
        )
    if reason_to_exclude is not None:
        warn(f"Optics method '{name}' will not be available. {reason_to_exclude}")

    def wrap(func: OPTICS_METHOD_SIGNATURE) -> OPTICS_METHOD_SIGNATURE:
        if reason_to_exclude is None:
            OPTICS_METHOD_REGISTRY[name] = func
        return func

    return wrap
