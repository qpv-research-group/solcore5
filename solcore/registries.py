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
from typing import Callable, Dict, Optional, Any
from warnings import warn

from .solar_cell import SolarCell
from .structure import Junction
from .state import State


def generic_register(
    name: str,
    registrator_name: str,
    registry: Dict[str, Callable],
    signature: Any,
    overwrite: bool = False,
    reason_to_exclude: Optional[str] = None,
) -> Callable:
    """Generic register that can be used by the other specific ones.

    Args:
        name: Name of the function to register
        registrator_name (str): Name of the action of the specifric register, eg.
        "Optics solver".
        registry: Registry in which to store the registered function.
        signature: Signature of that function.
        overwrite: If the method should overwrite an existing one with
            the same name. Defaults to False.
        reason_to_exclude: If there is any reason to exclude
            this method from the registry. If not None, the method will be excluded.
            Defaults to None.

    Raises:
        ValueError: If the name of the method exist already in the registry and
            overwrite is False.

    Returns:
        Callable: The inner decorator that will actually register the function.
    """
    if name in registry and not overwrite:
        raise ValueError(
            f"{registrator_name.capitalize()} '{name}' already exist in the registry. "
            "Give it another name or set `overwrite = True`."
        )
    if reason_to_exclude is not None:
        warn(
            f"{registrator_name.capitalize()} '{name}' will not be available. "
            f"{reason_to_exclude}"
        )

    def wrap(func: signature) -> signature:
        if reason_to_exclude is None:
            registry[name] = func
        return func

    return wrap


ACTIONS_SIGNATURE = Callable[[SolarCell, State], None]
ACTIONS_REGISTRY: Dict[str, ACTIONS_SIGNATURE] = {}


def register_action(
    name: str, overwrite: bool = False, reason_to_exclude: Optional[str] = None
) -> Callable:
    """Registers a global action to be executed over the whole cell.

    Examples of these actions could include soling the optics of the cell, calculate the
    IV curve or the quantum efficiency.

    Args:
        name (str): Name of the action.
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
    return generic_register(
        name=name,
        registrator_name="Action",
        registry=ACTIONS_REGISTRY,
        signature=ACTIONS_SIGNATURE,
        overwrite=overwrite,
        reason_to_exclude=reason_to_exclude,
    )


OPTICS_METHOD_SIGNATURE = Callable[[SolarCell, Any], None]
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
    return generic_register(
        name=name,
        registrator_name="Optics solver",
        registry=OPTICS_METHOD_REGISTRY,
        signature=OPTICS_METHOD_SIGNATURE,
        overwrite=overwrite,
        reason_to_exclude=reason_to_exclude,
    )


SHORT_CIRCUIT_SOLVER_SIGNATURE = Callable[[Junction, Any], None]
SHORT_CIRCUIT_SOLVER_REGISTRY: Dict[str, SHORT_CIRCUIT_SOLVER_SIGNATURE] = {}


def register_short_circuit_solver(
    name: str, overwrite: bool = False, reason_to_exclude: Optional[str] = None
) -> Callable:
    """Registers a function that solves a junction under short circuit conditions.

    The solver must accept as first argument a Junction object and can also have as
    input a variable number of parameters needed to perform the calculation, as well as
    a generic **kwargs.

    After running the function, the input Junction object will be updated with the
    bandstructure and recombination profile at short circuit will be available in a new
    `short_circuit_data` attribute of the junction object as a State object.

    Args:
        name (str): Name of the solver.
        overwrite (bool, optional): If the method should overwrite an existing one with
            the same name. Defaults to False.
        reason_to_exclude (Optional[str], optional): If there is any reason to exclude
            this solver from the registry. If not None, the method will be excluded.
            Defaults to None.

    Raises:
        ValueError: If the name of the solver exists already in the registry and
            overwrite is False.

    Returns:
        Callable: The inner decorator that will actually register the function.
    """
    return generic_register(
        name=name,
        registrator_name="Short circuit solver",
        registry=SHORT_CIRCUIT_SOLVER_REGISTRY,
        signature=SHORT_CIRCUIT_SOLVER_SIGNATURE,
        overwrite=overwrite,
        reason_to_exclude=reason_to_exclude,
    )


EQUILIBRIUM_SOLVER_SIGNATURE = Callable[[Junction, Any], None]
EQUILIBRIUM_SOLVER_REGISTRY: Dict[str, EQUILIBRIUM_SOLVER_SIGNATURE] = {}


def register_equilibrium_solver(
    name: str, overwrite: bool = False, reason_to_exclude: Optional[str] = None
) -> Callable:
    """Registers a function that solves a junction under equilibrium conditions.

    The solver must accept as first argument a Junction object and can also have as
    input a variable number of parameters needed to perform the calculation, as well as
    a generic **kwargs.

    After running the function, the input Junction object will be updated with the
    bandstructure profile at equilibrium, available in a new `equilibrium_data`
    attribute of the junction object as a State object.

    Args:
        name (str): Name of the solver.
        overwrite (bool, optional): If the method should overwrite an existing one with
            the same name. Defaults to False.
        reason_to_exclude (Optional[str], optional): If there is any reason to exclude
            this solver from the registry. If not None, the method will be excluded.
            Defaults to None.

    Raises:
        ValueError: If the name of the solver exists already in the registry and
            overwrite is False.

    Returns:
        Callable: The inner decorator that will actually register the function.
    """
    return generic_register(
        name=name,
        registrator_name="Equilibrium solver",
        registry=EQUILIBRIUM_SOLVER_REGISTRY,
        signature=EQUILIBRIUM_SOLVER_SIGNATURE,
        overwrite=overwrite,
        reason_to_exclude=reason_to_exclude,
    )


IV_SOLVER_SIGNATURE = Callable[[Junction, Any], None]
IV_SOLVER_REGISTRY: Dict[str, EQUILIBRIUM_SOLVER_SIGNATURE] = {}


def register_iv_solver(
    name: str, overwrite: bool = False, reason_to_exclude: Optional[str] = None
) -> Callable:
    """Registers a function that solves the IV curve of an independent junction.

    The solver must accept as first argument a Junction object and can also have as
    input a variable number of parameters needed to perform the calculation, as well as
    a generic **kwargs.

    After running the function, the input Junction object will be updated with the
    voltage, the current, a linear interpolator to calculate the IV curve at any voltage
    and, potentialy, other auxiliary currents - related to different recombination
    mechanisms or correspoinding to different regions of the cell.

    Args:
        name (str): Name of the solver.
        overwrite (bool, optional): If the method should overwrite an existing one with
            the same name. Defaults to False.
        reason_to_exclude (Optional[str], optional): If there is any reason to exclude
            this solver from the registry. If not None, the method will be excluded.
            Defaults to None.

    Raises:
        ValueError: If the name of the solver exists already in the registry and
            overwrite is False.

    Returns:
        Callable: The inner decorator that will actually register the function.
    """
    return generic_register(
        name=name,
        registrator_name="IV solver",
        registry=IV_SOLVER_REGISTRY,
        signature=IV_SOLVER_SIGNATURE,
        overwrite=overwrite,
        reason_to_exclude=reason_to_exclude,
    )
