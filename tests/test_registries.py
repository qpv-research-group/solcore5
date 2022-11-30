import pytest


def test_generic_register():
    from solcore import registries
    from typing import Callable

    REGISTRY = {}
    SIGNATURE = Callable

    register_a_solver = registries.generic_register(
        "a-solver",
        registrator_name="Cool features",
        registry=REGISTRY,
        signature=SIGNATURE,
    )

    @register_a_solver
    def solver(*args, **kwargs):
        pass

    assert REGISTRY["a-solver"] == solver

    with pytest.raises(ValueError):

        registries.generic_register(
            "a-solver",
            registrator_name="Cool features",
            registry=REGISTRY,
            signature=SIGNATURE,
        )

    register_a_solver = registries.generic_register(
        "a-solver",
        registrator_name="Cool features",
        registry=REGISTRY,
        signature=SIGNATURE,
        overwrite=True,
    )

    @register_a_solver
    def second_solver(*args, **kwargs):
        pass

    assert REGISTRY["a-solver"] == second_solver


def test_register_action(mocker):
    from solcore import registries

    mock_gr = mocker.patch("solcore.registries.generic_register")
    name = "pre-process"
    overwrite = False
    reason_to_exclude = None

    @registries.register_action(
        name, overwrite=overwrite, reason_to_exclude=reason_to_exclude
    )
    def solver(*args, **kwargs):
        pass

    mock_gr.assert_called_once_with(
        name=name,
        registrator_name="Action",
        registry=registries.ACTIONS_REGISTRY,
        signature=registries.ACTIONS_SIGNATURE,
        overwrite=overwrite,
        reason_to_exclude=reason_to_exclude,
    )


def test_register_optics(mocker):
    from solcore import registries

    mock_gr = mocker.patch("solcore.registries.generic_register")
    name = "custom_optics"
    overwrite = False
    reason_to_exclude = None

    @registries.register_optics(
        name, overwrite=overwrite, reason_to_exclude=reason_to_exclude
    )
    def solver(*args, **kwargs):
        pass

    mock_gr.assert_called_once_with(
        name=name,
        registrator_name="Optics solver",
        registry=registries.OPTICS_METHOD_REGISTRY,
        signature=registries.OPTICS_METHOD_SIGNATURE,
        overwrite=overwrite,
        reason_to_exclude=reason_to_exclude,
    )


def test_register_short_circuit_solver(mocker):
    from solcore import registries

    mock_gr = mocker.patch("solcore.registries.generic_register")
    name = "custom_short_circuit"
    overwrite = False
    reason_to_exclude = None

    @registries.register_short_circuit_solver(
        name, overwrite=overwrite, reason_to_exclude=reason_to_exclude
    )
    def solver(*args, **kwargs):
        pass

    mock_gr.assert_called_once_with(
        name=name,
        registrator_name="Short circuit solver",
        registry=registries.SHORT_CIRCUIT_SOLVER_REGISTRY,
        signature=registries.SHORT_CIRCUIT_SOLVER_SIGNATURE,
        overwrite=overwrite,
        reason_to_exclude=reason_to_exclude,
    )


def test_register_equilibrium_solver(mocker):
    from solcore import registries

    mock_gr = mocker.patch("solcore.registries.generic_register")
    name = "custom_equilibrium"
    overwrite = False
    reason_to_exclude = None

    @registries.register_equilibrium_solver(
        name, overwrite=overwrite, reason_to_exclude=reason_to_exclude
    )
    def solver(*args, **kwargs):
        pass

    mock_gr.assert_called_once_with(
        name=name,
        registrator_name="Equilibrium solver",
        registry=registries.EQUILIBRIUM_SOLVER_REGISTRY,
        signature=registries.EQUILIBRIUM_SOLVER_SIGNATURE,
        overwrite=overwrite,
        reason_to_exclude=reason_to_exclude,
    )


def test_register_iv_solver(mocker):
    from solcore import registries

    mock_gr = mocker.patch("solcore.registries.generic_register")
    name = "custom_iv_solver"
    overwrite = False
    reason_to_exclude = None

    @registries.register_iv_solver(
        name, overwrite=overwrite, reason_to_exclude=reason_to_exclude
    )
    def solver(*args, **kwargs):
        pass

    mock_gr.assert_called_once_with(
        name=name,
        registrator_name="IV solver",
        registry=registries.IV_SOLVER_REGISTRY,
        signature=registries.IV_SOLVER_SIGNATURE,
        overwrite=overwrite,
        reason_to_exclude=reason_to_exclude,
    )
