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

    assert "a-solver" in REGISTRY
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


def test_register_action():
    from solcore import registries

    @registries.register_action("pre-process")
    def pre_process_cell(*args, **kwargs):
        pass

    assert "pre-process" in registries.ACTIONS_REGISTRY

    with pytest.raises(ValueError):

        @registries.register_action("pre-process")
        def custom_pre_process_cell(*args, **kwargs):
            pass

    @registries.register_action("pre-process", overwrite=True)
    def another_pre_process_cell(*args, **kwargs):
        pass

    assert registries.ACTIONS_REGISTRY["pre-process"] == another_pre_process_cell


def test_register_optics():
    from solcore import registries

    @registries.register_optics("approximate")
    def approximate_optics(*args, **kwargs):
        pass

    assert "approximate" in registries.OPTICS_METHOD_REGISTRY

    with pytest.raises(ValueError):

        @registries.register_optics("approximate")
        def custom_approximate_optics(*args, **kwargs):
            pass

    @registries.register_optics("approximate", overwrite=True)
    def another_approximate_optics(*args, **kwargs):
        pass

    assert (
        registries.OPTICS_METHOD_REGISTRY["approximate"] == another_approximate_optics
    )

    @registries.register_optics("final_approximate", reason_to_exclude="Doesn't work")
    def final_approximate_optics(*args, **kwargs):
        pass

    assert "final_approximate" not in registries.OPTICS_METHOD_REGISTRY
