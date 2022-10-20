import pytest


def test_register_action():
    from solcore import registries

    @registries.register_action("pre-process")
    def pre_process_cell(*args, **kwargs):
        pass

    assert "pre-process" in registries.ACTIONS_REGISTRY

    with pytest.raises(ValueError):

        @registries.register_action("pre-process")
        def custom_solve_optics(*args, **kwargs):
            pass

    @registries.register_action("optics", overwrite=True)
    def custom_solve_optics(*args, **kwargs):
        pass

    assert registries.ACTIONS_REGISTRY["optics"] == custom_solve_optics
