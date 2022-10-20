import pytest


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

    @registries.register_optics("final_approximate", available=False)
    def final_approximate_optics(*args, **kwargs):
        pass

    assert "final_approximate" not in registries.OPTICS_METHOD_REGISTRY
