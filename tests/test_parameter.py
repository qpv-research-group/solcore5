from typing import Tuple

from pytest import mark, raises
from unittest.mock import MagicMock

from solcore.parameter import Parameter, ParameterSource


class TestParameter:
    def test_description(self):
        from solcore.parameter import Parameter

        d = "Great parameter"
        var = Parameter(42, "k", description=d)
        assert var.description == d
        assert var.d == d

        var = Parameter(42, "k")
        assert var.description == ""
        assert var.d == ""

    def test_reference(self):
        from solcore.parameter import Parameter

        ref = "Some paper"
        var = Parameter(42, "k", reference=ref)
        assert var.reference == (ref,)
        assert var.r == (ref,)

        ref = ("Some paper", "another paper")
        var = Parameter(42, "k", reference=ref)
        assert var.reference == ref
        assert var.r == ref

        var = Parameter(42, "k")
        assert var.reference == ()
        assert var.r == ()


class TestParameterManager:
    def test_gather_sources(self, parameter_manager):
        assert len(parameter_manager._known_sources) == 0
        parameter_manager.gather_sources()
        assert len(parameter_manager._known_sources) > 0

    def test_initialize(self, parameter_manager):
        class DummySource:
            load_source = MagicMock(return_value="A source")

        parameter_manager._known_sources["my source"] = DummySource
        assert len(parameter_manager.sources) == 0
        parameter_manager.initialize()
        assert parameter_manager.sources["my source"] == "A source"

    def test_add_source(self, parameter_manager):
        class Dummy:
            cache_clear = MagicMock()

        ps = parameter_manager
        ps._normalise_source = Dummy()
        ps._validate_source = Dummy()

        ps.add_source("my source", None)
        assert ps._known_sources["my source"] is None
        ps._normalise_source.cache_clear.assert_called()
        ps._validate_source.cache_clear.assert_called()

        with raises(ValueError):
            ps.add_source("my source", None)

    def test_known_sources(self, parameter_manager):

        sources = ("source 1", "source 2")
        ps = parameter_manager
        for s in sources:
            ps._known_sources[s] = None
        assert ps.known_sources == sources

    def test_get_parameter(self, parameter_manager):
        from solcore.parameter import ParameterError

        ps = parameter_manager
        sources = ("source 1", "source 2")
        ps._normalise_source = MagicMock(return_value=sources)

        class DummySource:
            source_1 = MagicMock()

            def __init__(self, name):
                self.name = name

            def get_parameter(self, material, parameter, **kwargs):
                if parameter == "the answer" and self.name == "source 1":
                    self.source_1()
                    raise ParameterError
                elif parameter == "the answer" and self.name == "source 2":
                    return 42
                else:
                    raise ParameterError

        ps._load_source = MagicMock(side_effect=[DummySource(s) for s in sources])
        value = ps.get_parameter("dark matter", "the answer")
        DummySource.source_1.assert_called_once()
        assert value == 42

        ps._load_source = MagicMock(side_effect=[DummySource(s) for s in sources])
        with raises(ParameterError):
            ps.get_parameter("dark matter", "stupid question")

    def test_get_multiple_parameters(self, parameter_manager):
        ps = parameter_manager
        sources = ("source 1", "source 2")
        ps._normalise_source = MagicMock(return_value=sources)

        # include and exclude provided results in error
        with raises(ValueError):
            ps.get_multiple_parameters(
                "dark matter", include=("param 1",), exclude=("param 2")
            )

        # Only included
        ps.get_parameter = MagicMock(return_value=42)
        included = ("param 1", "param 2")
        expected = {"param 1": 42, "param 2": 42}
        actual = ps.get_multiple_parameters("dark matter", include=included)
        assert actual == expected

        # All with excluded
        class DummySource:
            source_1 = MagicMock()

            def __init__(self, name):
                self.name = name

            def parameters(self, material):
                if self.name == "source 1":
                    return "param 1", "param 2", "param 3"
                elif self.name == "source 2":
                    return tuple(("param 2",))

        ps._load_source = MagicMock(side_effect=[DummySource(s) for s in sources])
        excluded = ("param 1",)
        expected = {"param 2": 42, "param 3": 42}
        actual = ps.get_multiple_parameters("dark matter", exclude=excluded)
        assert actual == expected

    def test__validate_source(self, parameter_manager):
        from solcore.parameter import ParameterSourceError

        sources = ("source 1", "source 2")
        ps = parameter_manager
        for s in sources:
            ps._known_sources[s] = None

        ps._validate_source("source 1")
        with raises(ParameterSourceError):
            ps._validate_source("source 3")

    def test__load_source(self, parameter_manager):
        ps = parameter_manager
        ps.initialize = MagicMock()
        ps._validate_source = MagicMock()

        with raises(KeyError):
            ps._load_source("source 1")
            ps.initialize.assert_called_once()
            ps._validate_source.assert_called_once()

        ps.initialize.reset_mock()
        ps._validate_source.reset_mock()
        ps.sources["source 1"] = 42
        assert ps._load_source("source 1") == 42
        ps.initialize.assert_not_called()
        ps._validate_source.assert_called_once()

    def test__normalise_source(self, parameter_manager):
        ps = parameter_manager
        ps.initialize = MagicMock()
        ps._validate_source = MagicMock()

        class DummySource:
            def __init__(self, priority):
                self.priority = priority

        sources = ("source 1", "source 2", "source 3")
        for i, s in enumerate(sources):
            ps._known_sources[s] = DummySource
            ps.sources[s] = DummySource(i)

        # No specific source given -> all returned sorted by priority
        expected = ("source 3", "source 2", "source 1")
        assert ps._normalise_source(None) == expected

        # A single source provided, returned as a tuple
        expected = ("source 3",)
        assert ps._normalise_source("source 3") == expected

        # A couple of sources provided, returned iun the order requested by the user
        expected = ("source 1", "source 3")
        assert ps._normalise_source(("source 1", "source 3")) == expected


class TestParameterSourceBase:
    def test__init_subclass__(self, parameter_manager):
        from solcore.parameter import ParameterSourceBase

        parameter_manager.add_source = MagicMock()

        class NewSource(ParameterSourceBase):
            name = "my source"

        parameter_manager.add_source.assert_called_once_with("my source", NewSource)
        parameter_manager.add_source.reset_mock()

        class NewSource2(ParameterSourceBase):
            name = ["my source 1", "my source 2"]

        assert parameter_manager.add_source.call_count == 2

    def test_properties(self, parameter_manager):
        from solcore.parameter import ParameterSourceBase

        parameter_manager.add_source = MagicMock()

        class NewSource(ParameterSourceBase):
            name = "my class"
            _priority = 10

            @classmethod
            def load_source(cls, source_name: str = "") -> ParameterSource:
                pass

            @property
            def materials(self) -> Tuple[str, ...]:
                pass

            def parameters(self, material: str) -> Tuple[str, ...]:
                pass

            def get_parameter(self, material: str, parameter: str,
                              **kwargs) -> Parameter:
                pass

        assert NewSource().parman == parameter_manager
        assert NewSource().priority == 10


@mark.parametrize("bow", [-1, 1])
def test_alloy_parameter(bow):
    import numpy as np
    from solcore.parameter import alloy_parameter

    p0 = 2
    p1 = 5
    x = np.random.random()
    no_bow = alloy_parameter(p0, p1, x)
    with_bow = alloy_parameter(p0, p1, x, bow)
    assert with_bow <= no_bow if bow >= 0 else with_bow > no_bow
