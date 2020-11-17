from pytest import raises
from unittest.mock import MagicMock


def test_eg():
    assert False


def test_eg_gamma():
    assert False


def test_eg_x():
    assert False


def test_eg_l():
    assert False


def test_band_gap():
    assert False


def test_lowest_band():
    assert False


def test_eff_mass_split_off():
    assert False


def test_eff_mass_hh_z():
    assert False


def test_eff_mass_hh_110():
    assert False


def test_eff_mass_hh_111():
    assert False


def test_eff_mass_lh_z():
    assert False


def test_eff_mass_lh_110():
    assert False


def test_eff_mass_lh_111():
    assert False


def test_eff_mass_electron():
    assert False


def test_permittivity():
    assert False


class TestCalculableParameters:
    def test__getitem__(self):
        from solcore.parameter_sources import CalculableParameters
        from solcore.parameter import ParameterError

        cp = CalculableParameters()
        with raises(ParameterError):
            cp["param1"]

        assert cp._params["eg_gamma"] == cp["eg_gamma"]

    def test___getattr__(self):
        from solcore.parameter_sources import CalculableParameters
        from solcore.parameter import ParameterError

        cp = CalculableParameters()
        with raises(ParameterError):
            cp.param1

        assert cp._params["eg_gamma"] == cp.eg_gamma

    def test_load_source(self):
        from solcore.parameter_sources import CalculableParameters

        cp = CalculableParameters()
        assert cp == cp.load_source()

    def test_register_calculable(self):
        from solcore.parameter_sources import CalculableParameters
        from solcore.parameter import ParameterSourceError

        cp = CalculableParameters()
        assert len(cp._params) > 0
        assert tuple(cp._params.keys()) == tuple(cp._descriptions.keys())

        with raises(ParameterSourceError):

            @CalculableParameters.register_calculable
            def eg_gamma():
                pass

    def test_materials(self):
        from solcore.parameter_sources import CalculableParameters

        cp = CalculableParameters()
        assert not cp._warned
        assert cp.materials == ()
        assert cp._warned

    def test_parameters(self):
        from solcore.parameter_sources import CalculableParameters

        cp = CalculableParameters()
        assert cp.parameters() == tuple((p for p in cp._params))

    def test_list_arguments(self):
        from solcore.parameter_sources import CalculableParameters

        cp = CalculableParameters()
        expected = ("T", "eg0_gamma", "alpha_gamma", "beta_gamma")
        assert cp.list_arguments("eg_gamma") == expected

    def test_get_parameter(self):
        from solcore.parameter_sources import CalculableParameters
        from solcore.parameter import ParameterError, Parameter

        cp = CalculableParameters()

        @CalculableParameters.register_calculable(description="The answer")
        def dummy_param():
            return 42

        out = cp.get_parameter("Dark matter", "dummy_param")
        assert out.m == 42
        assert out.d == "The answer"
        assert out.r == ("Calculable",)

        @CalculableParameters.register_calculable(description="Twice the answer")
        def dummy_param_2(dummy_param):
            return dummy_param * 2

        out = cp.get_parameter("Dark matter", "dummy_param_2")
        assert out.m == 84
        assert out.d == "Twice the answer"
        assert out.r == ("Calculable",)

        def get_param(material, param, **kwargs):
            if param == "param_16":
                return Parameter(16, reference="SomePaper")
            else:
                raise ParameterError()

        cp.parman.get_parameter = get_param

        @CalculableParameters.register_calculable(description="Answer is 16")
        def dummy_param_16(param_16):
            return param_16

        out = cp.get_parameter("Dark matter", "dummy_param_16")
        assert out.m == 16
        assert out.d == "Answer is 16"
        assert set(out.r) == {"Calculable", "SomePaper"}

        with raises(ParameterError):
            cp.get_parameter("Dark matter", "other_dummy_param")

        # Let's get rid of all the dummies...
        dummies = (d for d in cp.parameters() if "dummy" in d)
        for d in dummies:
            cp._params.pop(d)
