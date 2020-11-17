from pytest import raises, mark, approx
from unittest.mock import patch, MagicMock


def test_eg():
    from pint import Quantity
    from solcore.parameter_sources.calculable_parameters import eg

    T = Quantity(298, "K")
    eg0 = Quantity(1.6, "eV")
    alpha = Quantity(0.1, "meV/K")
    beta = Quantity(600, "K")

    out = eg(T, eg0, alpha, beta)
    assert out.check("electron_volt")
    assert out < eg0


@mark.parametrize("fun", ["eg_gamma", "eg_x", "eg_l"])
def test_eg_at_point(fun):
    from pint import Quantity

    T = Quantity(298, "K")
    eg0 = Quantity(1.6, "eV")
    alpha = Quantity(0.1, "meV/K")
    beta = Quantity(600, "K")
    mock_eg = MagicMock()
    package = "solcore.parameter_sources.calculable_parameters"

    with patch(f"{package}.eg", mock_eg):
        from solcore.parameter_sources.calculable_parameters import eg_gamma, eg_x, eg_l

        f = {"eg_gamma": eg_gamma, "eg_x": eg_x, "eg_l": eg_l}[fun]
        f(T, eg0, alpha, beta)
        mock_eg.assert_called_once_with(T, eg0, alpha, beta)
        mock_eg.reset_mock()


def test_band_gap_and_lowest_band():
    from pint import Quantity
    from solcore.parameter_sources.calculable_parameters import band_gap, lowest_band

    eg_gamma = Quantity(1.4, "eV")
    eg_x = Quantity(1.6, "eV")
    eg_l = Quantity(1, "eV")

    with raises(ValueError):
        band_gap()

    gap = band_gap(eg_gamma, eg_x, eg_l)
    assert eg_l == gap
    assert Quantity("L", "dimensionless") == lowest_band(gap, eg_gamma, eg_x, eg_l)


def test_eff_mass_split_off():
    from pint import Quantity
    from solcore.parameter_sources.calculable_parameters import eff_mass_split_off

    g1 = Quantity(5.18)
    Ep = Quantity("18.7 eV")
    Delta_so = Quantity("0.676 eV")
    Eg = Quantity("1.42 eV")

    out = eff_mass_split_off(g1, Ep, Delta_so, Eg)
    assert out.u == "kilogram"


def test_eff_mass_z():
    from pint import Quantity
    from solcore.parameter_sources.calculable_parameters import (
        eff_mass_hh_z,
        eff_mass_lh_z,
    )
    from solcore.constants import electron_mass

    g1 = Quantity(5.18)
    g2 = Quantity(1.19)

    hh = eff_mass_hh_z(g1, g2)
    lh = eff_mass_lh_z(g1, g2)
    assert hh.u == "kilogram"
    assert lh.u == "kilogram"
    assert lh + hh == 2 * g1 * electron_mass
    assert lh - hh == 4 * g2 * electron_mass


def test_eff_mass_110():
    from pint import Quantity
    from solcore.parameter_sources.calculable_parameters import (
        eff_mass_hh_110,
        eff_mass_lh_110,
    )
    from solcore.constants import electron_mass

    g1 = Quantity(5.18)
    g2 = Quantity(1.19)
    g3 = Quantity(1.97)

    hh = eff_mass_hh_110(g1, g2, g3)
    lh = eff_mass_lh_110(g1, g2, g3)
    assert hh.u == "kilogram"
    assert lh.u == "kilogram"
    assert lh + hh == 2 * g1 * electron_mass
    assert round((lh - hh) / electron_mass, 6) == round(g2 + 3 * g3, 6)


def test_eff_mass_111():
    from pint import Quantity
    from solcore.parameter_sources.calculable_parameters import (
        eff_mass_hh_111,
        eff_mass_lh_111,
    )
    from solcore.constants import electron_mass

    g1 = Quantity(5.18)
    g3 = Quantity(1.97)

    hh = eff_mass_hh_111(g1, g3)
    lh = eff_mass_lh_111(g1, g3)
    assert hh.u == "kilogram"
    assert lh.u == "kilogram"
    assert lh + hh == 2 * g1 * electron_mass
    assert round(lh - hh, 6) == round(4 * g3 * electron_mass, 6)


def test_eff_mass_electron():
    from pint import Quantity
    from solcore.parameter_sources.calculable_parameters import eff_mass_electron
    from solcore.constants import electron_mass

    F = Quantity(-0.56)
    Ep = Quantity("18.7 eV")
    Delta_so = Quantity("0.676 eV")
    Eg = Quantity("1.42 eV")

    out = eff_mass_electron(F, Ep, Delta_so, Eg)
    assert out.u == "kilogram"

    out = eff_mass_electron(F, Ep * 0, Delta_so, Eg)
    assert out == (1 + 2 * F) * electron_mass


def test_permittivity():
    from pint import Quantity
    from solcore.parameter_sources.calculable_parameters import permittivity
    from solcore.constants import vacuum_permittivity

    er = Quantity(12.5)
    out = permittivity(er)
    assert out == vacuum_permittivity * 12.5


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
