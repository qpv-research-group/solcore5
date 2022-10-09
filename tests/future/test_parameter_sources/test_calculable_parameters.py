from pytest import raises, mark, approx
from unittest.mock import patch, MagicMock


def test_eg():
    from pint import Quantity
    from solcore.future.parameter_sources.calculable_parameters import eg

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
    package = "solcore.future.parameter_sources.calculable_parameters"

    with patch(f"{package}.eg", mock_eg):
        from solcore.future.parameter_sources.calculable_parameters import (
            eg_gamma,
            eg_x,
            eg_l,
        )

        f = {"eg_gamma": eg_gamma, "eg_x": eg_x, "eg_l": eg_l}[fun]
        f(T, eg0, alpha, beta)
        mock_eg.assert_called_once_with(T, eg0, alpha, beta)
        mock_eg.reset_mock()


def test_band_gap_and_lowest_band():
    from pint import Quantity
    from solcore.future.parameter_sources.calculable_parameters import (
        band_gap,
        lowest_band,
    )

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
    from solcore.future.parameter_sources.calculable_parameters import (
        eff_mass_split_off,
    )

    g1 = Quantity(5.18)
    Ep = Quantity("18.7 eV")
    Delta_so = Quantity("0.676 eV")
    Eg = Quantity("1.42 eV")

    out = eff_mass_split_off(g1, Ep, Delta_so, Eg)
    assert out.u == "kilogram"


def test_eff_mass_z():
    from pint import Quantity
    from solcore.future.parameter_sources.calculable_parameters import (
        eff_mass_hh_z,
        eff_mass_lh_z,
    )
    from solcore.future.constants import electron_mass_

    g1 = Quantity(5.18)
    g2 = Quantity(1.19)

    hh = eff_mass_hh_z(g1, g2)
    lh = eff_mass_lh_z(g1, g2)
    assert hh.u == "kilogram"
    assert lh.u == "kilogram"
    assert (electron_mass_ / lh + electron_mass_ / hh).m == approx(2 * g1.m)
    assert (electron_mass_ / lh - electron_mass_ / hh).m == approx(4 * g2.m)


def test_eff_mass_110():
    from pint import Quantity
    from solcore.future.parameter_sources.calculable_parameters import (
        eff_mass_hh_110,
        eff_mass_lh_110,
    )
    from solcore.future.constants import electron_mass_

    g1 = Quantity(5.18)
    g2 = Quantity(1.19)
    g3 = Quantity(1.97)

    hh = eff_mass_hh_110(g1, g2, g3)
    lh = eff_mass_lh_110(g1, g2, g3)
    assert hh.u == "kilogram"
    assert lh.u == "kilogram"
    assert (electron_mass_ / lh + electron_mass_ / hh).m == approx(2 * g1.m)
    assert (electron_mass_ / lh - electron_mass_ / hh).m == approx(g2.m + 3 * g3.m)


def test_eff_mass_111():
    from pint import Quantity
    from solcore.future.parameter_sources.calculable_parameters import (
        eff_mass_hh_111,
        eff_mass_lh_111,
    )
    from solcore.future.constants import electron_mass_

    g1 = Quantity(5.18)
    g3 = Quantity(1.97)

    hh = eff_mass_hh_111(g1, g3)
    lh = eff_mass_lh_111(g1, g3)
    assert hh.u == "kilogram"
    assert lh.u == "kilogram"
    assert electron_mass_ / lh + electron_mass_ / hh == approx(2 * g1.m)
    assert (electron_mass_ / lh - electron_mass_ / hh).m == approx(4 * g3.m)


def test_eff_mass_electron():
    from pint import Quantity
    from solcore.future.parameter_sources.calculable_parameters import eff_mass_electron
    from solcore.future.constants import electron_mass_

    F = Quantity(-0.56)
    Ep = Quantity("18.7 eV")
    Delta_so = Quantity("0.676 eV")
    Eg = Quantity("1.42 eV")

    out = eff_mass_electron(F, Ep, Delta_so, Eg)
    assert out.u == "kilogram"

    out = eff_mass_electron(F, Ep * 0, Delta_so, Eg)
    assert electron_mass_ / out == (1 + 2 * F)


def test_permittivity():
    from pint import Quantity
    from solcore.future.parameter_sources.calculable_parameters import permittivity
    from solcore.future.constants import vacuum_permittivity_

    er = Quantity(12.5)
    out = permittivity(er)
    assert out == vacuum_permittivity_ * 12.5


class TestCalculableParameters:
    def test__getitem__(self):
        from solcore.future.parameter_sources import CalculableParameters
        from solcore.future.parameter import ParameterMissing

        cp = CalculableParameters()
        with raises(ParameterMissing):
            cp["param1"]

        assert cp._params["eg_gamma"] == cp["eg_gamma"]

    def test___getattr__(self):
        from solcore.future.parameter_sources import CalculableParameters
        from solcore.future.parameter import ParameterMissing

        cp = CalculableParameters()
        with raises(ParameterMissing):
            cp.param1

        assert cp._params["eg_gamma"] == cp.eg_gamma

    def test_load_source(self):
        from solcore.future.parameter_sources import CalculableParameters

        cp = CalculableParameters()
        assert cp == cp.load_source()

    def test_register_calculable(self):
        from solcore.future.parameter_sources import CalculableParameters
        from solcore.future.parameter import ParameterSourceError

        cp = CalculableParameters()
        assert len(cp._params) > 0
        assert tuple(cp._params.keys()) == tuple(cp._descriptions.keys())

        with raises(ParameterSourceError):

            @CalculableParameters.register_calculable
            def eg_gamma():
                pass

    def test_materials(self):
        from solcore.future.parameter_sources import CalculableParameters

        cp = CalculableParameters()
        cp._warned = False
        assert cp.materials == ()
        assert cp._warned

    def test_parameters(self):
        from solcore.future.parameter_sources import CalculableParameters

        cp = CalculableParameters()
        assert cp.parameters() == tuple((p for p in cp._params))

    def test_list_arguments(self):
        from solcore.future.parameter_sources import CalculableParameters

        cp = CalculableParameters()
        expected = ("T", "eg0_gamma", "alpha_gamma", "beta_gamma")
        assert cp.list_arguments("eg_gamma") == expected

    def test_get_parameter(self):
        from solcore.future.parameter_sources import CalculableParameters
        from solcore.future.parameter import (
            InputArgumentMissing,
            Parameter,
            ParameterMissing,
            ParameterManager,
        )

        ParameterManager().add_source("Calculable", CalculableParameters)
        ParameterManager().initialize()
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
                raise ParameterMissing("A source", "gold", "value")

        cp.parman.get_parameter = get_param

        @CalculableParameters.register_calculable(description="Answer is 16")
        def dummy_param_16(param_16):
            return param_16

        out = cp.get_parameter("Dark matter", "dummy_param_16")
        assert out.m == 16
        assert out.d == "Answer is 16"
        assert set(out.r) == {"Calculable", "SomePaper"}

        with raises(ParameterMissing):
            cp.get_parameter("Dark matter", "other_dummy_param")

        @CalculableParameters.register_calculable(description="Answer is 16")
        def dummy_param_with_T(param_16, T):
            return param_16 * T

        with raises(InputArgumentMissing):
            cp.get_parameter("Dark matter", "dummy_param_with_T")

        # Let's get rid of all the dummies...
        dummies = (d for d in cp.parameters() if "dummy" in d)
        for d in dummies:
            cp._params.pop(d)


def test_electron_affinity():
    from pint import Quantity
    from solcore.future.parameter_sources.calculable_parameters import electron_affinity

    valence_band_offset = Quantity("-0.8 eV")
    Eg = Quantity("1.42 eV")
    electron_affinity_InSb = Quantity("0 eV")
    Eg_InSb = Quantity("0 eV")

    x = electron_affinity(valence_band_offset, Eg, Eg_InSb, electron_affinity_InSb)
    assert x == -valence_band_offset - Eg


def test_density_states():
    from pint import Quantity
    from solcore.future.parameter_sources.calculable_parameters import density_states
    from solcore.future.constants import electron_mass_

    T = Quantity(298, "K")
    mass = 0.1 * electron_mass_

    ds = density_states(T, mass)
    assert density_states(1.1 * T, mass).m == approx(ds.m * 1.1 ** (3 / 2))
    assert density_states(0.9 * T, mass).m == approx(ds.m * 0.9 ** (3 / 2))


def test_nc():
    from pint import Quantity
    from solcore.future.constants import electron_mass_

    T = Quantity(298, "K")
    mass = 0.1 * electron_mass_
    mock_density_states = MagicMock()
    package = "solcore.future.parameter_sources.calculable_parameters"

    with patch(f"{package}.density_states", mock_density_states):
        from solcore.future.parameter_sources.calculable_parameters import Nc

        Nc(T, mass)
        assert mock_density_states.call_count == 1


def test_nv():
    from pint import Quantity
    from solcore.future.constants import electron_mass_

    T = Quantity(298, "K")
    mass = 0.1 * electron_mass_
    mock_density_states = MagicMock()
    package = "solcore.future.parameter_sources.calculable_parameters"

    with patch(f"{package}.density_states", mock_density_states):
        from solcore.future.parameter_sources.calculable_parameters import Nv

        Nv(T, 2 * mass, mass)
        assert mock_density_states.call_count == 2


def test_ni():
    from pint import Quantity
    from solcore.future.parameter_sources.calculable_parameters import ni

    T = Quantity(298, "K")
    Nc = Quantity(1e18, "1/cm**3")
    Nv = Quantity(1e17, "1/cm**3")
    band_gap = Quantity("1.42 eV")

    ni0 = ni(T, Nc, Nv, band_gap)
    assert ni(0.9 * T, Nc, Nv, band_gap) < ni0
    assert ni(1.1 * T, Nc, Nv, band_gap) > ni0
    assert ni(T, Nc, Nv, 0.9 * band_gap) > ni0
    assert ni(T, Nc, Nv, 1.1 * band_gap) < ni0
