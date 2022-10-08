from pytest import raises, mark, approx
from pint import Quantity
from hypothesis import given
import hypothesis.strategies as st


def test_mobility_low_field():
    from solcore.future.parameter_sources.mobility_parameters import mobility_low_field

    T = Quantity(298, "K")
    N = Quantity(1e18, "1/cm**3")
    mu_min = Quantity(500, "cm**2/V/s")
    mu_max = Quantity(9200, "cm**2/V/s")
    Nref = Quantity(6e16, "1/cm**3")
    ll = 0.394
    t1 = 1.1
    t2 = 3.0

    out = mobility_low_field(T, N, mu_min, mu_max, Nref, ll, t1, t2)
    out_larger = mobility_low_field(T, N / 10, mu_min, mu_max, Nref, ll, t1, t2)
    assert out_larger > out

    smallT = mobility_low_field(T * 0.01, N, mu_min, mu_max, Nref, ll, t1, t2)
    largeT = mobility_low_field(T * 100, N, mu_min, mu_max, Nref, ll, t1, t2)
    assert smallT < out
    assert largeT < out

    t1 = 2.1
    smallT = mobility_low_field(T * 0.01, N, mu_min, mu_max, Nref, ll, t1, t2)
    largeT = mobility_low_field(T * 100, N, mu_min, mu_max, Nref, ll, t1, t2)
    assert out < smallT
    assert largeT < out


class TestSotoodehMobilitySource:
    def test_load_source(self):
        from solcore.future.parameter_sources import SotoodehMobilitySource

        sm = SotoodehMobilitySource.load_source()
        assert sm == SotoodehMobilitySource()
        assert len(sm.reference) > 0
        assert len(sm._descriptions) > 0
        assert len(sm._data) > 0

    def test_materials(self):
        from solcore.future.parameter_sources import SotoodehMobilitySource

        sm = SotoodehMobilitySource()
        assert sm.materials == tuple((m for m in sm._data if len(m) <= 4))

    def test_parameters(self):
        from solcore.future.parameter_sources import SotoodehMobilitySource

        sm = SotoodehMobilitySource()

        assert sm.parameters("Dark matter") == ()
        assert sm.parameters(sm.materials[0]) == ("electron_mobility", "hole_mobility")

    @mark.parametrize("carrier", ["electron", "hole"])
    @mark.parametrize(
        "kwargs",
        [
            {},
            {"T": Quantity(298, "K")},
            {"T": Quantity(298, "K"), "Nd": Quantity(1e18, "1/cm**3")},
            {"T": Quantity(298, "K"), "Na": Quantity(1e18, "1/cm**3")},
        ],
    )
    def test_get_parameter(self, carrier, kwargs):
        from solcore.future.parameter_sources import SotoodehMobilitySource
        from solcore.future.parameter import InputArgumentMissing

        sm = SotoodehMobilitySource()
        material = sm.materials[0]

        if "T" not in kwargs or ("Nd" not in kwargs and "Na" not in kwargs):
            with raises(InputArgumentMissing):
                sm.get_parameter(material, f"{carrier}_mobility", **kwargs)
        else:
            out = sm.get_parameter(material, f"{carrier}_mobility", **kwargs)
            assert out.d == sm._descriptions[f"{carrier}_mobility"]
            assert out.r == (sm.name,)
            assert out.u == "meter ** 2 / second / volt"

    @mark.xfail(reason="Method not implemented, yet.")
    def test__get_parameter_alloy(self):
        assert False


@given(
    egg=st.floats(0.2, 2),
    egx=st.floats(0.2, 2),
    egl=st.floats(0.2, 2),
    mg=st.floats(0.001, 3),
    mx=st.floats(0.001, 3),
    ml=st.floats(0.001, 3),
    T=st.floats(20, 500),
)
def test_fraction_electrons_direct_valley(egg, egx, egl, mg, mx, ml, T):
    from solcore.future.parameter_sources.mobility_parameters import (
        fraction_electrons_direct_valley,
    )

    egg = Quantity(egg, "eV")
    egx = Quantity(egx, "eV")
    egl = Quantity(egl, "eV")
    mg = Quantity(mg)
    mx = Quantity(mx)
    ml = Quantity(ml)
    T = Quantity(T, "K")

    fraction = fraction_electrons_direct_valley(egg, egx, egl, mg, mx, ml, T)
    assert 0 <= fraction <= 1

    egg = 0.9 * egg
    new_fraction = fraction_electrons_direct_valley(egg, egx, egl, mg, mx, ml, T)
    assert fraction <= new_fraction


@given(
    egx=st.floats(0.2, 2),
    egl=st.floats(0.2, 2),
    mx=st.floats(0.001, 3),
    ml=st.floats(0.001, 3),
    T=st.floats(20, 500),
)
def test_effective_indirect_mass(egx, egl, mx, ml, T):
    from solcore.future.parameter_sources.mobility_parameters import (
        effective_indirect_mass,
    )

    egx = Quantity(egx, "eV")
    egl = Quantity(egl, "eV")
    mx = Quantity(mx)
    ml = Quantity(ml)
    T = Quantity(T, "K")

    m = effective_indirect_mass(egx, egl, mx, ml, T)
    assert abs(m - mx).m + abs(m - ml).m == approx(abs(mx - ml).m)

    egl = egl / 2
    new_m = effective_indirect_mass(egx, egl, mx, ml, T)
    if not mx.m == approx(ml.m):
        assert abs(m - ml) >= abs(new_m - ml)
    else:
        assert m.m == approx(new_m.m)


@given(e1=st.floats(4, 20), e2=st.floats(4, 20), x=st.floats(0, 1))
def test_interpolate_epsilon(e1, e2, x):
    from solcore.future.parameter_sources.mobility_parameters import interpolate_epsilon

    ratio1 = (e1 - 1) / (e1 + 2)
    ratio2 = (e2 - 1) / (e2 + 2)
    expected = x * ratio1 + (1 - x) * ratio2
    eo = interpolate_epsilon(e1, e2, x)
    actual = (eo - 1) / (eo + 2)
    assert actual == approx(expected)


@given(n1=st.floats(1e12, 1e22), n2=st.floats(1e12, 1e22), x=st.floats(0, 1))
def test_interpolate_concentration(n1, n2, x):
    from solcore.future.parameter_sources.mobility_parameters import (
        interpolate_concentration,
    )
    import numpy as np

    expected = x * np.log10(n1) + (1.0 - x) * np.log10(n2)
    assert np.log10(interpolate_concentration(n1, n2, x)) == approx(expected)
