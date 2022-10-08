from pytest import approx


def test_registry():
    from solcore.future.junction_base import JUNCTIONS_REGISTRY, JunctionBase

    class MyJunction(JunctionBase):
        pass

    assert MyJunction.__name__ in JUNCTIONS_REGISTRY


def test_not_overwritten_behavior():
    from solcore.future.junction_base import JunctionBase

    class MyJunction(JunctionBase):
        def solve_equilibrium(self):
            return None

        def solve_short_circuit(self):
            return None

        def solve_iv(self):
            return None

        def solve_qe(self):
            return None

    junc = MyJunction()

    assert junc.total_width is None
    assert junc.widths is None
    assert junc.nk(None) is None
    assert junc.absorptivity(None, None) is None


def test_iv_parameters():
    from solcore.future.junction_base import iv_parameters, JunctionBase
    import numpy as np

    s = 1.5
    voltage = np.arange(-1, 2, 0.05)
    current = s - voltage

    expected = dict(
        Jsc=s, Voc=s, Pmmp=(s / 2) ** 2, Jmpp=(s / 2), Vmpp=(s / 2), FF=1 / 4
    )
    actual1 = iv_parameters(voltage, current)
    actual2 = JunctionBase.iv_parameters(voltage, current)
    for k, v in expected.items():
        assert actual1[k] == approx(v, abs=0.025)
        assert actual2[k] == approx(v, abs=0.025)
