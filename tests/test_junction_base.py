from pytest import raises


def test_registry():
    from solcore.junction_base import JUNCTIONS_REGISTRY, JunctionBase

    class MyJunction(JunctionBase):
        pass

    assert MyJunction.__name__ in JUNCTIONS_REGISTRY


def test_not_overriden_behavior():
    from solcore.junction_base import JunctionBase

    class MyJunction(JunctionBase):
        def total_width(self):
            return None

        def widths(self):
            return None

        def solve_iv(self):
            return None

        def solve_qe(self):
            return None

    junc = MyJunction()
    with raises(NotImplementedError):
        junc.solve_equilibrium()

    with raises(NotImplementedError):
        junc.solve_short_circuit(None, None)

    assert junc.nk(None) is None
    assert junc.absorptivity(None, None) is None
