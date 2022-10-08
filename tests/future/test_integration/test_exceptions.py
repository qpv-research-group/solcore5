from hypothesis import given
from hypothesis.strategies import sampled_from


def gather_known_sources():
    from solcore.future.parameter import ParameterManager

    pm = ParameterManager()
    pm.initialize()

    return tuple(pm.sources.keys())


def gather_known_materials():
    from solcore.future.parameter import ParameterManager

    pm = ParameterManager()
    pm.initialize()

    materials = set()
    for source in gather_known_sources():
        materials = materials | set(pm.sources[source].materials)

    return tuple(materials)


def gather_known_parameters():
    from solcore.future.parameter import ParameterManager, MaterialMissing

    pm = ParameterManager()
    pm.initialize()

    parameters = set()
    for m in gather_known_materials():
        for source in pm.sources.values():
            try:
                parameters = parameters | set(source.parameters(m))
            except MaterialMissing:
                pass

    exclude = ("parent0", "parent1", "x")
    return tuple((p for p in parameters if p not in exclude))


def get_parameter(material, parameter, attempt=0, **kwargs):
    from solcore.future.parameter import (
        InputArgumentMissing,
        ParameterMissing,
        MaterialMissing,
        ParameterManager,
    )
    from pint import Quantity

    if attempt > 4:
        raise RecursionError(
            f"Too many attempts to retrieve '{parameter}' for '{material}'"
        )

    pm = ParameterManager()
    try:
        pm.get_parameter(material, parameter, **kwargs)

    except MaterialMissing as err:
        assert (
            pm.sources[err.source].materials == ()
            or err.material not in pm.sources[err.source].materials
        )

    except ParameterMissing as err:
        if isinstance(err.source, str):
            assert err.parameter not in pm.sources[err.source].parameters(err.material)
        else:
            for s in err.source:
                assert err.parameter not in pm.sources[s].parameters(err.material)

    except InputArgumentMissing as err:
        a = err.argument
        if a == "T":
            assert a not in kwargs
            kwargs[a] = Quantity(300, "K")
        elif a == "Nd or Na":
            assert "Nd" not in kwargs and "Na" not in kwargs
            kwargs["Nd"] = Quantity(1e17, "1/cm**3")
        else:
            assert a not in kwargs.get("comp", {})
            if "comp" not in kwargs:
                kwargs["comp"] = {a: 0.2}
            else:
                kwargs["comp"][a] = 0.2

        attempt += 1
        get_parameter(material, parameter, attempt, **kwargs)


@given(
    material=sampled_from(gather_known_materials()),
    parameter=sampled_from(gather_known_parameters()),
)
def test_parameter_exceptions(material, parameter):
    print(material, parameter)
    get_parameter(material, parameter)
