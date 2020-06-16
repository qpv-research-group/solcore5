import sys
from pytest import mark

@mark.skipif(sys.platform != "linux", reason="Only works under linux")
def test_calculate_rat_rcwa():
    import numpy as np

    from solcore import material, si
    from solcore.structure import Layer, Structure
    from solcore.absorption_calculator.rigorous_coupled_wave import calculate_rat_rcwa

    T = 300
    Air = material("Air")(T=T)
    TiO2 = material("TiO2", sopra=True)(T=T)  # for the nanoparticles
    GaAs = material("GaAs")(T=T)

    NP_layer = Layer(
        si("50nm"),
        Air,
        geometry=[{"type": "circle", "mat": TiO2, "center": (200, 200), "radius": 50}],
    )
    np_struct = Structure([NP_layer])
    wl = np.linspace(300, 1000, 150)

    rat_np = calculate_rat_rcwa(
        np_struct, size=((400, 0), (0, 400)), orders=10, wavelength=wl, substrate=GaAs, incidence=Air
    )

    assert sorted(list(rat_np.keys())) == ["A", "A_per_layer", "R", "T"]
    for v in rat_np.values():
        assert len(v) == len(wl)


@mark.skipif(sys.platform != "linux", reason="Only works under linux")
def test_calculate_rat_rcwa():
    import numpy as np

    from solcore import material, si
    from solcore.structure import Layer, Structure
    from solcore.absorption_calculator.rigorous_coupled_wave import calculate_rat_rcwa

    T = 300
    Air = material("Air")(T=T)
    TiO2 = material("TiO2", sopra=True)(T=T)  # for the nanoparticles
    GaAs = material("GaAs")(T=T)

    NP_layer = Layer(
        si("50nm"),
        Air,
        geometry=[{"type": "circle", "mat": TiO2, "center": (200, 200), "radius": 50}],
    )
    np_struct = Structure([NP_layer])
    wl = np.linspace(300, 1000, 150)

    rat_np = calculate_rat_rcwa(
        np_struct, size=((400, 0), (0, 400)), orders=10, wavelength=wl, substrate=GaAs, incidence=Air
    )

    assert sorted(list(rat_np.keys())) == ["A", "A_per_layer", "R", "T"]
    for v in rat_np.values():
        assert len(v) == len(wl)


def test_calculate_absorption_profile_rcwa():
    import numpy as np

    from solcore import material, si
    from solcore.structure import Layer, Structure
    from solcore.absorption_calculator.rigorous_coupled_wave import calculate_rat_rcwa, calculate_absorption_profile_rcwa

    T = 300
    Air = material("Air")(T=T)
    TiO2 = material("TiO2", sopra=True)(T=T)  # for the nanoparticles
    GaAs = material("GaAs")(T=T)
    th = 50
    NP_layer = Layer(
        si(th, "nm"),
        Air,
        geometry=[{"type": "circle", "mat": TiO2, "center": (200, 200), "radius": 50}],
    )
    np_struct = Structure([NP_layer])
    wl = np.linspace(300, 1000, 10)

    rat_np = calculate_rat_rcwa(
        np_struct, size=((400, 0), (0, 400)), orders=10, wavelength=wl, substrate=GaAs, incidence=Air
    )

    A_output = rat_np['A']
    step_size = 2
    dist = np.arange(0, th, step_size)
    result = calculate_absorption_profile_rcwa(np_struct, size=((400, 0), (0, 400)), orders=10, wavelength=wl,
                                               rat_output_A=A_output, steps_size=step_size)

    parallel_result = calculate_absorption_profile_rcwa(np_struct, size=((400, 0), (0, 400)), orders=10, wavelength=wl,
                                               rat_output_A=A_output, parallel=True, steps_size=step_size)

    assert sorted(list(result.keys())) == ["absorption", "position"]
    assert sorted(list(parallel_result.keys())) == ["absorption", "position"]

    assert len(result["position"]) == len(dist)
    assert len(parallel_result["position"]) == len(dist)

    assert np.all(result["absorption"] == parallel_result["absorption"])