import sys
from pytest import mark, approx

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

    assert sorted(list(rat_np.keys())) == ["A", "A_per_layer", "A_pol", "R", "T"]
    for v in rat_np.values():
        assert v.shape[0] == len(wl)


@mark.skipif(sys.platform != "linux", reason="Only works under linux")
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
        geometry=[{"type": "ellipse", "mat": TiO2, "center": (200, 200), "halfwidths": [70, 100], "angle": 30}],
    )
    np_struct = Structure([NP_layer])
    wl = np.linspace(300, 1000, 10)

    rat_np = calculate_rat_rcwa(
        np_struct, size=((400, 0), (0, 400)), orders=10, wavelength=wl, substrate=GaAs, incidence=Air
    )

    A_output = rat_np['A_pol']
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

    assert approx(result["absorption"] == parallel_result["absorption"])


@mark.skipif(sys.platform != "linux", reason="Only works under linux")
def test_pol_rcwa():
    import numpy as np

    from solcore import material, si
    from solcore.structure import Layer, Structure
    from solcore.absorption_calculator.rigorous_coupled_wave import calculate_rat_rcwa, \
        calculate_absorption_profile_rcwa

    T = 300
    Air = material("Air")(T=T)
    TiO2 = material("TiO2", sopra=True)(T=T)  # for the nanoparticles
    GaAs = material("GaAs")(T=T)
    th = 50
    NP_layer = Layer(
        si(th, "nm"),
        Air,
        geometry=[{"type": "rectangle", "mat": TiO2, "center": (200, 200), "halfwidths": [100, 70], "angle": 40}],
    )
    np_struct = Structure([NP_layer])
    wl = np.linspace(300, 1000, 10)

    step_size = 2

    rat_np_u = calculate_rat_rcwa(
        np_struct, size=((400, 0), (0, 400)), orders=10, wavelength=wl, substrate=GaAs, incidence=Air, pol='u'
    )

    A_output_u = rat_np_u['A_pol']
    result_u = calculate_absorption_profile_rcwa(np_struct, size=((400, 0), (0, 400)), orders=10, wavelength=wl,
                                                 rat_output_A=A_output_u, parallel=True, steps_size=step_size,
                                                 pol='u')

    rat_np_s = calculate_rat_rcwa(
        np_struct, size=((400, 0), (0, 400)), orders=10, wavelength=wl, substrate=GaAs, incidence=Air, pol='s'
    )

    A_output_s = rat_np_s['A_pol']
    result_s = calculate_absorption_profile_rcwa(np_struct, size=((400, 0), (0, 400)), orders=10, wavelength=wl,
                                                 rat_output_A=A_output_s, parallel=True, steps_size=step_size,
                                                 pol='s')

    rat_np_p = calculate_rat_rcwa(
        np_struct, size=((400, 0), (0, 400)), orders=10, wavelength=wl, substrate=GaAs, incidence=Air, pol='p'
    )

    A_output_p = rat_np_p['A_pol']
    result_p = calculate_absorption_profile_rcwa(np_struct, size=((400, 0), (0, 400)), orders=10, wavelength=wl,
                                                 rat_output_A=A_output_p, parallel=True, steps_size=step_size,
                                                 pol='p')

    assert approx(A_output_u == 0.5*(A_output_p + A_output_s))

    assert approx(result_u["absorption"] == 0.5*(result_s["absorption"] + result_p["absorption"]))


@mark.skipif(sys.platform != "linux", reason="Only works under linux")
def test_arbitrary_pol_rcwa():
    import numpy as np

    from solcore import material, si
    from solcore.structure import Layer, Structure
    from solcore.absorption_calculator.rigorous_coupled_wave import calculate_rat_rcwa, \
        calculate_absorption_profile_rcwa

    T = 300
    Air = material("Air")(T=T)
    TiO2 = material("TiO2", sopra=True)(T=T)  # for the nanoparticles
    GaAs = material("GaAs")(T=T)
    th = 50
    NP_layer = Layer(
        si(th, "nm"),
        Air,
        geometry=[{"type": "ellipse", "mat": TiO2, "center": (200, 200), "halfwidths": [100, 70], "angle": 40}],
    )
    np_struct = Structure([NP_layer])
    wl = np.linspace(300, 1000, 10)

    step_size = 2

    rat_np = calculate_rat_rcwa(
        np_struct, size=((400, 0), (0, 400)), orders=10, wavelength=wl, substrate=GaAs, incidence=Air,
        pol=[1, 0]
    )

    A_output = rat_np['A_pol']
    result = calculate_absorption_profile_rcwa(np_struct, size=((400, 0), (0, 400)), orders=10, wavelength=wl,
                                                 rat_output_A=A_output, parallel=True, steps_size=step_size,
                                                 pol=[1, 0])

    rat_np_s = calculate_rat_rcwa(
        np_struct, size=((400, 0), (0, 400)), orders=10, wavelength=wl, substrate=GaAs, incidence=Air, pol='s'
    )

    A_output_s = rat_np_s['A_pol']
    result_s = calculate_absorption_profile_rcwa(np_struct, size=((400, 0), (0, 400)), orders=10, wavelength=wl,
                                                 rat_output_A=A_output_s, parallel=True, steps_size=step_size,
                                                 pol='s')

    assert approx(A_output == A_output_s)

    assert approx(result["absorption"] == result_s["absorption"])


@mark.skipif(sys.platform != "linux", reason="Only works under linux")
def test_rcwa_polygon():
    import numpy as np

    from solcore import material, si
    from solcore.structure import Layer, Structure
    from solcore.absorption_calculator.rigorous_coupled_wave import calculate_rat_rcwa, \
        calculate_absorption_profile_rcwa

    T = 300
    Air = material("Air")(T=T)
    TiO2 = material("TiO2", sopra=True)(T=T)  # for the nanoparticles
    GaAs = material("GaAs")(T=T)
    th = 50

    NP_layer_square = Layer(
        si(th, "nm"),
        Air,
        geometry=[{"type": "rectangle", "mat": TiO2, "center": (200, 200), "halfwidths": [50, 50], "angle": 0}],
    )
    np_struct_square = Structure([NP_layer_square])

    NP_layer_polygon = Layer(
        si(th, "nm"),
        Air,
        geometry=[{"type": "polygon", "mat": TiO2, "center": (200, 200),
                   "vertices": ((150, 150), (250,150), (250, 250), (150, 250)),
                                "angle": 0}],
    )
    np_struct_polygon = Structure([NP_layer_polygon])

    wl = np.linspace(300, 1000, 10)

    step_size = 2

    rat_np_sq = calculate_rat_rcwa(
        np_struct_square, size=((400, 0), (0, 400)), orders=10, wavelength=wl, substrate=GaAs, incidence=Air,
        pol=[1, 0]
    )

    A_output_sq = rat_np_sq['A_pol']
    result_square = calculate_absorption_profile_rcwa(np_struct_square, size=((400, 0), (0, 400)), orders=10, wavelength=wl,
                                                 rat_output_A=A_output_sq, parallel=True, steps_size=step_size,
                                                 pol=[1, 0])

    rat_np_pol = calculate_rat_rcwa(
        np_struct_polygon, size=((400, 0), (0, 400)), orders=10, wavelength=wl, substrate=GaAs, incidence=Air,
        pol=[1, 0]
    )

    A_output_pol = rat_np_pol['A_pol']
    result_polygon = calculate_absorption_profile_rcwa(np_struct_polygon, size=((400, 0), (0, 400)), orders=10, wavelength=wl,
                                                 rat_output_A=A_output_pol, parallel=True, steps_size=step_size,
                                                 pol=[1, 0])



    assert approx(A_output_sq == A_output_pol)

    assert approx(result_polygon["absorption"] == result_square["absorption"])