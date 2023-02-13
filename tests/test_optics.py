from .conftest import skip_s4_test
from pytest import mark, approx


@mark.skipif(skip_s4_test(), reason="Only works if S4 installed")
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


@mark.skipif(skip_s4_test(), reason="Only works if S4 installed")
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

    assert result["absorption"] == approx(parallel_result["absorption"], abs=1e-5)


@mark.skipif(skip_s4_test(), reason="Only works if S4 installed")
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

    assert np.mean(A_output_u, 1) == approx(0.5*(A_output_p + A_output_s))

    assert result_u["absorption"] == approx(0.5*(result_s["absorption"] + result_p["absorption"]))


@mark.skipif(skip_s4_test(), reason="Only works if S4 installed")
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

    assert A_output == approx(A_output_s)

    assert result["absorption"] == approx(result_s["absorption"])


@mark.skipif(skip_s4_test(), reason="Only works if S4 installed")
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



    assert A_output_sq == approx(A_output_pol)

    assert result_polygon["absorption"] == approx(result_square["absorption"], abs=0.001)


@mark.skipif(skip_s4_test(), reason="Only works if S4 installed")
@mark.parametrize("pol", ['s', 'p', 'u'])
@mark.parametrize("angle", [0, 60])
def test_tmm_rcwa_structure_comparison(pol, angle):
    import numpy as np
    from solcore import si, material
    from solcore.structure import Layer
    from solcore.solar_cell import SolarCell
    from solcore.optics.tmm import calculate_rat, OptiStack
    from solcore.optics.rcwa import calculate_rat_rcwa

    InGaP = material('GaInP')(In=0.5)
    GaAs = material('GaAs')()
    Ge = material('Ge')()
    Ag = material('Ag')()
    Air = material('Air')()

    Al2O3 = material('Al2O3')()

    wavelengths = np.linspace(250, 1900, 500)

    size = ((100, 0), (0, 100))

    ARC = [Layer(si('80nm'), Al2O3)]

    solar_cell = SolarCell(ARC + [Layer(material=InGaP, width=si('400nm')),
                                  Layer(material=GaAs, width=si('4000nm')),
                                  Layer(material=Ge, width=si('3000nm'))], substrate=Ag)

    solar_cell_OS = OptiStack(solar_cell, no_back_reflection=False, substrate=solar_cell.substrate)


    rcwa_result = calculate_rat_rcwa(solar_cell, size, 2, wavelengths, Air, Ag, theta=angle, phi=angle,
                                     pol=pol, parallel=True)
    tmm_result = calculate_rat(solar_cell_OS, wavelengths, angle, pol, no_back_reflection=False)

    assert tmm_result['A_per_layer'][1:-1] == approx(rcwa_result['A_per_layer'].T)
    assert tmm_result['R'] == approx(rcwa_result['R'])
    assert tmm_result['T'] == approx(rcwa_result['T'])

    assert np.sum(tmm_result['A_per_layer'][1:-1].T, 1) + tmm_result['R'] + tmm_result['T'] == approx(1)
    assert np.sum(rcwa_result['A_per_layer'], 1) + rcwa_result['R'] + rcwa_result['T'] == approx(1)


@mark.skipif(skip_s4_test(), reason="Only works if S4 installed")
@mark.parametrize("pol", ['s', 'p', 'u'])
@mark.parametrize("angle", [0, 60])
def test_tmm_rcwa_structure_profile_comparison(pol, angle):
    import numpy as np
    from solcore import si, material
    from solcore.structure import Layer
    from solcore.solar_cell import SolarCell
    from solcore.optics.tmm import calculate_rat, calculate_absorption_profile, OptiStack
    from solcore.optics.rcwa import calculate_rat_rcwa, calculate_absorption_profile_rcwa

    InGaP = material('GaInP')(In=0.5)
    GaAs = material('GaAs')()
    Ge = material('Ge')()
    Ag = material('Ag')()
    Air = material('Air')()

    wavelengths = np.linspace(250, 1900, 8)

    size = ((100, 0), (0, 100))

    solar_cell = SolarCell([Layer(material=InGaP, width=si('400.5nm')),
                                  Layer(material=GaAs, width=si('500nm')),
                                  Layer(material=Ge, width=si('1000nm'))], substrate=Ag)

    solar_cell_OS = OptiStack(solar_cell, no_back_reflection=False, substrate=solar_cell.substrate)

    rcwa_result = calculate_rat_rcwa(solar_cell, size, 2, wavelengths, Air, Ag, theta=angle, phi=angle,
                                     pol=pol, parallel=True)
    tmm_result = calculate_rat(solar_cell_OS, wavelengths, angle, pol, no_back_reflection=False)

    tmm_profile = calculate_absorption_profile(solar_cell_OS, wavelengths, no_back_reflection=False,
                                               angle=angle, pol=pol, RAT_out=tmm_result, steps_size=2)

    rcwa_profile = calculate_absorption_profile_rcwa(solar_cell, size, 2, wavelengths, rcwa_result['A_pol'],
                                                     theta=angle, phi=angle, pol=pol, incidence=Air,
                                                     substrate=solar_cell.substrate, parallel=True, steps_size=2)


    assert tmm_profile['position'] == approx(rcwa_profile['position'], abs=2e-6)
    assert tmm_profile['absorption'] == approx(rcwa_profile['absorption'], abs=2e-6)