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
        np_struct, size=(400, 400), orders=10, wavelength=wl, substrate=GaAs
    )

    assert sorted(list(rat_np.keys())) == ["A", "R", "T"]
