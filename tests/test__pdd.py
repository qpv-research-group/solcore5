""" Poisson_drift_diffusion related tests
"""
from pytest import approx


def test_light_iv(AlGaAs, light_source):
    import numpy as np
    from solcore.solar_cell import SolarCell, default_GaAs
    from solcore import material
    from solcore.solar_cell_solver import solar_cell_solver

    expected = np.array(
        [
            142.68025180227374,
            2.519346556870366,
            0.9169672186977382,
            329.61395441947565,
            2.347826086956522,
            140.3911287342211,
            0.3294918264376029,
        ]
    )

    T = AlGaAs.T
    Vin = np.linspace(-2, 2.61, 201)
    V = np.linspace(0, 2.6, 300)

    substrate = material("GaAs")(T=T)

    my_solar_cell = SolarCell(
        [AlGaAs, default_GaAs(T)], T=T, R_series=0, substrate=substrate
    )
    solar_cell_solver(
        my_solar_cell,
        "iv",
        user_options={
            "T_ambient": T,
            "db_mode": "boltzmann",
            "voltages": V,
            "light_iv": True,
            "wavelength": light_source.x,
            "optics_method": "BL",
            "mpp": True,
            "internal_voltages": Vin,
            "light_source": light_source,
        },
    )

    output = [
        my_solar_cell.iv.Isc,
        my_solar_cell.iv.Voc,
        my_solar_cell.iv.FF,
        my_solar_cell.iv.Pmpp,
        my_solar_cell.iv.Vmpp,
        my_solar_cell.iv.Impp,
        my_solar_cell.iv.Eta,
    ]

    assert np.array(output) == approx(expected, rel=1e-2)


def test_quantum_efficiency(AlGaAs, light_source):
    import numpy as np
    from solcore.solar_cell import SolarCell, default_GaAs
    from solcore import material
    from solcore.solar_cell_solver import solar_cell_solver

    expected = np.array(
        [
            0.9866334968497021,
            2.1512408472022467e-14,
            0.9779769012349702,
            0.03506561338387434,
        ]
    )
    T = AlGaAs.T
    Vin = np.linspace(-2, 2.61, 201)
    V = np.linspace(0, 2.6, 300)

    substrate = material("GaAs")(T=T)

    my_solar_cell = SolarCell(
        [AlGaAs, default_GaAs(T)], T=T, R_series=0, substrate=substrate
    )

    solar_cell_solver(
        my_solar_cell,
        "qe",
        user_options={
            "T_ambient": T,
            "db_mode": "boltzmann",
            "voltages": V,
            "light_iv": True,
            "wavelength": light_source.x,
            "optics_method": "BL",
            "mpp": True,
            "internal_voltages": Vin,
            "light_source": light_source,
        },
    )

    output = [
        my_solar_cell[0].eqe(500e-9),
        my_solar_cell[0].eqe(800e-9),
        my_solar_cell[1].eqe(700e-9),
        my_solar_cell[1].eqe(900e-9),
    ]
    assert np.array(output) == approx(expected, abs=1e-3)
