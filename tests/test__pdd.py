""" Poisson_drift_diffusion related tests
"""
from pytest import approx, mark


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


def test_calculate_iv(mocker):
    from . import mock_ddModel

    root = "solcore.poisson_drift_diffusion.DriftDiffusionUtilities"
    mocker.patch(f"{root}.dd", new=mock_ddModel)
    mock_sc = mocker.patch(f"{root}.short_circuit_pdd")
    mock_eq = mocker.patch(f"{root}.equilibrium_pdd")
    mock_bandstructure = mocker.patch(f"{root}.DumpBandStructure")
    mock_dump_iv = mocker.patch(f"{root}.DumpIV")

    from solcore.structure import Junction
    from solcore.poisson_drift_diffusion.DriftDiffusionUtilities import calculate_iv

    junction = Junction()
    vlimit = 1
    vstep = 0.2
    output_iv = 1

    # We check all relevant functions are called when running in the dark.
    out = calculate_iv(junction, vlimit, vstep, light_iv=False, output_iv=output_iv)
    mock_eq.assert_called()
    mock_sc.assert_not_called()
    mock_ddModel.runiv.assert_called_with(vlimit, vstep, output_iv, 0)
    mock_bandstructure.assert_called()
    mock_dump_iv.assert_called()
    assert "Bandstructure" in out and "IV" in out

    # And with illumination
    mock_eq.reset_mock()
    mock_sc.reset_mock()
    out = calculate_iv(junction, vlimit, vstep, light_iv=True, output_iv=output_iv)
    mock_eq.assert_not_called()
    mock_sc.assert_called()

    # And if the vlimit is negative, so it will be the step, internally
    out = calculate_iv(junction, -vlimit, vstep, light_iv=True, output_iv=output_iv)
    mock_ddModel.runiv.assert_called_with(-vlimit, -vstep, output_iv, 0)


def test_find_minimum_bandgap():
    from solcore.structure import Layer, Junction
    from solcore import material
    from solcore.constants import q
    from solcore.poisson_drift_diffusion.DriftDiffusionUtilities import (
        find_minimum_bandgap,
    )

    InGaAs = material("InGaAs")
    junction = Junction([Layer(1, InGaAs(In=x)) for x in (0, 0.1, 0.2)])
    bandgap = find_minimum_bandgap(junction)
    assert bandgap * q == junction[-1].material.band_gap


def kbT(T: float = 300) -> float:
    from solcore.constants import q, kb

    return 3 * kb * T / q


@mark.parametrize(
    ["bandgap", "p_on_n", "vmax", "vmin", "exp_vmax", "exp_vmin"],
    [
        [1, True, 1.5, -1.5, 1 + kbT(), -1.5],
        [1, False, 1.5, -1.5, 1.5, -1 - kbT()],
        [1.6, True, 1.5, -1.5, 1.5, -1.5],
        [1.6, False, 1.5, -1.5, 1.5, -1.5],
    ],
)
def test_find_voltage_limits(bandgap, p_on_n, vmax, vmin, exp_vmax, exp_vmin):
    from solcore.poisson_drift_diffusion.DriftDiffusionUtilities import (
        find_voltage_limits,
    )

    act_vmax, act_vmin = find_voltage_limits(bandgap, vmax, vmin, p_on_n, T=300)
    assert act_vmax == exp_vmax
    assert act_vmin == exp_vmin


def test_consolidate_iv():
    from solcore.poisson_drift_diffusion.DriftDiffusionUtilities import (
        consolidate_iv,
    )

    import numpy as np

    vmax = 5.1
    vmin = -5.1
    positive = {
        "V": np.arange(0, vmax, 1),
        "I": np.arange(0, vmax, 1),
    }
    negative = {
        "V": np.arange(0, vmin, -1),
        "I": np.arange(0, vmin, -1),
    }
    total = {
        "V": np.arange(min(negative["V"]), vmax, 1),
        "I": np.arange(min(negative["I"]), vmax, 1),
    }

    # Only positive
    actual = consolidate_iv(positive, {})
    for mag in actual.keys():
        assert actual[mag] == approx(positive[mag])

    # Only negative
    actual = consolidate_iv({}, negative)
    for mag in actual.keys():
        assert actual[mag] == approx(negative[mag])

    # Both
    actual = consolidate_iv(positive, negative)
    for mag in actual.keys():
        assert actual[mag] == approx(total[mag])
