""" Poisson_drift_diffusion related tests
"""
from pytest import approx, mark

import numpy as np
import os
import tempfile

import solcore.poisson_drift_diffusion as PDD
from solcore.solar_cell import SolarCell, default_GaAs
from solcore.structure import Layer, Junction
from solcore import si
from solcore import material
from solcore.light_source import LightSource
from solcore.solar_cell_solver import solar_cell_solver

T = 298

substrate = material("GaAs")(T=T)


def AlGaAs(T):
    # We create the other materials we need for the device
    window = material("AlGaAs")(T=T, Na=5e24, Al=0.8)
    p_AlGaAs = material("AlGaAs")(T=T, Na=1e24, Al=0.4)
    n_AlGaAs = material("AlGaAs")(T=T, Nd=8e22, Al=0.4)
    bsf = material("AlGaAs")(T=T, Nd=2e24, Al=0.6)

    output = Junction(
        [
            Layer(width=si("30nm"), material=window, role="Window"),
            Layer(width=si("150nm"), material=p_AlGaAs, role="Emitter"),
            Layer(width=si("1000nm"), material=n_AlGaAs, role="Base"),
            Layer(width=si("200nm"), material=bsf, role="BSF"),
        ],
        sn=1e6,
        sp=1e6,
        T=T,
        kind="PDD",
    )

    return output


Vin = np.linspace(-2, 2.61, 201)
V = np.linspace(0, 2.6, 300)
wl = np.linspace(350, 1000, 301) * 1e-9
light_source = LightSource(
    source_type="standard",
    version="AM1.5g",
    x=wl,
    output_units="photon_flux_per_m",
    concentration=1,
)


def test_light_iv():
    answer = [
        142.68025180227374,
        2.519346556870366,
        0.9169672186977382,
        329.61395441947565,
        2.347826086956522,
        140.3911287342211,
        0.3294918264376029,
    ]
    with tempfile.TemporaryDirectory(
        prefix="tmp", suffix="_sc3TESTS"
    ) as working_directory:
        filename = os.path.join(working_directory, "solcore_log.txt")
        PDD.log(filename)

        my_solar_cell = SolarCell(
            [AlGaAs(T), default_GaAs(T)], T=T, R_series=0, substrate=substrate
        )
        solar_cell_solver(
            my_solar_cell,
            "iv",
            user_options={
                "T_ambient": T,
                "db_mode": "boltzmann",
                "voltages": V,
                "light_iv": True,
                "wavelength": wl,
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

    for i in range(len(output)):
        assert output[i] == approx(answer[i], rel=1e-2)


def test_quantum_efficiency():
    answer = [
        0.9866334968497021,
        2.1512408472022467e-14,
        0.9779769012349702,
        0.03506561338387434,
    ]
    with tempfile.TemporaryDirectory(
        prefix="tmp", suffix="_sc3TESTS"
    ) as working_directory:
        filename = os.path.join(working_directory, "solcore_log.txt")
        PDD.log(filename)

        my_solar_cell = SolarCell(
            [AlGaAs(T), default_GaAs(T)], T=T, R_series=0, substrate=substrate
        )

        solar_cell_solver(
            my_solar_cell,
            "qe",
            user_options={
                "T_ambient": T,
                "db_mode": "boltzmann",
                "voltages": V,
                "light_iv": True,
                "wavelength": wl,
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

    for i in range(len(output)):
        assert output[i] == approx(answer[i], abs=1e-3)
