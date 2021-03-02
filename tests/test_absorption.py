""" Absorption calculator related tests
"""
from pytest import approx, mark

from solcore import material, si
from solcore.structure import Structure, Layer
from solcore.absorption_calculator import create_adachi_alpha
from solcore.absorption_calculator import (
    calculate_rat,
    calculate_ellipsometry,
    calculate_absorption_profile,
)
from solcore.state import State
from solcore.solar_cell_solver import prepare_solar_cell
from solcore.absorption_calculator.dielectric_constant_models import (
    DielectricConstantModel,
    Drude,
)

from solcore.solar_cell_solver import solar_cell_solver
from solcore.solar_cell import SolarCell
from solcore.material_system import create_new_material
from solcore.absorption_calculator import download_db, search_db
from solcore.absorption_calculator.nk_db import nkdb_load_n
from solcore.config_tools import add_source
from solcore.optics import solve_tmm

import os
from pathlib import Path

import numpy as np


def test_adachi_absorption():
    material_name = "InGaAs"
    solcore_material = material(material_name)(T=300, In=0.1)

    E, nn, kk, adachi_alpha_data = create_adachi_alpha(solcore_material, T=300)

    E /= 1.6e-19

    idx = np.argmin(abs(E - 2))

    out = [E[idx], nn[idx], kk[idx]]

    data = [2.0003267755918639, 3.6288603412514253, 0.30754355994545018]

    assert all([d == approx(o) for d, o in zip(data, out)])


def test_TMM_rat():
    GaAs = material("GaAs")(T=300)

    my_structure = Structure([Layer(si(3000, "nm"), material=GaAs)])

    wavelength = np.linspace(450, 1100, 300)
    idx = np.argmin(abs(wavelength - 800))

    out = calculate_rat(
        my_structure, wavelength, coherent=True, no_back_reflexion=False
    )
    out = (out["R"][idx], out["A"][idx], out["T"][idx])

    data = (0.33328918841743332, 0.65996607786373396, 0.0067447337188326914)

    assert all([d == approx(o) for d, o in zip(data, out)])


def test_TMM_absorption_profile():
    GaAs = material("GaAs")(T=300)

    my_structure = Structure(
        [Layer(si(3000, "nm"), material=GaAs), Layer(si(300, "um"), material=GaAs)]
    )

    out = calculate_absorption_profile(
        my_structure, np.array([800]), z_limit=3000, steps_size=20
    )

    data_path = Path(__file__).parent / "data" / "absorption_profile.txt"
    data = tuple(np.loadtxt(data_path))

    assert all([d == approx(o) for d, o in zip(data, out["absorption"][0])])


def test_TMM_ellipsometry():
    GaAs = material("GaAs")(T=300)

    my_structure = Structure(
        [Layer(si(3000, "nm"), material=GaAs), Layer(si(300, "um"), material=GaAs)]
    )

    wavelength = np.linspace(450, 1100, 300)
    idx = np.argmin(abs(wavelength - 800))

    angles = [60, 65, 70]
    out = calculate_ellipsometry(my_structure, wavelength, angle=angles)

    data = (
        22.2849089096,
        181.488417672,
        16.4604621886,
        182.277656469,
        9.10132195668,
        184.509752582,
    )

    for i in range(len(angles)):
        assert data[2 * i] == approx(out["psi"][idx, i])
        assert data[2 * i + 1] == approx(out["Delta"][idx, i])


def test_TMM_dielectric_model():
    drud = Drude(An=24.317000, Brn=0.125740)

    model = DielectricConstantModel(e_inf=3.4837, oscillators=[drud])

    wavelength = 2 * np.logspace(3, 4, 10)
    n = model.n_and_k(wavelength)

    data = (
        0.37377710 + 2.0726883j,
        0.53555835 + 3.03640188j,
        0.81383835 + 4.12781828j,
        1.25563998 + 5.39395751j,
        1.92396319 + 6.8504966j,
        2.88306008 + 8.48061105j,
        4.17692618 + 10.23820188j,
        5.81068443 + 12.06785126j,
        7.75184851 + 13.93546615j,
        9.95389660 + 15.84816722j,
    )

    assert all([d == approx(o) for d, o in zip(data, n)])


def test_sopra_absorption():
    from solcore.absorption_calculator import sopra_database

    # Import material constant data for Gallium Arsenide :: Do this by placing the
    # material name as the sole argument...
    SOPRA_Material = sopra_database("GaAs")

    # Can also load alpha data...
    GaAs_alpha = SOPRA_Material.load_alpha()

    out = GaAs_alpha[1][10]

    data = 163666134.03339368

    assert data == approx(out)


def test_substrate_presence_A():
    wavelength = np.linspace(300, 800, 3) * 1e-9

    GaAs = material("GaAs")(T=300)

    my_structure = SolarCell([Layer(si(700, "nm"), material=GaAs)], substrate=GaAs)

    solar_cell_solver(
        my_structure,
        "optics",
        user_options={
            "wavelength": wavelength,
            "optics_method": "TMM",
            "no_back_reflexion": False,
        },
    )

    z_pos = np.linspace(0, my_structure.width, 10)

    A_subs = my_structure[0].layer_absorption

    my_structure = SolarCell([Layer(si(700, "nm"), material=GaAs)])

    solar_cell_solver(
        my_structure,
        "optics",
        user_options={
            "wavelength": wavelength,
            "optics_method": "TMM",
            "no_back_reflexion": False,
        },
    )

    A_nosubs = my_structure[0].layer_absorption

    A = np.vstack((A_subs, A_nosubs))

    A_data = np.array(
        [[0.56610281, 0.62692985, 0.41923175], [0.56610281, 0.62711355, 0.37837737]]
    )

    assert all([d == approx(o) for d, o in zip(A, A_data)])


def test_BL_correction():

    wl = np.linspace(800, 950, 4) * 1e-9

    GaAs = material("GaAs")()

    thick_cell = SolarCell([Layer(material=GaAs, width=si("20um"))])

    opts = State()
    opts.position = None
    prepare_solar_cell(thick_cell, opts)
    position = np.arange(0, thick_cell.width, 1e-9)
    opts.position = position
    opts.recalculate_absorption = True
    opts.no_back_reflexion = False

    opts.BL_correction = False
    opts.wavelength = wl
    solve_tmm(thick_cell, opts)

    no_corr = thick_cell.absorbed

    opts.BL_correction = True

    solve_tmm(thick_cell, opts)

    with_corr = thick_cell.absorbed

    assert with_corr == approx(
        np.array([6.71457872e-01, 6.75496354e-01, 2.09738887e-01, 0])
    )
    assert no_corr == approx(
        np.array([6.71457872e-01, 6.75496071e-01, 2.82306407e-01, 0])
    )


def test_substrate_presence_profile():
    wavelength = np.linspace(300, 800, 3) * 1e-9

    GaAs = material("GaAs")(T=300)

    my_structure = SolarCell([Layer(si(700, "nm"), material=GaAs)], substrate=GaAs)

    solar_cell_solver(
        my_structure,
        "optics",
        user_options={
            "wavelength": wavelength,
            "optics_method": "TMM",
            "no_back_reflection": False,
        },
    )

    z_pos = np.linspace(0, my_structure.width, 10)

    profile_subs = my_structure[0].absorbed(z_pos)

    my_structure = SolarCell([Layer(si(700, "nm"), material=GaAs)])

    solar_cell_solver(
        my_structure,
        "optics",
        user_options={
            "wavelength": wavelength,
            "optics_method": "TMM",
            "no_back_reflection": False,
        },
    )

    profile_nosubs = my_structure[0].absorbed(z_pos)
    profile = np.vstack((profile_subs, profile_nosubs))

    data_path = Path(__file__).parent / "data" / "substrate_presence_profile.csv"
    expected = np.loadtxt(data_path, delimiter=",")

    assert profile.shape == expected.shape
    assert profile == approx(expected)


def test_inc_coh_tmm():
    GaInP = material("GaInP")(In=0.5)
    GaAs = material("GaAs")()
    Ge = material("Ge")()

    optical_struct = SolarCell(
        [
            Layer(material=GaInP, width=si("5000nm")),
            Layer(material=GaAs, width=si("200nm")),
            Layer(material=GaAs, width=si("5um")),
            Layer(material=Ge, width=si("50um")),
        ]
    )

    wl = np.linspace(400, 1200, 5) * 1e-9

    options = State()
    options.wavelength = wl
    options.optics_method = "TMM"
    options.no_back_reflection = False
    options.BL_correction = True
    options.recalculate_absorption = True

    c_list = [
        ["c", "c", "c", "c"],
        ["c", "c", "c", "i"],
        ["c", "i", "i", "c"],
        ["i", "i", "i", "i"],
    ]

    results = []
    for i1, cl in enumerate(c_list):
        options.coherency_list = cl
        solar_cell_solver(optical_struct, "optics", options)
        results.append(optical_struct.absorbed)

    A_calc = np.stack(results)
    A_data = np.array(
        [
            [0.5742503, 0.67956899, 0.73481184, 0.725372, 0.76792856],
            [0.5742503, 0.67956899, 0.73481184, 0.725372, 0.76792856],
            [0.5742503, 0.67956899, 0.73474943, 0.70493469, 0.70361194],
            [0.5742503, 0.67956899, 0.70927724, 0.71509221, 0.71592772]
        ]
    )
    assert A_calc == approx(A_data)


def test_define_material():
    from solcore import ParameterSystem, MaterialSystem
    this_dir = os.path.split(__file__)[0]
    create_new_material(
        "SiGeSn",
        os.path.join(this_dir, "data", "SiGeSn_n.txt"),
        os.path.join(this_dir, "data", "SiGeSn_k.txt"),
        os.path.join(this_dir, "data", "SiGeSn_params.txt"),
    )
    assert "SiGeSn" in ParameterSystem().database.sections()
    assert "SiGeSn".lower() in MaterialSystem().sources

    SiGeSn = material("SiGeSn")()
    assert SiGeSn.n(400e-9) == approx(4.175308391752484)
    assert SiGeSn.k(400e-9) == approx(2.3037424963866306)


@mark.skip(reason="Flaky test. Need review.")
def test_database_materials():
    download_db(confirm=True)
    wl, n = nkdb_load_n(2683)  # Should be carbon, from Phillip

    data_path = Path(__file__).parent / "data" / "database_materials.txt"
    n_data = np.loadtxt(data_path)

    assert all([d == approx(o) for d, o in zip(n, n_data)])
