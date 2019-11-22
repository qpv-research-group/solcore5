import os
from pathlib import Path

import numpy as np
from pytest import mark, approx

from solcore.absorption_calculator import download_db


def test_define_material(user_config_file):
    from solcore.config_tools import add_source
    from solcore.materials import create_new_material

    home_folder = user_config_file
    custom_nk_path = os.path.join(home_folder, "Solcore/custommats")
    param_path = os.path.join(home_folder, "Solcore/custom_params.cfg")

    add_source("Others", "custom_mats", custom_nk_path)
    add_source("Parameters", "custom", param_path)

    this_dir = os.path.split(__file__)[0]
    create_new_material(
        "SiGeSn",
        os.path.join(this_dir, "data", "SiGeSn_n.txt"),
        os.path.join(this_dir, "data", "SiGeSn_k.txt"),
        os.path.join(this_dir, "data", "SiGeSn_params.txt"),
    )


@mark.xfail(reason="materials database not reloaded after adding new material.")
def test_use_custom_material():
    from solcore import material

    SiGeSn = material("SiGeSn")()
    assert SiGeSn.n(400e-9) == approx(4.175308391752484)
    assert SiGeSn.k(400e-9) == approx(2.3037424963866306)


def test_database_materials(user_config_file):
    from solcore.config_tools import add_source
    from solcore.absorption_calculator.nk_db import nkdb_load_n

    nk_db_path = os.path.join(user_config_file, "NK.db")

    add_source("Others", "nk", nk_db_path)
    download_db(confirm=True, outputfolder=user_config_file)
    wl, n = nkdb_load_n(2683)  # Should be carbon, from Phillip

    data_path = Path(__file__).parent / "data" / "database_materials.txt"
    n_data = np.loadtxt(data_path)

    assert all([d == approx(o) for d, o in zip(n, n_data)])
