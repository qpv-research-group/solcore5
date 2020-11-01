from pathlib import Path

import numpy as np
import pytest


def test_get_built_in_nk_data():
    from solcore.material_data import (
        get_built_in_nk_data,
        MaterialNKDatabaseError,
    )
    import xarray as xr

    GaAs = get_built_in_nk_data("GaAs")
    AlGaAs = get_built_in_nk_data("AlGaAs", composition={"Al": 0.5})

    assert isinstance(GaAs, xr.DataArray)
    assert isinstance(AlGaAs, xr.DataArray)

    with pytest.raises(MaterialNKDatabaseError):
        get_built_in_nk_data("GaAs", composition={"Al": 0.5})


def test_get_sopra_nk_data():
    from solcore.material_data.sopra_db import (
        get_sopra_nk_data,
        MaterialNKDatabaseError,
    )
    import xarray as xr

    GaAs = get_sopra_nk_data("GaAs")
    AlGaAs = get_sopra_nk_data("AlGaAs", composition={"Al": 0.5})

    assert isinstance(GaAs, xr.DataArray)
    assert isinstance(AlGaAs, xr.DataArray)

    with pytest.raises(MaterialNKDatabaseError):
        get_sopra_nk_data("GaAs", composition={"Al": 0.5})


def test_get_refractiveindex_info_nk_data():
    from solcore.material_data import (
        get_refractiveindex_info_nk_data,
        MaterialNKDatabaseError,
        download_db,
    )
    import xarray as xr

    download_db(confirm=True)
    GaAs = get_refractiveindex_info_nk_data("GaAs")
    assert isinstance(GaAs, xr.DataArray)

    with pytest.raises(MaterialNKDatabaseError):
        get_refractiveindex_info_nk_data("AlGaAs")

    with pytest.raises(MaterialNKDatabaseError):
        get_refractiveindex_info_nk_data("GaAs", pageid=-1)


def test_sopra_absorption():
    from solcore.material_data.sopra_db import sopra_database

    # Import material constant data for Gallium Arsenide :: Do this by placing the
    # material name as the sole argument...
    SOPRA_Material = sopra_database("GaAs")

    # Can also load alpha data...
    GaAs_alpha = SOPRA_Material.load_alpha()

    out = GaAs_alpha[1][10]

    data = 163666134.03339368

    assert data == pytest.approx(out)
