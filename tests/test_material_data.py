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
