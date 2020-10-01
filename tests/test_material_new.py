import pytest


def test_factory():
    from solcore.material_new import Material
    import xarray as xr

    # No parameters nor nk data, just user params
    mat = Material.factory(parametric=False, weight=42)
    assert "weight" in mat.params

    # Get some parameters from the database
    mat = Material.factory(name="GaAs", weight=42)
    assert "band_gap" in mat.params
    assert "weight" in mat.params

    # And now an optical material
    mat = Material.factory(name="GaAs", nk="built-in", parametric=False)
    assert "band_gap" not in mat.params
    assert isinstance(mat.nk, xr.DataArray)

    nk = xr.DataArray([1, 2, 3])
    with pytest.raises(ValueError):
        Material.factory(nk=nk, parametric=False)

    nk = xr.DataArray([1, 2, 3], dims=["wavelength"], coords={"wavelength": [0, 1, 2]})
    mat = Material.factory(nk=nk, parametric=False)
    xr.testing.assert_equal(nk, mat.nk)


def test_from_dict():
    from solcore.material_new import Material
    import xarray as xr

    nk = xr.DataArray([1, 2, 3], dims=["wavelength"], coords={"wavelength": [0, 1, 2]})
    data = dict(name="my mat", T=300, composition={"Fe": 0.1}, nk=nk, weight=42)

    mat = Material.from_dict(data)
    assert mat.name == "my mat"
    assert mat.composition["Fe"] == 0.1
    xr.testing.assert_equal(nk, mat.nk)
    assert "weight" in mat.params

    nk = xr.DataArray([1, 2, 3])
    data = dict(name="my mat", T=300, composition={"Fe": 0.1}, nk=nk, weight=42)
    with pytest.raises(ValueError):
        Material.from_dict(data)


def test_from_dataframe():
    from solcore.material_new import Material
    import xarray as xr
    import pandas as pd

    data = pd.DataFrame.from_dict(
        dict(name="my mat", T=300, composition=[{"Fe": 0.1}], weight=42)
    )

    # No nk data is retrieved
    mat = Material.from_dataframe(data)
    assert mat.name == "my mat"
    assert mat.T == 300
    assert "weight" in mat.params
    assert mat.composition["Fe"] == 0.1
    assert mat.nk is None

    # We ask for nk data, but it is not there
    with pytest.raises(ValueError):
        Material.from_dataframe(data, nk_cols=True)

    # nk data columns are guessed
    nk = xr.DataArray([1, 2, 3], dims=["wavelength"], coords={"wavelength": [0, 1, 2]})
    b = pd.DataFrame.from_dict(
        {"wavelength my mat": nk.wavelength.data, "nk my mat": nk.data}
    )
    data1 = pd.concat((data, b), axis=1)
    mat = Material.from_dataframe(data1, nk_cols=True)
    xr.testing.assert_equal(nk, mat.nk)

    # Custom names are given to the nk columns
    b = pd.DataFrame.from_dict({"wl": nk.wavelength.data, "nk": nk.data})
    data2 = pd.concat((data, b), axis=1)
    mat = Material.from_dataframe(data2, nk_cols=("wl", "nk"))
    xr.testing.assert_equal(nk, mat.nk)


def test_get_attribute():
    from solcore.material_new import Material

    mat = Material.factory(name="GaAs", weight=42)
    assert mat.band_gap
    assert mat.weight == 42


def test_material_str():
    from solcore.material_new import Material

    mat = Material.factory(name="FeO", composition={"Fe": 0.1}, parametric=False)
    assert mat.material_str == "Fe0.1O"


def test_to_dict():
    from solcore.material_new import Material
    import xarray as xr

    nk = xr.DataArray([1, 2, 3], dims=["wavelength"], coords={"wavelength": [0, 1, 2]})
    data = dict(name="my mat", T=300, composition={"Fe": 0.1}, nk=nk, weight=42)
    mat = Material.from_dict(data)

    new_dict = mat.to_dict()
    for k in data.keys():
        if k == "nk":
            continue
        assert data[k] == new_dict[k]

    xr.testing.assert_equal(data["nk"], new_dict["nk"])


def test_to_dataframe():
    from solcore.material_new import Material

    mat = Material.factory(
        name="InGaAs", composition={"In": 0.1}, nk="built-in", weight=42
    )
    df = mat.to_dataframe(include_nk=True)

    for k in ("name", "T", "Na", "Nd", "composition"):
        assert mat.__dict__[k] == df.loc[0, k]

    assert mat.weight == df.loc[0, "weight"]
    assert df.loc[:, "wavelength In0.1GaAs"].values == pytest.approx(
        mat.nk.wavelength.data
    )
    assert df.loc[:, "nk In0.1GaAs"].values == pytest.approx(mat.nk.data)
