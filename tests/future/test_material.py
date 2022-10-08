class TestMaterial:
    def test_validate_args(self):
        from solcore.future.material import Material
        from pint import Quantity as Q_

        with_units = Material._validate_args(T=300, Na=1e23, Nd=1e22, band_gap=1.42)

        assert with_units["T"] == Q_(300, "K")
        assert with_units["Na"] == Q_(1e23, "1/m**3")
        assert with_units["Nd"] == Q_(1e22, "1/m**3")
        assert with_units["band_gap"] == Q_(1.42, "dimensionless")

    def test_factory(self):
        from solcore.future.material import Material
        from pint import Quantity as Q_
        import xarray as xr

        nk = xr.DataArray(
            [1, 2, 3], dims=["wavelength"], coords={"wavelength": [0, 1, 2]}
        )
        mat = Material.factory(
            name="GaAs", include=("band_gap",), T=Q_(300, "K"), nk=nk
        )
        assert "band_gap" in mat._params
        assert "T" in mat._params
        xr.testing.assert_equal(nk, mat._nk)

    def test_params(self):
        from solcore.future.material import Material
        from pint import Quantity as Q_

        mat = Material.factory(name="GaAs", include=("band_gap",), T=Q_(300, "K"))
        assert all(p in mat.params for p in mat._params.keys())

    def test_get_attribute(self):
        from solcore.future.material import Material
        from pint import Quantity as Q_

        mat = Material.factory(name="GaAs", T=Q_(300, "K"))
        assert mat.T == Q_(300, "K")
        assert "band_gap" not in mat._params
        mat.band_gap
        assert "band_gap" in mat._params

    def test_nk(self):
        from solcore.future.material import Material

        mat = Material.factory(name="GaAs")
        assert mat._nk.shape == ()
        mat.nk
        assert mat._nk.shape != ()

    def test_material_str(self):
        from solcore.future.material import Material

        mat = Material.factory(name="FeO", comp={"Fe": 0.1})
        assert mat.material_str == "Fe0.1O"


def test_material():
    from solcore.future.material import material, Material
    from functools import partial

    InGaAs = material("InGaAs")
    assert isinstance(InGaAs, partial)
    n_InGaAs = InGaAs(T=300, In=0.2, Na=1e23)
    other_InGaAs = Material.factory("InGaAs", comp={"In": 0.2}, Na=1e23)

    assert n_InGaAs.name == other_InGaAs.name
    assert n_InGaAs.comp == other_InGaAs.comp
    assert n_InGaAs.Na == other_InGaAs.Na
