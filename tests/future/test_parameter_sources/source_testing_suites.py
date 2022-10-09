from pytest import raises


class NKSourceTestSuite:
    @property
    def source(self):
        return self._source._instance

    def test_load_source(self, *args, **kwargs):
        self._source._instance = None
        assert self.source is None
        source = self._source.load_source()
        assert source is not None
        assert source == self.source

    def test_materials(self):
        assert len(self.source.materials) > 0

    def test_parameters(self):
        assert self.source.parameters("") == ()
        mat = self.source.materials[0]
        assert self.source.parameters(mat) == ("nk",)

    def test_get_parameter(self):
        from solcore.future.parameter import ParameterMissing

        with raises(ParameterMissing):
            self.source.get_parameter("", "")

    def test_get_nk(self):
        mat = self.source.materials[0]
        nk = self.source.get_nk(mat)
        assert nk.dtype == complex
        assert "wavelength" in nk.dims
        assert self.source.name == nk.attrs["reference"]
