from pytest import raises, approx
from unittest.mock import patch

from source_testing_suites import NKSourceTestSuite
from solcore.future.parameter_sources import BuiltinNKSource


class TestBuiltinNKSource(NKSourceTestSuite):
    _source = BuiltinNKSource

    def test_get_nk(self):
        from solcore.future.parameter import InputArgumentMissing

        nk = self.source.get_nk("GaAs")
        assert nk.dtype == complex
        assert "wavelength" in nk.dims
        assert self.source.name == nk.attrs["reference"]

        with raises(InputArgumentMissing):
            self.source.get_nk("InGaAs")

        nk = self.source.get_nk("InGaAs", comp={"In": 0.1})
        assert nk.dtype == complex
        assert "wavelength" in nk.dims
        assert self.source.name == nk.attrs["reference"]

    def test_load(self, tmp_path):
        import numpy as np

        data = np.random.random((10, 2))
        path = tmp_path / "n.txt"
        np.savetxt(path, data)
        (wl, n) = self.source._load(path.parent, "n")
        assert wl == approx(data[:, 0])
        assert n == approx(data[:, 1])

    def test_load_alloy(self, tmp_path):
        import numpy as np

        data = np.random.random((10, 2))

        def mock_load(*args, **kwargs):
            return {"nk": (data[:, 0], None)}, None

        def mock_interpolate(*args, **kwargs):
            return None, data[:, 1], None

        package = "solcore.future.parameter_sources.tools.critical_point_interpolate"
        with patch(f"{package}.load_data_from_directory", mock_load), patch(
            f"{package}.critical_point_interpolate", mock_interpolate
        ):

            (wl, n) = self.source._load_alloy(tmp_path, "n", 0.1)
            assert wl == approx(data[:, 0])
            assert n == approx(data[:, 1])
