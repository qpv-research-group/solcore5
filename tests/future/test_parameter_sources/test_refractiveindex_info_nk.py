from pytest import raises, mark
from contextlib import contextmanager
import os

from source_testing_suites import NKSourceTestSuite
from solcore.future.parameter_sources import RefractiveindexInfoNKSource


@contextmanager
def working_directory(path):
    """A context manager for temporary changing the working directory.

    Taken from:
        https://gist.github.com/nottrobin/3d675653244f8814838a

    """
    prev_cwd = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(prev_cwd)


@mark.xfail
class TestRefractiveindexInfoNKSource(NKSourceTestSuite):
    _source = RefractiveindexInfoNKSource

    def test_download_db(self, nk_database):
        assert not nk_database.is_file()
        source = self._source.download_db(nk_database)
        assert self.source == source
        assert nk_database.is_file()
        self._source._instance = None

    def test_select_db(self, nk_database):
        with raises(FileNotFoundError):
            self._source.select_db("")

        source = self._source.select_db(nk_database)
        assert self.source == source
        self._source._instance = None

    def test_load_source(self, nk_database):
        import os

        super().test_load_source()
        assert self.source._db is None
        self._source._instance = None

        with working_directory(str(nk_database.parent)):
            super().test_load_source()
            assert self.source._db is not None
            assert self.source._db.db_path == str(nk_database)
            self._source._instance = None

        os.environ["SOLCORE_REFRACTIVE_INDEX_DB"] = str(nk_database)
        super().test_load_source()
        assert self.source._db is not None
        assert self.source._db.db_path == str(nk_database)

    def test_materials(self):
        super().test_materials()

    def test_parameters(self):
        super().test_parameters()

    def test_get_parameter(self):
        super().test_get_parameter()

    def test_get_nk(self):
        from solcore.future.parameter import InputArgumentMissing

        mat = self.source.materials[0]

        with raises(InputArgumentMissing):
            self.source.get_nk(mat)

        pageid = self.source._db.get_available_for_material(mat).index[0]
        nk = self.source.get_nk(mat, pageid=pageid)
        assert nk.dtype == complex
        assert "wavelength" in nk.dims
        assert self.source.name == nk.attrs["reference"]
