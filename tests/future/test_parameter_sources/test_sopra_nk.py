from pytest import raises

from source_testing_suites import NKSourceTestSuite
from solcore.future.parameter_sources import SopraNKSource


class TestSopraNKSource(NKSourceTestSuite):
    _source = SopraNKSource

    def test_load_alloy(self):
        with raises(NotImplementedError):
            self.source._load_alloy("", 0)
