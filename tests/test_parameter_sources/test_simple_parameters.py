from pytest import raises
from unittest.mock import MagicMock, patch


def test_app_dir(app_local_path):
    from solcore.parameter_system.simple_parameters import app_dir
    from pathlib import Path
    import os

    assert Path.home() != os.path.expanduser("~")
    with raises(OSError):
        app_dir()

    app_path = Path.home() / app_local_path
    os.makedirs(str(app_path))
    assert app_dir() == app_path


def test_locate_source_files_builtin():
    from solcore.parameter_system.simple_parameters import locate_source_files_builtin

    assert len(list(locate_source_files_builtin())) > 0


def test_locate_source_files_in_solcore_app_dir(app_temp_path):
    from solcore.parameter_system.simple_parameters import (
        locate_source_files_builtin,
        locate_source_files_in_solcore_app_dir,
    )
    import shutil

    assert len(list(locate_source_files_in_solcore_app_dir())) == 0
    builtin = list(locate_source_files_builtin())[0]
    shutil.copy(builtin, app_temp_path)

    in_app_dir = list(locate_source_files_in_solcore_app_dir())
    assert len(in_app_dir) == 1
    assert in_app_dir[0].stem == builtin.stem


def test_locate_source_files_in_pwd(tmp_path):
    from solcore.parameter_system.simple_parameters import (
        locate_source_files_builtin,
        locate_source_files_in_pwd,
    )
    import shutil
    import os

    pwd = os.getcwd()
    os.chdir(tmp_path)

    assert len(list(locate_source_files_in_pwd())) == 0
    builtin = list(locate_source_files_builtin())[0]
    shutil.copy(builtin, tmp_path)

    in_app_dir = list(locate_source_files_in_pwd())
    assert len(in_app_dir) == 1
    assert in_app_dir[0].stem == builtin.stem

    os.chdir(pwd)


def test_locate_source_files():
    builtin = MagicMock()
    app_dir = MagicMock()
    pwd = MagicMock()

    package = "solcore.parameter_system.simple_parameters"
    with patch(f"{package}.locate_source_files_builtin", builtin), patch(
        f"{package}.locate_source_files_in_solcore_app_dir", app_dir
    ), patch(f"{package}.locate_source_files_in_pwd", pwd):
        from solcore.parameter_system.simple_parameters import locate_source_files

        files = locate_source_files()

    assert len(files) == 3
    builtin.assert_called_once()
    app_dir.assert_called_once()
    pwd.assert_called_once()


def test_populate_sources():
    from solcore.parameter_system.simple_parameters import populate_sources
    from itertools import chain
    from pathlib import Path

    locations = (
        (Path("file1.json"), Path("file2.json")),
        (Path("file3.json"),),
        (Path("file4.json"),),
    )
    name, path, priority = populate_sources(locations)

    expected = list((n.stem for n in chain.from_iterable(locations)))
    assert name == expected
    expected = {n: Path(n + ".json") for n in expected}
    assert path == expected
    expected = {n.stem: 5 * i for i, loc in enumerate(locations) for n in loc}
    assert priority == expected


class TestSimpleSource:
    def test_load_source(self, simple_data, simple_data_file):
        from solcore.parameter_system import SimpleSource

        SimpleSource._path = {}
        with raises(ValueError, match="'path' for 'source 1' is 'None'"):
            SimpleSource.load_source("source 1")

        SimpleSource._path = {"source 1": simple_data_file}
        with raises(ValueError, match="'priority' for 'source 1' is 'None'"):
            SimpleSource.load_source("source 1")

        SimpleSource._path = {"source 1": simple_data_file}
        SimpleSource._priority = {"source 1": 10}

        ss = SimpleSource.load_source("source 1")

        assert ss.name == "source 1"
        assert ss._priority == 10
        assert ss.reference == simple_data["reference"]
        assert ss._descriptions == simple_data["descriptions"]
        assert ss._data == {
            k: v
            for k, v in simple_data.items()
            if k not in ["reference", "descriptions"]
        }

    def test_materials(self):
        assert False

    def test_parameters(self):
        assert False

    def test_get_parameter(self):
        assert False

    def test__get_parameter_alloy(self):
        assert False

    def test_to_param(self):
        assert False
