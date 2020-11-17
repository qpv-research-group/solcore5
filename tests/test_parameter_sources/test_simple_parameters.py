from pytest import raises
from unittest.mock import MagicMock, patch


def test_app_dir(app_local_path):
    from solcore.parameter_sources.simple_parameters import app_dir
    from pathlib import Path
    import os

    assert Path.home() != os.path.expanduser("~")
    with raises(OSError):
        app_dir()

    app_path = Path.home() / app_local_path
    os.makedirs(str(app_path))
    assert app_dir() == app_path


def test_locate_source_files_builtin():
    from solcore.parameter_sources.simple_parameters import locate_source_files_builtin

    assert len(list(locate_source_files_builtin())) > 0


def test_locate_source_files_in_solcore_app_dir(app_temp_path):
    package = "solcore.parameter_sources.simple_parameters"
    with patch(f"{package}.app_dir", lambda: app_temp_path):
        from solcore.parameter_sources.simple_parameters import (
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
    from solcore.parameter_sources.simple_parameters import (
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

    package = "solcore.parameter_sources.simple_parameters"
    with patch(f"{package}.locate_source_files_builtin", builtin), patch(
        f"{package}.locate_source_files_in_solcore_app_dir", app_dir
    ), patch(f"{package}.locate_source_files_in_pwd", pwd):
        from solcore.parameter_sources.simple_parameters import locate_source_files

        files = locate_source_files()

    assert len(files) == 3
    builtin.assert_called_once()
    app_dir.assert_called_once()
    pwd.assert_called_once()


def test_populate_sources():
    from solcore.parameter_sources.simple_parameters import populate_sources
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
        from solcore.parameter_sources import SimpleSource

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

    def test_materials(self, simple_param_source):
        ss = simple_param_source
        assert ss.materials == tuple((m for m in ss._data.keys()))

    def test_parameters(self, simple_param_source):
        from solcore.parameter import ParameterError

        ss = simple_param_source

        with raises(ParameterError):
            ss.parameters("ice")

        actual = ss.parameters("Dark matter")
        assert actual == tuple(ss._data["Dark matter"].keys())

    def test_get_parameter(self, simple_param_source):
        from solcore.parameter import ParameterError

        ss = simple_param_source
        ss._get_parameter_alloy = MagicMock(return_value=84)
        ss.to_param = MagicMock(return_value=42)

        msg = f"Material 'ice' not in '{ss.name}' source."
        with raises(ParameterError, match=msg):
            ss.get_parameter(material="ice", parameter="param1")

        msg = f"Parameter 'mass' not available for 'Dark matter' in source '{ss.name}'."
        with raises(ParameterError, match=msg):
            ss.get_parameter(material="Dark matter", parameter="mass")

        raw = ss._data["Dark matter"]["param1"]
        out = ss.get_parameter(material="Dark matter", parameter="param1")
        assert out == 42
        ss.to_param.assert_called_with(raw, "param1")

        ss.to_param.reset_mock()
        ss._data["Dark matter"]["x"] = "Lp"
        out = ss.get_parameter(material="Dark matter", parameter="param1")
        assert out == 84
        ss._get_parameter_alloy.assert_called_with("Dark matter", "param1")
        ss.to_param.assert_not_called()

    def test__get_parameter_alloy(self, simple_param_source):
        from solcore.parameter import ParameterError

        ss = simple_param_source
        ss.get_parameter = MagicMock(side_effect=[2, 3])
        ss.to_param = MagicMock(side_effect=[0, 42])

        ss._data["Dark matter"]["x"] = "Lp"
        msg = "Ternary alloys must have 'parent0' and 'parent1'"
        with raises(ParameterError, match=msg):
            ss._get_parameter_alloy(material="Dark matter", parameter="param1")

        ss.get_parameter = MagicMock(side_effect=[2, 3])
        ss.to_param = MagicMock(side_effect=[0, 42])
        ss._data["Dark matter"]["x"] = "Lp"
        ss._data["Dark matter"]["parent0"] = "mamma"
        ss._data["Dark matter"]["parent1"] = "papa"
        msg = f"Composition for element {ss._data['Dark matter']['x']} is required"
        with raises(KeyError, match=msg):
            ss._get_parameter_alloy(material="Dark matter", parameter="param1")

        ss.get_parameter = MagicMock(side_effect=[2, 3])
        ss.to_param = MagicMock(side_effect=[0, 42])
        comp = {"Lp": 0.3}
        out = ss._get_parameter_alloy(
            material="Dark matter", parameter="param1", comp=comp
        )
        assert out == 42

    def test_to_param(self, simple_param_source):
        from solcore.parameter import Parameter

        ss = simple_param_source
        out = ss.to_param(42, "the answer")
        assert out.m == 42

        out = ss.to_param("42 eV", "the answer")
        assert out.m == 42
        assert out.u == "electron_volt"

        out = ss.to_param(Parameter("42 eV"), "the answer")
        assert out.m == 42
        assert out.u == "electron_volt"

        with raises(NameError):
            ss.to_param("42*T eV", "the answer")

        out = ss.to_param("42*T eV", "the answer", T=2)
        assert out.m == 84
        assert out.u == "electron_volt"
