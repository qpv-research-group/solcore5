import os
import sys
from pathlib import Path

from pytest import fixture


@fixture(autouse=True)
def set_userdir(monkeypatch, tmp_path):
    """Sets the user path to a temp directory."""
    monkeypatch.setattr(Path, "home", lambda: tmp_path)


@fixture
def app_local_path():
    if sys.platform == "win32":
        path = Path("AppData") / "Local" / "solcore"
    elif sys.platform == "darwin":
        path = Path("Library") / "ApplicationSupport" / "solcore"
    else:
        path = Path(".solcore")
    return path


@fixture
def app_temp_path(app_local_path):
    path = Path.home() / app_local_path
    os.makedirs(path)
    return path


@fixture
def parameter_manager():
    from solcore.future.parameter import ParameterManager

    if ParameterManager._instance is not None:
        ParameterManager._instance = None
    return ParameterManager()


@fixture
def simple_data():
    return {
        "reference": "A paper",
        "descriptions": {
            "param1": "The first parameter",
            "param2": "The second parameter",
            "param3": "The third parameter",
        },
        "Dark matter": {"param1": 42, "param2": "23 eV", "param3": "2*T K"},
        "Light matter": {"param1": 84, "param2": "42 eV"},
    }


@fixture
def simple_data_file(tmp_path, simple_data):
    import json

    path = tmp_path / "data.json"
    with path.open("a") as f:
        json.dump(simple_data, f)
    return path


@fixture
def simple_param_source(simple_data_file):
    from solcore.future.parameter_sources import SimpleSource

    NewSource = type(
        f"Source1{SimpleSource.__name__}",
        (SimpleSource,),
        {"name": "Source1", "_path": simple_data_file, "_priority": 1},
    )

    return NewSource.load_source()


@fixture(scope="session")
def nk_database(tmp_path_factory):
    return tmp_path_factory.mktemp("nk_folder") / "nk.db"
