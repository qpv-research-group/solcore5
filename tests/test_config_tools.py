from pytest import fixture, raises
from pathlib import Path


@fixture
def user_config(tmpdir):
    return tmpdir / "solcore_config.txt"


@fixture
def default_config():
    return Path(__file__).parent.parent / "solcore" / "solcore_config.txt"


@fixture
def config(default_config, user_config):
    from solcore.config_tools import SolcoreConfig
    return SolcoreConfig(default_config, user_config)


def test_solcore_config_class(config):
    from solcore import SOLCORE_ROOT
    from random import choice
    import os

    assert os.path.isfile(config.user_config)
    assert config.sections == config._default_data.sections()

    s = choice(config.sections)
    o = choice(config.sources(s))
    assert config[s, o] == config._user_data[s][o].replace("SOLCORE_ROOT", SOLCORE_ROOT)
    assert config[s] == config._user_data[s]
    with raises(KeyError):
        assert config["cat"] == config._user_data[s]


def test_add_remove_source(config):
    from random import choice
    from configparser import ConfigParser

    # Added correctly
    s = choice(config.sections)
    config[s, "cat"] = "/kittens/pictures"
    assert "cat" in config.sources(s)
    uconf = ConfigParser()
    uconf.read(config.user_config)
    assert "cat" in uconf.options(s)

    # Removed correctly
    config.remove_source(s, "cat")
    assert "cat" not in config.sources(s)
    uconf = ConfigParser()
    uconf.read(config.user_config)
    assert "cat" not in uconf.options(s)

    # Removing a default Solcore source: it is disabled, instead
    o = choice(config.sources(s))
    config.remove_source(s, o)
    assert config[s, o] == ""
    uconf = ConfigParser()
    uconf.read(config.user_config)
    assert uconf[s][o] == ""


def test_restore_default_source(config):
    from random import choice

    s = choice(config.sections)
    o = choice(config.sources(s))
    value = config[s, o]

    # Wrong source: value unchanged
    config[s, "night"] = "owl"
    assert value == config[s, o]

    # Correct source: value changed
    config[s, o] = "owl"
    assert value != config[s, o]

    # Restored
    config.restore_default_source(s, o)
    assert value == config[s, o]


def test_restore_defaults(config):
    from random import choice

    s = choice(config.sections)
    o = choice(config.sources(s))
    value = config[s, o]

    config[s, "night"] = "owl"
    config[s, o] = "owl"
    assert config[s, "night"] == "owl"
    assert value != config[s, o]

    config.restore_defaults()
    assert config[s, "night"] == "owl"
    assert value == config[s, o]


def test_register_observer(config):
    from random import choice

    def dummy():
        pass

    s = choice(config.sections)
    config.register_observer(s, dummy)

    assert s in config._observers
    assert len(config._observers[s]) == 1
    assert config._observers[s][0] == dummy


def test_notify_observer(config):
    from random import choice

    sources = {}

    def dummy(source, value):
        sources[source] = value

    s = choice(config.sections)
    o = choice(config.sources(s))
    config.register_observer(s, dummy)

    config[s, o] = "spider"
    assert sources[o] == config[s, o]


def test_welcome_message(config):
    config.welcome_message(False)
    assert config.welcome_message() is False


def test_verbose_loading(config):
    config.verbose_loading(False)
    assert config.verbose_loading() is False


def test_spice(config):
    config.spice("neither/here/nor/there")
    assert config.spice() == "neither/here/nor/there"


def test_smarts(config):
    config.smarts("neither/here/nor/there")
    assert config.smarts() == "neither/here/nor/there"


def test_units(config):
    from random import choice

    assert config.units() == config.sources("Units")
    o = choice(config.sources("Units"))
    config.units(o, "neither/here/nor/there")
    assert config.units(o) == "neither/here/nor/there"


def test_parameters(config):
    from random import choice

    assert config.parameters() == config.sources("Parameters")
    o = choice(config.sources("Parameters"))
    config.parameters(o, "neither/here/nor/there")
    assert config.parameters(o) == "neither/here/nor/there"


def test_materials(config):
    from random import choice

    assert config.materials() == config.sources("Materials")
    o = choice(config.sources("Materials"))
    config.materials(o, "neither/here/nor/there")
    assert config.materials(o) == "neither/here/nor/there"
