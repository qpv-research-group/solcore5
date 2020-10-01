from glob import glob
from pathlib import Path
from pytest import fixture, mark, skip
import os
from .conftest import patch_plots

cwd = os.getcwd()


@fixture(scope="session")
def examples_directory(tmp_path_factory):
    from shutil import copytree

    exdir = Path(__file__).parent.parent / "examples"
    target = tmp_path_factory.getbasetemp() / "examples"
    copytree(exdir, target)
    return target


def example_scripts():
    exdir = Path(__file__).parent.parent / "examples"
    scripts = sorted([os.path.basename(f) for f in glob(str(exdir / "*.py"))])
    return scripts


def example_notebooks():
    exdir = Path(__file__).parent.parent / "examples" / "notebooks"
    scripts = sorted([os.path.basename(f) for f in glob(str(exdir / "*.ipynb"))])
    return scripts


@patch_plots
@mark.parametrize("example", example_scripts())
def test_example_scripts(example, examples_directory, monkeypatch):
    from solcore.spice import SpiceSolverError
    from solcore.light_source import SmartsSolverError
    from solcore.absorption_calculator import RCWASolverError
    import sys

    os.chdir(examples_directory)
    monkeypatch.setattr("builtins.input", lambda *x, **y: "y")

    script = str(examples_directory / Path(example))
    try:
        with open(script) as fp:
            code = compile(fp.read(), script, "exec")

        # try to emulate __main__ namespace as much as possible
        # What trace python module uses https://docs.python.org/3.7/library/trace.html
        globs = {
            "__file__": script,
            "__name__": "__main__",
            "__package__": None,
            "__cached__": None,
        }
        exec(code, globs, globs)
        os.chdir(cwd)

    except OSError as err:
        os.chdir(cwd)
        skip("Cannot run file %r because: %s" % (sys.argv[0], err))

    except SystemExit:
        os.chdir(cwd)
        pass

    except SpiceSolverError:
        os.chdir(cwd)
        skip("No SPICE solver found.")

    except SmartsSolverError:
        os.chdir(cwd)
        skip("No SMARTS solver found.")

    except RCWASolverError:
        os.chdir(cwd)
        skip("No RCWA solver found.")


@patch_plots
@mark.parametrize("example", example_notebooks())
def test_example_notebooks(example, examples_directory, monkeypatch):
    from solcore.spice import SpiceSolverError
    from solcore.light_source import SmartsSolverError
    from solcore.absorption_calculator import RCWASolverError
    import nbconvert
    import nbformat
    import sys

    notebooks_folder = examples_directory / "notebooks"
    os.chdir(notebooks_folder)
    monkeypatch.setattr("builtins.input", lambda *x, **y: "y")

    script = str(notebooks_folder / Path(example))
    exporter = nbconvert.PythonExporter()
    try:
        nbcode = nbformat.reads(open(script).read(), as_version=4)
        (code, _) = exporter.from_notebook_node(nbcode)
        with open(script) as fp:
            code = compile(code, script, "exec")

        # try to emulate __main__ namespace as much as possible
        # What trace python module uses https://docs.python.org/3.7/library/trace.html
        globs = {
            "__file__": script,
            "__name__": "__main__",
            "__package__": None,
            "__cached__": None,
        }
        exec(code, globs, globs)
        os.chdir(cwd)

    except OSError as err:
        os.chdir(cwd)
        skip("Cannot run file %r because: %s" % (sys.argv[0], err))

    except SystemExit:
        os.chdir(cwd)
        pass

    except SpiceSolverError:
        os.chdir(cwd)
        skip("No SPICE solver found.")

    except SmartsSolverError:
        os.chdir(cwd)
        skip("No SMARTS solver found.")

    except RCWASolverError:
        os.chdir(cwd)
        skip("No RCWA solver found.")
