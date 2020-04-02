from glob import glob
from pathlib import Path
from pytest import fixture, mark, skip
import os
from .conftest import patch_plots

# nbdir = exdir / "notebooks"
# ndfiles = glob(str(nbdir) + "*.ipynb")

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
    scripts = [os.path.basename(f) for f in glob(str(exdir / "*.py"))]
    scripts.remove("MJ_solar_cell_using_DA.py")
    scripts.remove("MJ_solar_cell_tutorial.py")
    return scripts


@patch_plots
@mark.parametrize("example", example_scripts())
def test_example_scripts(example, examples_directory):
    from solcore.spice import SpiceSolverError
    from solcore.light_source import SmartsSolverError
    import sys
    os.chdir(examples_directory)

    script = examples_directory / Path(example)
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
        sys.exit("Cannot run file %r because: %s" % (sys.argv[0], err))

    except SystemExit:
        os.chdir(cwd)
        pass

    except SpiceSolverError:
        os.chdir(cwd)
        skip("No SPICE solver found.")

    except SmartsSolverError:
        os.chdir(cwd)
        skip("No SMARTS solver found.")
