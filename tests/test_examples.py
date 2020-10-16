from glob import glob
from pathlib import Path
from pytest import fixture, mark, skip
import os
import subprocess
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
def test_example_scripts(example, examples_directory):

    env = os.environ.copy()
    env["PYTHONWARNINGS"] = "ignore"
    with open(examples_directory / example) as fp:
        script = "import matplotlib\nmatplotlib.use('Agg')\n" + fp.read()
        process = subprocess.run(
            ["python", "-c", script],
            capture_output=True,
            cwd=str(examples_directory),
            env=env,
        )

    if b"SpiceSolverError" in process.stderr:
        skip("No SPICE solver found.")
    elif b"SmartsSolverError" in process.stderr:
        skip("No SMARTS solver found.")
    elif b"RCWASolverError" in process.stderr:
        skip("No RCWA solver found.")
    elif process.stderr != b"":
        raise Exception(process.stderr.decode())


@patch_plots
@mark.parametrize("example", example_notebooks())
def test_example_notebooks(example, examples_directory):
    import nbconvert
    import nbformat

    exporter = nbconvert.PythonExporter()
    notebooks_folder = examples_directory / "notebooks"
    env = os.environ.copy()
    env["PYTHONWARNINGS"] = "ignore"
    with open(notebooks_folder / example) as fp:
        nbcode = nbformat.reads(fp.read(), as_version=4)
        (code, _) = exporter.from_notebook_node(nbcode)
        script = "import matplotlib\nmatplotlib.use('Agg')\n" + code
        process = subprocess.run(
            ["python", "-c", script],
            capture_output=True,
            cwd=str(notebooks_folder),
            env=env,
        )

    if b"SpiceSolverError" in process.stderr:
        skip("No SPICE solver found.")
    elif b"SmartsSolverError" in process.stderr:
        skip("No SMARTS solver found.")
    elif b"RCWASolverError" in process.stderr:
        skip("No RCWA solver found.")
    elif process.stderr != b"":
        raise Exception(process.stderr.decode())
