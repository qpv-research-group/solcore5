import click
import sys
from devpy import util
from devpy.cmds.util import get_site_packages, set_pythonpath, run


@click.command()
@click.option(
    "--build-dir", default="build", help="Build directory; default is `$PWD/build`"
)
@click.argument("codecov_args", nargs=-1)
def codecov(build_dir, codecov_args):
    """🔧 Run codecov in the build directory
    CODECOV_ARGS are passed through directly to codecov, e.g.:
    ./dev.py codecov -- -v
    """

    site_path = get_site_packages(build_dir)
    set_pythonpath(build_dir)

    run(
        ["codecov"] + list(codecov_args),
        cwd=site_path,
    )