from ..registries import register_optics
from .beer_lambert import solve_beer_lambert  # noqa
from .external_optics import solve_external_optics  # noqa
from .rcwa import rcwa_options, solve_rcwa  # noqa
from .tmm import solve_tmm  # noqa


@register_optics(name=None)
def no_optics(*args, **kwargs) -> None:
    """Dummy function that does not calculate any optics."""
    return
