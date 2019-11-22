from .beer_lambert import solve_beer_lambert
from .external_optics import solve_external_optics
from .tmm import solve_tmm

try:
    from .rcwa import solve_rcwa, rcwa_options
except ImportError:
    rcwa_options = {}
    solve_rcwa = None
