from .absorption_QW import calc_alpha, calc_emission
from .adachi_alpha import create_adachi_alpha
from .Custom_CPPB import Custom_CPPB
from .dielectric_constant_models import DielectricConstantModel
from .nk_db import create_nk_txt, download_db, search_db
from .sopra_db import sopra_database
from .tmm_core_vec import (
    ellips,
    find_in_structure,
    inc_tmm,
    position_resolved,
    unpolarized_RT,
)
from .transfer_matrix import (
    OptiStack,
    calculate_absorption_profile,
    calculate_ellipsometry,
    calculate_rat,
)

try:
    from .rigorous_coupled_wave import (
        calculate_absorption_profile_rcwa,
        calculate_rat_rcwa,
        rcwa_rat,
        initialise_S,
        necessary_materials,
        update_epsilon,
        rcwa_position_resolved,
    )
except Exception:
    print(
        "WARNING: The RCWA solver will not be available because an S4 installation has not been found."
    )
