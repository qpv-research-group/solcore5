from .adachi_alpha import create_adachi_alpha
from .absorption_QW import calc_alpha, calc_emission
from .transfer_matrix import calculate_absorption_profile, calculate_rat, calculate_ellipsometry, OptiStack
from .dielectric_constant_models import DielectricConstantModel
from .sopra_db import sopra_database
from .nk_db import download_db, search_db, create_nk_txt
from .cppm import Custom_CPPB
from .tmm_core_vec import inc_tmm, position_resolved, inc_tmm, unpolarized_RT, ellips, find_in_structure


class RCWASolverError(Exception):
    pass


try:
    from .rigorous_coupled_wave import calculate_absorption_profile_rcwa, calculate_rat_rcwa, \
        rcwa_rat, initialise_S, necessary_materials, update_epsilon, rcwa_position_resolved
except RCWASolverError:
    print('WARNING: The RCWA solver will not be available because an S4 installation has not been found.')
