from .kp_QW import solve_bandstructure_QW
from .heterostructure_alignment import VBO_align
from .structure_utilities import structure_to_potentials
from .potential_utilities import potentials_to_wavefunctions_energies
from . import graphics
from solcore.absorption_calculator.absorption_QW import calc_alpha
from solcore.interpolate import interp1d


def schrodinger(structure, plot_bands=False, kpoints=40, krange=1e9, num_eigenvalues=10, symmetric=True,
                quasiconfined=0.0, return_qw_boolean_for_layer=False, Efield=0, blur=False, blurmode="even",
                mode='kp4x4', step_size=None, minimum_step_size=0, smallest_feature_steps=20, filter_strength=0,
                periodic=False, offset=0, graphtype=[], calculate_absorption=False, alpha_params=None):
    """ Solves the Schrodinger equation of a 1 dimensional structure. Depending on the inputs, the method for solving
    the problem is more or less sophisticated. In all cases, the output includes the band structure and effective
    masses around k.txt=0 as a function of the position, the energy levels of the QW (electrons and holes) and the
    wavefunctions, although for the kp4x4 and kp6x6 modes, this is provided as a function of the k.txt value, and therefore
    there is much more information.

    :param structure: The strucutre to solve
    :param plot_bands: (False) If the bands should be plotted
    :param kpoints: (30) The number of points in the k.txt space
    :param krange: (1e-9) The range in the k.txt space
    :param num_eigenvalues: (10) Maximum number of eigenvalues to calculate
    :param symmetric: (True) If the structure is symmetric, in which case the calculation can be speed up
    :param quasiconfined: (0.0) Energy above the band edges that an energy level can have before rejecting it
    :param return_qw_boolean_for_layer: (False) Return an boolean array indicating which positions are inside the QW
    :param Efield: (0) Electric field.
    :param blur: (False) If the potentials and effective masses have to be blurred
    :param blurmode: ('even') Other values are 'right' and 'left'
    :param mode: ('kp4x4') The mode of calculating the bands and effective masses. See 'structure_utilities.structure_to_potentials'
    :param step_size: (None) The discretization step of the structure. If none, it is estimated based on the smallest feature.
    :param minimum_step_size: (0) The minimum step size.
    :param smallest_feature_steps: (20) The number of steps in the smallest feature
    :param filter_strength: (0) If > 0, defines the fraction of the wavefunction that has to be inside the QW in order to consider it 'confined'
    :param periodic: (False) If the strucuture is periodic. Affects the boundary conditions.
    :param offset: (0) Energy offset used in the calculation of the energy levels in the case of the 'bulk' solvers
    :param graphtype: [] If 'potential', the band profile and wavefunctions are ploted
    :return: A dictionary containing the band structure and wavefunctions as a function of the position and k.txt
    """

    if Efield != 0:
        symmetric = False

    aligned_structure = VBO_align(structure)
    potentials = structure_to_potentials(aligned_structure, return_qw_boolean_for_layer=return_qw_boolean_for_layer,
                                         Efield=Efield, blur=blur, blurmode=blurmode, mode=mode, step_size=step_size,
                                         minimum_step_size=minimum_step_size,
                                         smallest_feature_steps=smallest_feature_steps)

    # Now that we have the potential and effective masses ready, we solve the problem.
    if mode in ['kp4x4', 'kp6x6']:
        bands = solve_bandstructure_QW(potentials, num=num_eigenvalues, kpoints=kpoints, krange=krange,
                                       symmetric=symmetric, quasiconfined=quasiconfined, plot_bands=plot_bands)

        result_band_edge = {
            "x": bands['x'],
            "potentials": {key: potentials[key] for key in potentials.keys() if key[0] in "Vx"},
            "effective_masses": {key: potentials[key] for key in potentials.keys() if key[0] in "mx"},
            "wavefunctions": {key: bands[key][:, 0] for key in bands.keys() if 'psi' in key},
            "E": {key: bands[key][:, 0] for key in bands.keys() if key[0] in 'E'},
        }

    else:
        bands = potentials_to_wavefunctions_energies(structure=structure, num_eigenvalues=num_eigenvalues,
                                                     filter_strength=filter_strength, offset=offset, periodic=periodic,
                                                     **potentials)

        result_band_edge = {
            "x": bands['x'],
            "potentials": {key: potentials[key] for key in potentials.keys() if key[0] in "Vx"},
            "effective_masses": {key: potentials[key] for key in potentials.keys() if key[0] in "mx"},
            "wavefunctions": {key: bands[key] for key in bands.keys() if 'psi' in key},
            "E": {key: bands[key] for key in bands.keys() if key[0] in "E"},
        }

    if "potentials" in graphtype:
        schrodinger_plt = graphics.split_schrodinger_graph(result_band_edge)
        schrodinger_plt.draw()

    if calculate_absorption:
        result_band_edge["alpha"] = calc_alpha(result_band_edge, **alpha_params)
        result_band_edge["alphaE"] = interp1d(x=result_band_edge["alpha"][0], y=result_band_edge["alpha"][1])

    return result_band_edge, bands
