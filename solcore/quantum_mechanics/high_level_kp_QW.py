from solcore.quantum_mechanics.kp_QW import solve_bandstructure_QW
from solcore.quantum_mechanics.heterostructure_alignment import VBO_align
from solcore.quantum_mechanics.structure_utilities import structure_to_potentials
from solcore.quantum_mechanics.potential_utilities import potentials_to_wavefunctions_energies
from solcore.quantum_mechanics import graphics
from solcore.absorption_calculator.absorption_QW import calc_alpha
from solcore.interpolate import interp1d
from solcore.constants import q, vacuum_permittivity
import numpy as np


def schrodinger(structure, plot_bands=False, kpoints=40, krange=1e9, num_eigenvalues=10, symmetric=True,
                quasiconfined=0.0, return_qw_boolean_for_layer=False, Efield=0, blur=False, blurmode="even",
                mode='kp8x8_bulk', step_size=None, minimum_step_size=0, smallest_feature_steps=20, filter_strength=0,
                periodic=False, offset=0, graphtype=[], calculate_absorption=False, alpha_params=None, **kwargs):
    """ Solves the Schrodinger equation of a 1 dimensional structure. Depending on the inputs, the method for solving
    the problem is more or less sophisticated. In all cases, the output includes the band structure and effective
    masses around k=0 as a function of the position, the energy levels of the QW (electrons and holes) and the
    wavefunctions, although for the kp4x4 and kp6x6 modes, this is provided as a function of the k value, and therefore
    there is much more information.

    :param structure: The strucutre to solve
    :param plot_bands: (False) If the bands should be plotted
    :param kpoints: (30) The number of points in the k space
    :param krange: (1e-9) The range in the k space
    :param num_eigenvalues: (10) Maximum number of eigenvalues to calculate
    :param symmetric: (True) If the structure is symmetric, in which case the calculation can be speed up
    :param quasiconfined: (0.0 eV) Energy above the band edges that an energy level can have before rejecting it
    :param return_qw_boolean_for_layer: (False) Return an boolean array indicating which positions are inside the QW
    :param Efield: (0) Electric field.
    :param blur: (False) If the potentials and effective masses have to be blurred
    :param blurmode: ('even') Other values are 'right' and 'left'
    :param mode: ('kp4x4') The mode of calculating the bands and effective masses. See 'structure_utilities.structure_to_potentials'
    :param step_size: (None) The discretization step of the structure. If none, it is estimated based on the smallest feature.
    :param minimum_step_size: (0) The minimum step size.
    :param smallest_feature_steps: (20) The number of steps in the smallest feature
    :param filter_strength: (0) If > 0, defines the fraction of the wavefunction that has to be inside the QW in order to consider it 'confined'
    :param periodic: (False) If the structure is periodic. Affects the boundary conditions.
    :param offset: (0) Energy offset used in the calculation of the energy levels in the case of the 'bulk' solvers
    :param graphtype: [] If 'potential', the band profile and wavefunctions are ploted
    :return: A dictionary containing the band structure and wavefunctions as a function of the position and k
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
                                                     quasiconfined=quasiconfined, **potentials)

        result_band_edge = {
            "x": bands['x'],
            "potentials": {key: potentials[key] for key in potentials.keys() if key[0] in "Vx"},
            "effective_masses": {key: potentials[key] for key in potentials.keys() if key[0] in "mx"},
            "wavefunctions": {key: bands[key] for key in bands.keys() if 'psi' in key},
            "E": {key: np.array(bands[key]) for key in bands.keys() if key[0] in "E"},
        }

    if "potentials" is graphtype:
        schrodinger_plt = graphics.split_schrodinger_graph_potentials(result_band_edge, **kwargs)
        # schrodinger_plt.draw()

    if "potentialsLDOS" is graphtype:
        Ee, LDOSe, Eh, LDOSh = graphics.split_schrodinger_graph_LDOS(result_band_edge, **kwargs)
        result_band_edge["LDOS"] = {"x": bands['x'], 'Ee': Ee, 'LDOSe': LDOSe, 'Eh': Eh, 'LDOSh': LDOSh}

    if calculate_absorption:
        result_band_edge["alpha"] = calc_alpha(result_band_edge, **alpha_params)
        result_band_edge["alphaE"] = interp1d(x=result_band_edge["alpha"][0], y=result_band_edge["alpha"][1])

    return result_band_edge, bands


if __name__ == "__main__":
    from solcore import si, material
    from solcore.structure import Layer, Structure
    import matplotlib.pyplot as plt
    import numpy as np

    bulk = material("GaAs")(T=293)
    barrier = material("GaAsP")(T=293, P=0.1)

    bulk.strained = False
    barrier.strained = True

    top_layer = Layer(width=si("30nm"), material=bulk)
    inter = Layer(width=si("3nm"), material=bulk)
    barrier_layer = Layer(width=si("15nm"), material=barrier)
    bottom_layer = top_layer

    E = np.linspace(1.15, 1.5, 300) * q
    alfas = np.zeros((len(E), 6))

    alfas[:, 0] = E / q

    alpha_params = {
        "well_width": si("7.2nm"),
        "theta": 0,
        "eps": 12.9 * vacuum_permittivity,
        "espace": E,
        "hwhm": si("6meV"),
        "dimensionality": 0.16,
        "line_shape": "Gauss"
    }

    comp = [0.05, 0.10, 0.15, 0.20]
    colors = plt.cm.jet(np.linspace(0, 1, len(comp)))

    plt.figure(figsize=(6, 4.5))
    for j, i in enumerate(comp):
        QW = material("InGaAs")(T=293, In=i)
        QW.strained = True
        well_layer = Layer(width=si("7.2nm"), material=QW)

        # test_structure = Structure([top_layer, barrier_layer, inter] + 1 * [well_layer, inter, barrier_layer, inter] +
        #                            [bottom_layer])

        test_structure = Structure([barrier_layer, inter] + 1 * [well_layer, inter] +
                                   [barrier_layer])

        # test_structure = Structure([top_layer, barrier_layer] + 10 * [well_layer, barrier_layer] +
        #                            [bottom_layer])

        test_structure.substrate = bulk

        output = schrodinger(test_structure, quasiconfined=0,  mode='kp4x4', plot_bands=False,
                             num_eigenvalues=20, alpha_params=alpha_params, calculate_absorption=True)

        alfa = output[0]['alphaE'](E)
        plt.plot(1240 / (E / q), alfa / 100, label='{}%'.format(int(i * 100)))
        alfas[:, j + 1] = alfa / 100

    plt.xlim(826, 1100)
    plt.ylim(0, 23000)
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('$\\alpha$ cm$^{-1}$')
    plt.legend(loc='upper right', frameon=False)
    plt.tight_layout()

    import os

    root = os.path.expanduser('~')
    plt.savefig(root + '/Desktop/abs.pdf')
    plt.show()
