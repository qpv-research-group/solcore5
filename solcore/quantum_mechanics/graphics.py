import copy
from solcore.graphing import Graph, GraphData, graph_defaults
from solcore.graphing.graph_support import open_with_os
from solcore.constants import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import FigureCanvasPdf
import tempfile

defaults = {
    "edit": lambda x, y: (x * 1e9, y / q),
    "xlabel": "Depth (nm)",
    "ylabel": "Energy (eV)",
}


def L(x, centre, hwhm):
    return 1 / pi * (0.5 * hwhm) / ((x - centre) ** 2 + (0.5 * hwhm) ** 2)  # Lorenzian (area normalised to 1)


def structure_graph(structure, **kwargs):
    global defaults
    options = copy.copy(defaults)
    options.update(kwargs)

    current_x = 0
    x = []
    conduction = []
    valence = []
    for layer in structure:
        x.append(current_x)
        conduction.append(layer.Ec)
        conduction.append(layer.Ec)
        valence.append(layer.Ev)
        valence.append(layer.Ev)
        current_x += layer.width
        x.append(current_x)

    # x = np.array(x)
    # conduction = np.array(conduction)
    # valence = np.array(valence)
    g = Graph(
        GraphData(x, conduction, label="Conduction Band", linewidth=2, color="black"),
        GraphData(x, valence, label="Valence Band", linewidth=2, color="black"),
        **options
    )
    return g


def normalise_psi(p):
    return p / (max(p) - min(p)) * q


def prepare_wavefunction_data(x, E, Psi, trim_levels_beyond, linewidth, color, scale, suppress_invert=False,
                              square=False, alpha=1, label="untitled"):
    data = []
    for e, psi in zip(E, Psi):
        norm_psi = normalise_psi(psi if not square else psi ** 2)
        trim_to = norm_psi ** 2 / q ** 2 > trim_levels_beyond ** 2
        trimmed_x = x[trim_to]
        norm_psi = norm_psi[trim_to]
        invert = -1 if norm_psi[0] < 0 and not suppress_invert and not square else 1

        linekwargs = {"linewidth": linewidth, "alpha": alpha, "label": label}
        data.append(GraphData(trimmed_x, e + norm_psi * scale * invert, color=color, **linekwargs))
    return data


def wavefunctions_graph(
        x, Ee=None,
        psi_e=None,
        Ehh=None,
        psi_hh=None,
        Elh=None,
        psi_lh=None,
        trim_levels_beyond=1e-3,
        linewidth=1,
        scale=0.06,
        colour_by_band=True,
        **kwargs):
    global defaults
    options = copy.copy(defaults)
    options.update(kwargs)

    data = []

    # normalise_psi = lambda p:p/numpy.sqrt(numpy.trapz(x=x*1e9,y=p**2))*q

    if Ee is not None: data.extend(
        prepare_wavefunction_data(x, Ee, psi_e, trim_levels_beyond, linewidth, "blue", scale, label="e"))
    if Ehh is not None: data.extend(
        prepare_wavefunction_data(x, Ehh, psi_hh, trim_levels_beyond, linewidth, "green", scale, label="hh"))
    if Elh is not None: data.extend(
        prepare_wavefunction_data(x, Elh, psi_lh, trim_levels_beyond, linewidth, "red", scale, label="lh"))

    g = Graph(data, **options)
    return g


def potentials_graph(x, Ve=None, Vhh=None, Vlh=None, color="grey", **kwargs):
    global defaults
    options = copy.copy(defaults)
    options.update(kwargs)
    data = []

    normalise_psi = lambda p: p * q * 5e-6

    if Ve is not None: data.append(GraphData(x, Ve, color=color, linewidth=2, label="Ve"))
    if Vlh is not None: data.append(GraphData(x, Vlh, color=color, linewidth=2, dashes=[1, 1], label="Vlh"))
    if Vhh is not None: data.append(GraphData(x, Vhh, color=color, linewidth=2, label="Vhh"))

    g = Graph(data, **options)
    return g


def schrodinger_graph(schrodinger_result, **kwargs):
    potentials_kwargs = copy.copy(schrodinger_result["potentials"])
    potentials_kwargs.update(kwargs)
    potential_graph = potentials_graph(**potentials_kwargs)

    wavefunctions_kwargs = copy.copy(schrodinger_result["wavefunctions"])
    wavefunctions_kwargs.update(kwargs)
    wavefunctions_kwargs.update(schrodinger_result["E"])

    w = wavefunctions_graph(figure=potential_graph.figure, trim_levels_beyond=1e-3, **wavefunctions_kwargs)
    return w


def split_schrodinger_graph(schrodinger_result,
                            trim_levels_beyond=1e-2,
                            linewidth=1,
                            scale=0.03,
                            suppress_invert=False,
                            probability_density=False,
                            wfalpha=0.8,
                            potentialalpha=0.8,

                            **kwargs):
    options = copy.copy(defaults)
    options["square"] = False

    defaults.update(kwargs)
    potentials = schrodinger_result["potentials"]
    wavefunctions = schrodinger_result["wavefunctions"]
    energy_levels = schrodinger_result["E"]
    x = schrodinger_result["x"]

    # This is used when doing the 4x4 kp calculation rather than the single band calculation
    if 'EU' in energy_levels.keys():
        energy_levels['Ehh'] = energy_levels['EU']
        energy_levels['Elh'] = energy_levels['EU']
        wavefunctions["psi_hh"] = wavefunctions["psi_g1"]
        wavefunctions["psi_lh"] = wavefunctions["psi_g2"]

    conduction_data = [GraphData(x, potentials["Ve"], linewidth=2, color="grey", alpha=potentialalpha)]
    valence_data = [
        GraphData(x, potentials["Vlh"], linewidth=2, color="grey", alpha=potentialalpha, dashes=[1, 1], label="Vlh"),
        GraphData(x, potentials["Vhh"], linewidth=2, color="grey", alpha=potentialalpha, label="Vhh")
    ]
    # normalise_psi = lambda p:p/numpy.sqrt(numpy.trapz(x=x*1e9,y=p**2))*q

    conduction_data.extend(
        prepare_wavefunction_data(x, energy_levels["Ee"], wavefunctions["psi_e"], trim_levels_beyond, linewidth, "blue",
                                  scale, suppress_invert, alpha=wfalpha, square=probability_density))
    valence_data.extend(
        prepare_wavefunction_data(x, energy_levels["Ehh"], wavefunctions["psi_hh"], trim_levels_beyond, linewidth,
                                  "green",
                                  scale, suppress_invert, alpha=wfalpha, square=probability_density))
    valence_data.extend(
        prepare_wavefunction_data(x, energy_levels["Elh"], wavefunctions["psi_lh"], trim_levels_beyond, linewidth,
                                  "red",
                                  scale, suppress_invert, alpha=wfalpha, square=probability_density))

    g = Graph(valence_data, subplot=212, **options)
    del options["xlabel"]
    g.add_subplot(conduction_data, subplot=211, **options)
    # g.axis.get_xaxis().set_ticklabels([])
    return g


def split_schrodinger_graph_potentials(schrodinger_result,
                                       trim_levels_beyond=1e-2,
                                       linewidth=1,
                                       scale=0.3,
                                       suppress_invert=False,
                                       probability_density=False,
                                       wfalpha=0.8,
                                       potentialalpha=0.8,
                                       **kwargs):
    defaults = {'step': 0.002, 'margin': 0.02, 'pdf': False, 'show': False, 'dpi': 100, 'fontsize': 12,
                'figsize': (7, 6)}
    options = copy.copy(defaults)
    options["square"] = False

    defaults.update(kwargs)
    potentials = schrodinger_result["potentials"]
    wavefunctions = schrodinger_result["wavefunctions"]
    energy_levels = schrodinger_result["E"]
    x = schrodinger_result["x"]

    # This is used when doing the 4x4 kp calculation rather than the single band calculation
    if 'EU' in energy_levels.keys():
        energy_levels['Ehh'] = energy_levels['EU']
        energy_levels['Elh'] = energy_levels['EU']
        wavefunctions["psi_hh"] = wavefunctions["psi_g1"]
        wavefunctions["psi_lh"] = wavefunctions["psi_g2"]

    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=defaults['figsize'], dpi=defaults['dpi'])

    ax1.plot(x * 1e9, potentials["Ve"] / q, 'k', linewidth=2, label='Ve')
    ax1.set_ylabel('Energy (eV)', fontsize=defaults["fontsize"])
    ax1.tick_params(labelsize=defaults["fontsize"])
    ax1.grid(color='grey', linestyle='--', linewidth=0.5)

    ax2.plot(x * 1e9, potentials["Vlh"] / q, 'k--', linewidth=2, label="Vlh"),
    ax2.plot(x * 1e9, potentials["Vhh"] / q, 'k', linewidth=2, label="Vhh")
    ax2.set_ylabel('Energy (eV)', fontsize=defaults["fontsize"])
    ax2.set_xlabel('Possition (nm)', fontsize=defaults["fontsize"])
    ax2.tick_params(labelsize=defaults["fontsize"])
    ax2.grid(color='grey', linestyle='--', linewidth=0.5)

    e_data = prepare_wavefunction_data_only(x, energy_levels["Ee"] / q, wavefunctions["psi_e"], trim_levels_beyond,
                                            linewidth, "blue",
                                            0.03 / q, suppress_invert, alpha=wfalpha, square=probability_density)
    hh_data = prepare_wavefunction_data_only(x, energy_levels["Ehh"] / q, wavefunctions["psi_hh"], trim_levels_beyond,
                                             linewidth,
                                             "green", 0.03 / q, suppress_invert, alpha=wfalpha,
                                             square=probability_density)
    lh_data = prepare_wavefunction_data_only(x, energy_levels["Elh"] / q, wavefunctions["psi_lh"], trim_levels_beyond,
                                             linewidth,
                                             "red", 0.03 / q, suppress_invert, alpha=wfalpha,
                                             square=probability_density)

    for x, y in e_data:
        ax1.plot(x * 1e9, y, 'blue', linewidth=2, label='e')
    for x, y in hh_data:
        ax2.plot(x * 1e9, y, 'green', linewidth=2, label='hh')
    for x, y in lh_data:
        ax2.plot(x * 1e9, y, 'red', linewidth=2, label='lh')

    if defaults["show"]:
        plt.tight_layout()
        plt.show()

    if defaults["pdf"]:
        handle, path = tempfile.mkstemp(prefix="tmp_solcore_", suffix=".%s" % graph_defaults["format"])
        canvas = FigureCanvasPdf(fig)
        canvas.print_figure(path, dpi=defaults["dpi"], bbox_inches='tight')
        open_with_os(path)


def split_schrodinger_graph_LDOS(schrodinger_result, **kwargs):
    defaults = {'step': 0.001, 'margin': 0.02, 'pdf': False, 'show': False, 'dpi': 100, 'fontsize': 12,
                'figsize': (7, 6)}
    defaults.update(kwargs)

    effective_masses = schrodinger_result["effective_masses"]
    potentials = schrodinger_result["potentials"]
    wavefunctions = schrodinger_result["wavefunctions"]
    energy_levels = schrodinger_result["E"]
    x = schrodinger_result["x"]

    if 'EU' in energy_levels.keys():
        energy_levels['Ehh'] = energy_levels['EU']
        energy_levels['Elh'] = energy_levels['EU']
        wavefunctions["psi_hh"] = wavefunctions["psi_g1"]
        wavefunctions["psi_lh"] = wavefunctions["psi_g2"]

    Ee, LDOSe = LDOS1D_e(x, energy_levels, wavefunctions, effective_masses, defaults['step'], defaults['margin'])
    Eh, LDOSh = LDOS1D_h(x, energy_levels, wavefunctions, effective_masses, defaults['step'], defaults['margin'])

    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=defaults['figsize'], dpi=defaults['dpi'])

    ax1.contourf(x * 1e9, Ee / q, LDOSe, 100, cmap='gnuplot2_r', vmin=0, vmax=max(LDOSe.flatten()) * 1.2)
    ax1.plot(x * 1e9, potentials["Ve"] / q, 'k', linewidth=2, label='Ve')
    ax1.set_ylabel('Energy (eV)', fontsize=defaults["fontsize"])
    ax1.tick_params(labelsize=defaults["fontsize"])

    ax2.contourf(x * 1e9, Eh / q, LDOSh, 100, cmap='gnuplot2_r', vmin=0, vmax=max(LDOSh.flatten()) * 1.2)
    ax2.plot(x * 1e9, potentials["Vlh"] / q, 'k--', linewidth=2, label="Vlh"),
    ax2.plot(x * 1e9, potentials["Vhh"] / q, 'k', linewidth=2, label="Vhh")
    ax2.set_ylabel('Energy (eV)', fontsize=defaults["fontsize"])
    ax2.set_xlabel('Possition (nm)', fontsize=defaults["fontsize"])
    ax2.tick_params(labelsize=defaults["fontsize"])

    if defaults["show"]:
        plt.tight_layout()
        plt.show()

    if defaults["pdf"]:
        handle, path = tempfile.mkstemp(prefix="tmp_solcore_", suffix=".%s" % graph_defaults["format"])
        canvas = FigureCanvasPdf(fig)
        canvas.print_figure(path, dpi=defaults["dpi"], bbox_inches='tight')
        open_with_os(path)

    return Ee, LDOSe, Eh, LDOSh


def LDOS1D_e(x, E, psi, m, step=0.001, margin=0.02, broad=0.005):
    Emax = max(E['Ee']) + margin * q
    Emin = min(E['Ee']) - margin * q

    energy = np.arange(Emin, Emax, step * q)
    LDOS = np.zeros((len(energy), len(x)))

    for i, ee in enumerate(E['Ee']):
        m_plane = calculate_in_plane_masses(x, psi['psi_e'], m['me'])
        LDOS = LDOS + m_plane[i] / pi / hbar ** 2 * np.outer(L(energy, ee, broad * q), psi['psi_e'][i] ** 2)

    return energy, LDOS


def LDOS1D_h(x, E, psi, m, step=0.001, margin=0.02, broad=0.005):
    Emax = max(max(E['Ehh']), max(E['Elh'])) + margin * q
    Emin = min(min(E['Ehh']), min(E['Elh'])) - margin * q

    energy = np.arange(Emin, Emax, step * q)
    LDOS = np.zeros((len(energy), len(x)))

    for i, ee in enumerate(E['Ehh']):
        m_plane = calculate_in_plane_masses(x, psi['psi_hh'], m['mhh'])
        LDOS = LDOS + m_plane[i] / pi / hbar ** 2 * np.outer(L(energy, ee, broad * q), psi['psi_hh'][i] ** 2)

    for i, ee in enumerate(E['Elh']):
        m_plane = calculate_in_plane_masses(x, psi['psi_lh'], m['mlh'])
        LDOS = LDOS + m_plane[i] / pi / hbar ** 2 * np.outer(L(energy, ee, broad * q), psi['psi_lh'][i] ** 2)

    return energy, LDOS


def calculate_in_plane_masses(x, psi, m):
    """ Calculates the in-plane effective mass for each level, considering that the wavefunction leaks into the barriers."""

    m_out = []
    for ps in psi:
        m_out.append(np.trapz(ps ** 2 * m, x))

    return m_out


def prepare_wavefunction_data_only(x, E, Psi, trim_levels_beyond, linewidth, color, scale, suppress_invert=False,
                                   square=False, alpha=1, label="untitled"):
    data = []
    for e, psi in zip(E, Psi):
        norm_psi = normalise_psi(psi if not square else psi ** 2)
        trim_to = norm_psi ** 2 / q ** 2 > trim_levels_beyond ** 2
        trimmed_x = x[trim_to]
        norm_psi = norm_psi[trim_to]
        invert = -1 if norm_psi[0] < 0 and not suppress_invert and not square else 1

        data.append((trimmed_x, e + norm_psi * scale * invert))

    return data


"""


def LDOS_h(x, E, psi, m, step=0.002, margin=0.05):
    Emax = max(max(E['Ehh']), max(E['Elh'])) + margin * q
    Emin = min(min(E['Ehh']), min(E['Elh'])) - margin * q

    energy = np.arange(Emin, Emax, step * q)
    LDOS = np.zeros((len(energy), len(x)))

    for i, ee in enumerate(E['Ehh']):
        m_plane = calculate_in_plane_masses(x, psi['psi_hh'], m['mhh'])
        LDOS = LDOS + m_plane[i] / pi / hbar ** 2 * np.outer((energy <= ee), psi['psi_hh'][i] ** 2)

    for i, ee in enumerate(E['Elh']):
        m_plane = calculate_in_plane_masses(x, psi['psi_lh'], m['mlh'])
        LDOS = LDOS + m_plane[i] / pi / hbar ** 2 * np.outer((energy <= ee), psi['psi_lh'][i] ** 2)

    return energy, LDOS

def LDOS_e(x, E, psi, m, step=0.002, margin=0.05):
    Emax = max(E['Ee']) + margin * q
    Emin = min(E['Ee']) - margin * q

    energy = np.arange(Emin, Emax, step * q)
    LDOS = np.zeros((len(energy), len(x)))

    for i, ee in enumerate(E['Ee']):
        m_plane = calculate_in_plane_masses(x, psi['psi_e'], m['me'])
        LDOS = LDOS + m_plane[i] / pi / hbar ** 2 * np.outer((energy >= ee), psi['psi_e'][i] ** 2)

    # plt.contourf(x, energy/q, DOS, 50, cmap='gnuplot2_r', vmin=0, vmax=max(DOS.flatten())*1.2)
    # plt.clim(0, max(DOS.flatten())*1.2)
    # plt.show()

    return energy, LDOS
    


"""
