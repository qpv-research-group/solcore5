import copy
from solcore.graphing import Graph, GraphData
from solcore.constants import *

defaults = {
    "edit": lambda x, y: (x * 1e9, y / q),
    "xlabel": "Depth (nm)",
    "ylabel": "Energy (eV)",
}


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
