__author__ = 'diego'

""" Implementation of the mobility model by:

M. Sotoodeh, A. H. Khalid, and A. A. Rezazadeh,
“Empirical low-field mobility model for III–V compounds applicable in device simulation codes,”
J. Appl. Phys., vol. 87, no. 6, p. 2890, 2000.

"""

import json
import numpy as np
import os
from solcore.science_tracker import science_reference

# Constants
kb = 8.6173324e-5  # eV K-1
Log = lambda x: np.log10(x)

this_dir = os.path.split(__file__)[0]
parameters = os.path.join(this_dir, "mobility_parameters.json")
f = open(parameters, mode="r")
data = json.load(f)


def mobility_low_field(N, muMin, muMax, Nref, l, t1, t2, T=300):
    m = muMin + (muMax * (300 / T) ** t1 - muMin) / (1 + (N / (Nref * (300 / T) ** t2)) ** l)

    return m


def calculate_mobility(material, holes, N, x=0.0, y=0.0, T=300):
    """ Calculates the mobility using the model by Sotoodeh et al. If the material is not in the database, then the function returns the mobility for GaAs at that temperature, T, and impurity concentration, N.

    :param material: A string with the material name
    :param holes: If calculation should be done for electrons (holes=0) or holes (holes=1)
    :param N: Impurity concentration
    :param x: The fractional composition in the case of ternaries
    :param y: The other fractional composition in the case of quaternaries
    :param T: Temperature
    :return: The calculated mobility
    """
    science_reference('mobility calculator', 'M. Sotoodeh, A. H. Khalid, and A. A. Rezazadeh,'
                                             '“Empirical low-field mobility model for III–V compounds applicable in device simulation codes,"'
                                             'J. Appl. Phys., vol. 87, no. 6, p. 2890, 2000.')

    i = 1
    if holes: i = 2

    if material not in data.keys():
        print("Error: Material {0} not in the database for the mobility. Reverting to GaAs.".format(material))
        d = data['GaAs'][i]
    elif data[material][0] == 2:
        d = data[material][i]
    elif material == "InGaAs":
        d = calculate_InGaAs(x, i)
    elif material == "GaInP":
        d = calculate_InGaP(x, i, T)
    elif material == "AlGaAs":
        d = calculate_AlGaAs(x, i, T)
    elif material == "InAlAs":
        d = calculate_InAlAs(x, i, T)
    elif material == "InGaAsP":
        d = calculate_InGaAsP(x, y, i, T)
    else:
        d = calculate_General(material, x, i, T)

    muMin = d["muMin"]
    muMax = d["muMax"]
    Nref = d["Nref"]
    l = d["l"]
    t1 = d["t1"]
    t2 = d["t2"]

    m = mobility_low_field(N / 1e6, muMin, muMax, Nref, l, t1, t2, T) / 10000  # To convert it from cm2 to m2
    return m


def calculate_InGaAs(x, i):
    """ Calculates the parameters for an InGaAs alloy.

    :param x: Indium fraction
    :param i: If the data for electrons (1) or holes (2) should be calculated
    :return:
    """

    p0 = data["InAs"][i]
    p1 = data["InGaAs"][i]
    p2 = data["GaAs"][i]

    xi = data["InGaAs"][3]["x"]

    newData = {}
    newData["muMin"] = interpolate_parameter_quad(x, p0["muMin"], p1["muMin"], p2["muMin"], xi)
    newData["muMax"] = interpolate_parameter_quad(x, p0["muMax"], p1["muMax"], p2["muMax"], xi)
    newData["Nref"] = 10 ** (interpolate_parameter_quad(x, Log(p0["Nref"]), Log(p1["Nref"]), Log(p2["Nref"]), xi))
    newData["l"] = interpolate_parameter_quad(x, p0["l"], p1["l"], p2["l"], xi)
    newData["t1"] = interpolate_parameter_quad(x, p0["t1"], p1["t1"], p2["t1"], xi)
    newData["t2"] = interpolate_parameter_quad(x, p0["t2"], p1["t2"], p2["t2"], xi)

    return newData


def calculate_InGaP(x, i, T):
    """ Calculates the parameters for an InGaP alloy.

    :param x: Indium fraction
    :param i: If the data for electrons (1) or holes (2) should be calculated
    :return:
    """

    p0 = data["InP"][i]
    p1 = data["GaInP"][i]
    p2 = data["GaP"][i]

    xi = data["GaInP"][3]["x"]

    newData = {}
    newData["muMin"] = interpolate_parameter_quad(x, p0["muMin"], p1["muMin"], p2["muMin"], xi)
    newData["muMax"] = interpolate_parameter_quad(x, p0["muMax"], p1["muMax"], p2["muMax"], xi)
    newData["Nref"] = 10 ** (interpolate_parameter_quad(x, Log(p0["Nref"]), Log(p1["Nref"]), Log(p2["Nref"]), xi))
    newData["l"] = interpolate_parameter_quad(x, p0["l"], p1["l"], p2["l"], xi)
    newData["t1"] = interpolate_parameter_quad(x, p0["t1"], p1["t1"], p2["t1"], xi)
    newData["t2"] = interpolate_parameter_quad(x, p0["t2"], p1["t2"], p2["t2"], xi)

    if i == 1:
        # We have electrons, the above muMin and muMax need to be recalculated in an smarter way to take into
        # account the indirect bandgap above certain composition

        p0b = data["InP"][3]
        p1b = data["GaInP"][3]
        p2b = data["GaP"][3]

        # We calculate the band parameters for the alloy
        newBand = {}
        for key in p0b.keys():
            if key in ["es", "einf"]:
                newBand[key] = interpolate_epsilon(x, p0b[key], p2b[key])
            elif key in ["mnG", "mnX", "mnL"]:
                newBand[key] = interpolate_parameter_linear(x, p0b[key], p2b[key])
            else:
                newBand[key] = interpolate_parameter_linear(x, p0b[key], p2b[key], ABC=p1b["b{0}".format(key[2])])

        # Now we use these to calculate the direct and indirect mobilities (max and min)
        muDmax = (x - xi) / (1 - xi) * p0["muMax"] + (x - 1) / (xi - 1) * p1["muMax"]
        muDmin = (x - xi) / (1 - xi) * p0["muMin"] + (x - 1) / (xi - 1) * p1["muMin"]

        mind_alloy = mind(newBand["EgX"], newBand["EgL"], newBand["mnX"], newBand["mnL"], T)
        C = (p2b["mnX"] / mind_alloy) ** 1.5 * (p2b["einf"] * (-1) - p2b["es"] * (-1)) / (
        newBand["einf"] * (-1) - newBand["es"] * (-1))

        muImax = C * p2["muMax"]
        muImin = C * p2["muMin"]

        f = Rd(newBand["EgG"], newBand["EgX"], newBand["EgL"], newBand["mnG"], newBand["mnX"], newBand["mnL"], T)

        # Finally
        newData["muMin"] = f * muDmin + (1 - f) * muImin
        newData["muMax"] = f * muDmax + (1 - f) * muImax

    return newData


def calculate_AlGaAs(x, i, T):
    """ Calculates the parameters for an AlGaAs alloy.

    :param x: Al fraction
    :param i: If the data for electrons (1) or holes (2) should be calculated
    :return:
    """

    p0 = data["AlAs"][i]
    p1 = data["AlGaAs"][i]
    p2 = data["GaAs"][i]

    xi = data["AlGaAs"][3]["x"]

    newData = {}
    if i == 2:
        # We have holes, which are easy: quadratic interpolation between AlAs, Al0.3GaAs and GaAs.
        newData["muMin"] = interpolate_parameter_quad(x, p0["muMin"], p1["muMin"], p2["muMin"], xi)
        newData["muMax"] = interpolate_parameter_quad(x, p0["muMax"], p1["muMax"], p2["muMax"], xi)
        newData["Nref"] = 10 ** (interpolate_parameter_quad(x, Log(p0["Nref"]), Log(p1["Nref"]), Log(p2["Nref"]), xi))
        newData["l"] = interpolate_parameter_quad(x, p0["l"], p1["l"], p2["l"], xi)
        newData["t1"] = interpolate_parameter_linear(x, p0["t1"], p2["t1"]) / (1 + x * (1 - x))
        newData["t2"] = p2["t2"]

    else:
        # We have electrons. Data is interpolated linearly between AlAs and GaAs
        newData["Nref"] = 10 ** (interpolate_parameter_linear(x, Log(p0["Nref"]), Log(p2["Nref"])))
        newData["l"] = interpolate_parameter_linear(x, p0["l"], p2["l"])
        newData["t1"] = interpolate_parameter_linear(x, p0["t1"], p2["t1"]) / (1 + x * (1 - x))
        newData["t2"] = interpolate_parameter_linear(x, p0["t2"], p2["t2"])

        # For muMin and muMax need to be recalculated in an smarter way to take into
        # account the indirect bandgap above certain composition

        p0b = data["AlAs"][3]
        p1b = data["AlGaAs"][3]
        p2b = data["GaAs"][3]

        # We calculate the band parameters for the alloy
        newBand = {}
        for key in p0b.keys():
            if key in ["es", "einf"]:
                newBand[key] = interpolate_epsilon(x, p0b[key], p2b[key])
            elif key in ["mnG", "mnX", "mnL"]:
                newBand[key] = interpolate_parameter_linear(x, p0b[key], p2b[key])
            else:
                newBand[key] = interpolate_parameter_linear(x, p0b[key], p2b[key], ABC=p1b["b{0}".format(key[2])])

        # Now we use these to calculate the direct and indirect mobilities (max and min)
        C = (p2b["mnG"] / newBand["mnG"]) ** 1.5 * (p2b["einf"] * (-1) - p2b["es"] * (-1)) / (
        newBand["einf"] * (-1) - newBand["es"] * (-1))

        muDmax = C * p2["muMax"]
        muDmin = C * p2["muMin"]

        muImax = p1["muMax"]
        muImin = p1["muMin"]

        f = Rd(newBand["EgG"], newBand["EgX"], newBand["EgL"], newBand["mnG"], newBand["mnX"], newBand["mnL"], T)

        # Finally
        newData["muMin"] = f * muDmin + (1 - f) * muImin
        newData["muMax"] = f * muDmax + (1 - f) * muImax

    return newData


def calculate_InAlAs(x, i, T):
    """ Calculates the parameters for an InAlAs alloy.

    :param x: Al fraction
    :param i: If the data for electrons (1) or holes (2) should be calculated
    :return:
    """

    p0 = data["InAs"][i]
    p1 = data["InAlAs"][i]
    p2 = data["AlAs"][i]

    xi = data["InAlAs"][3]["x"]

    newData = {}
    if i == 2:
        # We have holes, which are easy
        newData["muMin"] = interpolate_parameter_linear(x, p0["muMin"], p2["muMin"]) / (1 + x * (1 - x))
        newData["muMax"] = interpolate_parameter_linear(x, p0["muMax"], p2["muMax"])
        newData["Nref"] = 10 ** (interpolate_parameter_linear(x, Log(p0["Nref"]), Log(p2["Nref"])))
        newData["l"] = interpolate_parameter_linear(x, p0["l"], p2["l"])
        newData["t1"] = interpolate_parameter_linear(x, p0["t1"], p2["t1"]) / (1 + x * (1 - x))
        newData["t2"] = interpolate_parameter_linear(x, p0["t2"], p2["t2"])

    else:
        # We have electrons
        newData["Nref"] = 10 ** (interpolate_parameter_quad(x, Log(p0["Nref"]), Log(p1["Nref"]), Log(p2["Nref"]), xi))
        newData["l"] = interpolate_parameter_quad(x, p0["l"], p1["l"], p2["l"], xi)
        newData["t1"] = interpolate_parameter_linear(x, p0["t1"], p2["t1"]) / (1 + x * (1 - x))
        newData["t2"] = interpolate_parameter_linear(x, p0["t2"], p2["t2"])

        # For muMin and muMax need to be recalculated in an smarter way to take into
        # account the indirect bandgap above certain composition

        p0b = data["InAs"][3]
        p1b = data["InAlAs"][3]
        p2b = data["AlAs"][3]

        # We calculate the band parameters for the alloy
        newBand = {}
        for key in p0b.keys():
            if key in ["es", "einf"]:
                newBand[key] = interpolate_epsilon(x, p0b[key], p2b[key])
            elif key in ["mnG", "mnX", "mnL"]:
                newBand[key] = interpolate_parameter_linear(x, p0b[key], p2b[key])
            else:
                newBand[key] = interpolate_parameter_linear(x, p0b[key], p2b[key], ABC=p1b["b{0}".format(key[2])])

        # Now we use these to calculate the direct and indirect mobilities (max and min)
        f1 = (x - xi) / (1 - xi)
        f2 = (x - 1) / (xi - 1)
        muDmax = p0["muMax"] ** f1 * p1["muMax"] ** f2
        muDmin = p0["muMin"] ** f1 * p1["muMin"] ** f2

        mind_alloy = mind(newBand["EgX"], newBand["EgL"], newBand["mnX"], newBand["mnL"], T)
        C = (p2b["mnX"] / mind_alloy) ** 1.5 * (p2b["einf"] * (-1) - p2b["es"] * (-1)) / (
        newBand["einf"] * (-1) - newBand["es"] * (-1))

        muImax = C * p2["muMax"]
        muImin = C * p2["muMin"]

        f = Rd(newBand["EgG"], newBand["EgX"], newBand["EgL"], newBand["mnG"], newBand["mnX"], newBand["mnL"], T)

        # Finally
        newData["muMin"] = f * muDmin + (1 - f) * muImin
        newData["muMax"] = f * muDmax + (1 - f) * muImax

    return newData


def calculate_InGaAsP(x, y, i, T):
    """ Calculates the parameters for an InGaAsP alloy. The calculation is based on a interpolation scheme between
    InGaP and InGaAs using data of compositions lattice matched to InP. Results for compositions away from this might
    not be very accurate.

    :param x: Indium fraction
    :param y: Phosphorus fraction
    :param i: If the data for electrons (1) or holes (2) should be calculated
    :return:
    """

    p0 = calculate_InGaP(x, i, T)
    p1 = calculate_InGaAs(x, i)

    newData = {}
    newData["muMax"] = interpolate_parameter_linear(y, p0["muMax"], p1["muMax"]) / (1 + 6 * y * (1 - y))
    newData["Nref"] = 10 ** (interpolate_parameter_linear(y, Log(p0["Nref"]), Log(p1["Nref"])))
    newData["l"] = interpolate_parameter_linear(y, p0["l"], p1["l"])
    newData["t1"] = interpolate_parameter_linear(y, p0["t1"], p1["t1"]) / (1 + y * (1 - y))
    newData["t2"] = interpolate_parameter_linear(y, p0["t2"], p1["t2"])

    if i == 2:
        # We have holes
        newData["muMin"] = interpolate_parameter_linear(y, p0["muMin"], p1["muMin"])
    else:
        # We have electrons
        newData["muMin"] = interpolate_parameter_linear(y, p0["muMin"], p1["muMin"]) / (1 + 6 * y * (1 - y))

    return newData


def calculate_General(material, x, i, T):
    """ Calculates the parameters for a general alloy of the materials in the database assuming a simple linear
    interpolation. Only ternaries are supported this way.

    :param material: Material to calculate, which must be in the database
    :param x: Main fraction
    :param i: If the data for electrons (1) or holes (2) should be calculated
    :return:
    """

    parent1 = data[material][4]
    parent2 = data[material][5]

    p0 = data[parent1][i]
    p1 = data[parent2][i]

    newData = {}
    newData["muMin"] = interpolate_parameter_linear(x, p0["muMin"], p1["muMin"])
    newData["muMax"] = interpolate_parameter_linear(x, p0["muMax"], p1["muMax"])
    newData["Nref"] = 10 ** (interpolate_parameter_linear(x, Log(p0["Nref"]), Log(p1["Nref"])))
    newData["l"] = interpolate_parameter_linear(x, p0["l"], p1["l"])
    newData["t1"] = interpolate_parameter_linear(x, p0["t1"], p1["t1"])
    newData["t2"] = interpolate_parameter_linear(x, p0["t2"], p1["t2"])

    return newData


def Rd(EgG, EgX, EgL, mnG, mnX, mnL, T):
    """ Calculates the fraction of electrons in the direct valley.

    :param EgG: Gamma-valley bandgap
    :param EgX: X-valley bandgap
    :param EgL: L-valley bandgap
    :param mnG: Gamma-valley effective mass
    :param mnX: X-valley effective mass
    :param mnL: L-valley effective mass
    :param T: The temperature
    :return: The fraction.
    """

    RX = (mnX / mnG) ** 1.5 * np.exp((EgG - EgX) / (kb * T))
    RL = (mnL / mnG) ** 1.5 * np.exp((EgG - EgL) / (kb * T))

    fraction = 1. / (1 + RX + RL)
    return fraction


def mind(EgX, EgL, mnX, mnL, T):
    RX = np.exp((EgL - EgX) / (kb * T))

    fraction = 1. / (1 + RX)

    return mnL * fraction + mnX * (1 - fraction)


def interpolate_parameter_linear(x, AC, BC, ABC=0):
    return x * AC + (1 - x) * BC - ABC * x * (1 - x)


def interpolate_parameter_quad(x, y0, y1, y2, x1, x0=1, x2=0):
    dem0 = (x0 - x1) * (x0 - x2)
    dem1 = (x1 - x0) * (x1 - x2)
    dem2 = (x2 - x0) * (x2 - x1)

    nom0 = (x - x1) * (x - x2)
    nom1 = (x - x0) * (x - x2)
    nom2 = (x - x0) * (x - x1)

    p = y0 * nom0 / dem0 + y1 * nom1 / dem1 + y2 * nom2 / dem2

    return p


def interpolate_epsilon(x, AC, BC):
    ratioAC = (AC - 1) / (AC + 2)
    ratioBC = (BC - 1) / (BC + 2)

    k = interpolate_parameter_linear(x, ratioAC, ratioBC)

    return (2 * k + 1) / (1 - k)


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    N = np.logspace(14, 21, 10)
    x = np.linspace(0, 1)
    # mo = calculate_mobility("InGaAsP", 0, N=1e17, x=x, y=0.5, T=300)
    mo = calculate_mobility("GaInP", 1, N=1e23, x=x, y=0.0, T=300)
    # mo = calculate_mobility("GaAs", 0, N=N, x=0, T=300)

    # plt.semilogx(N, mo)
    plt.plot(x, mo)
    plt.show()
