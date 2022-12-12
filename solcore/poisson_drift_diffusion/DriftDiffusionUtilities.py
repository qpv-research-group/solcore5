from typing import Dict, Tuple

import numpy as np
from numpy.typing import NDArray
from scipy.interpolate import interp1d

from .. import asUnit, constants
from ..light_source import LightSource
from ..registries import (
    register_equilibrium_solver,
    register_iv_solver,
    register_short_circuit_solver,
)
from ..state import State
from ..structure import Junction
from .DeviceStructure import CalculateAbsorptionProfile, CreateDeviceStructure

try:
    from .ddModel import driftdiffusion as dd

    reason_to_exclude = None
except ImportError:
    reason_to_exclude = (
        "The built-in Poisson-Drift-Difussion solver could not be found."
    )


Epsi0 = constants.vacuum_permittivity
q = constants.q
pi = constants.pi
h = constants.h
kb = constants.kb
m0 = constants.electron_mass

log = dd.log_file

pdd_options = State()

# Mesh control
pdd_options.meshpoints = -400
pdd_options.growth_rate = 0.7
pdd_options.coarse = 20e-9
pdd_options.fine = 1e-9
pdd_options.ultrafine = 0.2e-9

# Convergence control
pdd_options.clamp = 20
pdd_options.nitermax = 100
pdd_options.ATol = 1e-14
pdd_options.RTol = 1e-6

# Recombination control
pdd_options.srh = 1
pdd_options.rad = 1
pdd_options.aug = 0
pdd_options.sur = 1
pdd_options.gen = 0

# Output control
pdd_options.output_equilibrium = 1
pdd_options.output_sc = 1
pdd_options.output_iv = 1
pdd_options.output_qe = 1


def process_structure(
    junction: Junction, T: float = 298.0, meshpoints: int = -400, **options
) -> dict:
    """Dump the structure information to the Fotran module and initialise the structure.

    This function extracts the electrical and optical properties of the materials, and
    loads all that information into the Fortran variables. Finally, it initialises the
    device (in Fortran) calculating an initial mesh and all the properties as a function
    of the position.

    Args:
        junction: A junction object with layers.
        T (float, optional): Temperature of the cell. Defaults to 298.0.
        meshpoints (int, optional): Number of mesh points. If negative, that sets the
        initial number of points, but the calculation will dynamically adapt those.
            Defaults to -400.

    Returns:
        Dictionary with a "Properties" key containing the device structure properties as
        a function of the position. Additionally, the state of the fortran solver is
        updated with the structural information of the cell and the boundary conditions.
    """
    options = State(**options)

    SetConvergenceParameters(options)
    SetMeshParameters(options)
    SetRecombinationParameters(options)

    device = CreateDeviceStructure("Junction", T=T, layers=junction)

    print("Processing structure...")
    # First, we clean any previous data from the Fortran code
    dd.reset()
    output = {}

    i = 0
    while i < device["numlayers"]:
        layer = device["layers"][i]["properties"]
        args_list = [
            layer["width"],
            asUnit(layer["band_gap"], "eV"),
            asUnit(layer["electron_affinity"], "eV"),
            layer["electron_mobility"],
            layer["hole_mobility"],
            layer["Nc"],
            layer["Nv"],
            layer["electron_minority_lifetime"],
            layer["hole_minority_lifetime"],
            layer["relative_permittivity"] * Epsi0,
            layer["radiative_recombination"],
            layer["electron_auger_recombination"],
            layer["hole_auger_recombination"],
            layer["Na"],
            layer["Nd"],
        ]
        dd.addlayer(args=args_list)
        i = i + device["layers"][i]["numlayers"]

    # We set the surface recombination velocities.
    # This needs to be improved at some point to consider other boundary conditions
    dd.frontboundary("ohmic", device["sn"], device["sp"], 0)
    dd.backboundary("ohmic", device["sn"], device["sp"], 0)

    dd.initdevice(meshpoints)
    print("...done!\n")

    output["Properties"] = DumpInputProperties()
    return output


@register_equilibrium_solver("PDD", reason_to_exclude=reason_to_exclude)
def equilibrium_pdd(
    junction: Junction,
    T: float = 298.0,
    output_equilibrium: int = 1,
    meshpoints: int = -400,
    **options
) -> None:
    """Solves the PDD equations under equilibrium

    That is, in the dark with no external current and zero applied voltage.

    Args:
        junction: A junction object with layers.
        T (float, optional): Temperature of the cell. Defaults to 298.0.
        output_equilibrium:  If the ouput of the equilibrium calculation process
            should be shown. Defaults to 1.
        meshpoints (int, optional): Number of mesh points. If negative, that sets the
        initial number of points, but the calculation will dynamically adapt those.
            Defaults to -400.

    Return:
        None. The input Junction is updated with a new "equilibrium_data" attribute
        containing a State object with the "Properties" and the "Bandstructure" under
        equilibrium conditions as a function of the position.
    """
    output = process_structure(junction=junction, T=T, meshpoints=meshpoints, **options)
    dd.gen = 0

    print("Solving equilibrium...")
    dd.equilibrium(output_equilibrium)
    print("...done!\n")

    output["Bandstructure"] = DumpBandStructure()
    junction.equilibrium_data = State(**output)


@register_short_circuit_solver("PDD", reason_to_exclude=reason_to_exclude)
def short_circuit_pdd(
    junction: Junction,
    wavelength: NDArray,
    light_source: LightSource,
    output_sc: int = 1,
    **options
):
    """Solves the electronic properties of the junction at short circuit.

    Internally, it calls Equilibrium to get the initial bandstructure and then increases
    progrssively the intensity of the light - done in the Fortran code - until the
    desired spectrum.

    The junction object is expected to have a `absorbed` method taking as input the
    an array of mesh points in the junction and returning a 2D array showing the
    fraction of light absorbed at each wavelenght and depth position. If that were not
    the case, absorption is set to zero... which defeats the purpose of doing this
    calculation to start with.

    Args:
        junction: A junction object with layers.
        wavelength (NDArray): Array of the wavelenghts to solve the problem for.
        light_source (LightSource): Light source object defining the illumination
            specturm.
        output_sc (int, optional): If the ouput of the shortcirctuit calculation process
            should be shown. Defaults to 1.
    """

    equilibrium_pdd(junction, **options)

    _, ph = light_source.spectrum(x=wavelength, output_units="photon_flux_per_m")

    z = junction.equilibrium_data["Bandstructure"]["x"]
    absorption = junction.absorbed if hasattr(junction, "absorbed") else lambda x: 0
    abs_profile = CalculateAbsorptionProfile(z, wavelength, absorption)

    dd.set_generation(abs_profile)
    dd.gen = 1

    dd.illumination(ph)

    dd.lightsc(output_sc, 1)
    print("...done!\n")

    output = {
        "Bandstructure": DumpBandStructure(),
        "Optics": {"wavelengths": wavelength},
    }
    junction.short_circuit_data = State(**output)


def calculate_iv(
    junction: Junction,
    vlimit: float,
    vstep: float,
    light_iv: bool = False,
    output_iv: int = 1,
    **options
) -> Dict[str, Dict[str, NDArray]]:
    """Launches the actual IV calculation from 0V to the chosen limit.

    Args:
        junction: A junction object with layers.
        vlimit: Limit of the voltage range. Can be positive or negative.
        vstep: Voltage step. Its sign is adjusted internally to match that
        of the vlimit, becoming possitive if vlimit is possitive and negative
        otherwise.
        light_iv: If light IV should be calculated. Defaults to False.
        output_iv: If IV calculation information should be output. Defaults to 1.

    Returns:
        Returns a dictionary with the "Bandstructure" at the last voltage point in the
        calculation (a dictionary itself) and the "IV" data (another dictionary).
    """
    vstep = vlimit / abs(vlimit) * abs(vstep)
    if light_iv:
        short_circuit_pdd(junction, **options)
    else:
        equilibrium_pdd(junction, **options)

    dd.runiv(vlimit, vstep, output_iv, 0)
    return {"Bandstructure": DumpBandStructure(), "IV": DumpIV()}


@register_iv_solver("PDD", reason_to_exclude=reason_to_exclude)
def iv_pdd(
    junction: Junction,
    internal_voltages: NDArray,
    T: float = 298.0,
    light_iv: bool = False,
    output_iv: int = 1,
    **options
):
    """Calculates the IV curve of the device between 0 V and a given voltage.

    Depending on the options, the IV will be calculated in the dark (calling the
    equilibrium_pdd function) or under illumination (calling the short_circuit_pdd
    function). If the voltage range has positive and negative values, the problem is
    solved twice: from 0 V to the maximum positive and from 0 V to the maximum negative,
    concatenating the results afterwards.

    Args:
        junction: A junction object with layers.
        internal_voltages (NDArray): Array of internal voltages used for the calculation
        of the IV curve of the junction.
        T: Temperature of the cell. Defaults to 298.0.
        light_iv: If light IV should be calculated. Defaults to False.
        output_iv: If IV calculation information should be output. Defaults to 1.

    Returns:
        None. The junction object is updated with the follwing attributes:

            - voltage: Array with the internal voltages (ideltical to the input).
            - current: Array with the current at the internal voltages.
            - recombination_currents: State object with the recombination currents at
              the internal voltages (Jsrh, Jrad, Jaug, Jsur)
            - pdd_data: State object with the raw "positive_V" and "negative_V" data
              produced byt the Fotran solver.
            - iv: A method that takes as input a voltage and produce the correspoinding
              current by interpolating the voltage and current above.

        This is in addition to the attributes added by the equilibrium and the short
        circuit calculations, as relevant. See the description of those functions.
    """
    junction.voltage = internal_voltages
    R_shunt = min(junction.R_shunt, 1e14) if hasattr(junction, "R_shunt") else 1e14

    # Find the appropriate voltage limits
    min_bandgap = find_minimum_bandgap(junction)
    p_on_n = True if junction[0].material.Na >= junction[0].material.Nd else False
    vmax, vmin = find_voltage_limits(
        min_bandgap, max(junction.voltage), min(junction.voltage), p_on_n, T
    )
    vstep = junction.voltage[1] - junction.voltage[0]

    # Run the calculation for the positive and negative voltage ranges separately
    print("Solving IV...")
    output_pos = (
        calculate_iv(junction, vmax, vstep, light_iv, output_iv, **options)
        if vmax > 0
        else {}
    )
    output_neg = (
        calculate_iv(junction, vmin, vstep, light_iv, output_iv, **options)
        if vmin < 0
        else {}
    )
    print("...done!\n")

    # Now we need to put together the data for the positive and negative regions.
    total = consolidate_iv(output_pos.get("IV", {}), output_neg.get("IV", {}))

    # Finally, we calculate the currents at the desired voltages and update the junction
    junction.pdd_data = State({"positive_V": output_pos, "negative_V": output_neg})
    junction.current = (
        np.interp(junction.voltage, total["V"], total["J"]) + junction.voltage / R_shunt
    )
    junction.recombination_currents = State(
        {
            curr: np.interp(junction.voltage, total["V"], total[curr])
            for curr in ("Jrad", "Jsrh", "Jaug", "Jsur")
        }
    )

    # And an interpolator
    junction.iv = interp1d(
        junction.voltage,
        junction.current,
        kind="linear",
        bounds_error=False,
        assume_sorted=True,
        fill_value=(junction.current[0], junction.current[-1]),
    )


def qe_pdd(junction, options):
    """Calculates the quantum efficiency of the device at short circuit. Internally it
    calls ShortCircuit

    :param junction: A junction object
    :param options: Options to be passed to the solver
    :return: None
    """
    print("Solving quantum efficiency...")
    output_info = options.output_qe

    short_circuit_pdd(junction, **options)

    dd.runiqe(output_info)
    print("...done!\n")

    output = {}
    output["QE"] = DumpQE()
    output["QE"]["WL"] = options.wavelength

    z = junction.short_circuit_data["Bandstructure"]["x"]
    absorbed_per_wl = np.trapz(junction.absorbed(z), z, axis=0)

    # This is redundant but useful to keep the same format than the other solvers
    junction.qe = State(**output["QE"])
    junction.qe["IQE"] = junction.qe["EQE"] / np.maximum(absorbed_per_wl, 0.00001)

    # The EQE is actually the IQE inside the fortran solver due to an error in the
    # naming --> to be changed
    junction.eqe = interp1d(
        options.wavelength,
        output["QE"]["EQE"],
        kind="linear",
        bounds_error=False,
        assume_sorted=True,
        fill_value=(output["QE"]["EQE"][0], output["QE"]["EQE"][-1]),
    )


# ----
# Functions for dumping data from the fortran variables
def DumpInputProperties():
    output = {}
    output["x"] = dd.get("x")[0 : dd.m + 1]
    output["Xi"] = dd.get("xi")[0 : dd.m + 1]
    output["Eg"] = dd.get("eg")[0 : dd.m + 1]
    output["Nc"] = dd.get("nc")[0 : dd.m + 1]
    output["Nv"] = dd.get("nv")[0 : dd.m + 1]
    output["Nd"] = dd.get("nd")[0 : dd.m + 1]
    output["Na"] = dd.get("na")[0 : dd.m + 1]

    return output


def DumpBandStructure():
    output = {}
    output["x"] = dd.get("x")[0 : dd.m + 1]
    output["n"] = dd.get("n")[0 : dd.m + 1]
    output["p"] = dd.get("p")[0 : dd.m + 1]
    output["ni"] = dd.get("ni")[0 : dd.m + 1]
    output["Rho"] = dd.get("rho")[0 : dd.m + 1]
    output["Efe"] = dd.get("efe")[0 : dd.m + 1]
    output["Efh"] = dd.get("efh")[0 : dd.m + 1]
    output["potential"] = dd.get("psi")[0 : dd.m + 1]
    output["Ec"] = dd.get("ec")[0 : dd.m + 1]
    output["Ev"] = dd.get("ev")[0 : dd.m + 1]
    output["GR"] = dd.get("gr")[0 : dd.m + 1]
    output["G"] = dd.get("g")[0 : dd.m + 1]
    output["Rrad"] = dd.get("rrad")[0 : dd.m + 1]
    output["Rsrh"] = dd.get("rsrh")[0 : dd.m + 1]
    output["Raug"] = dd.get("raug")[0 : dd.m + 1]

    return output


def DumpIV(IV_info=False):
    # Depending of having PN or NP the calculation of the MPP is a bit different. We
    # move everithing to the 1st quadrant and then send it back to normal
    Nd = dd.get("nd")[0 : dd.m + 1][0]
    Na = dd.get("na")[0 : dd.m + 1][0]
    s = Nd > Na

    output = {}
    output["V"] = (-1) ** s * dd.get("volt")[1 : dd.nvolt + 1]
    output["J"] = dd.get("jtot")[1 : dd.nvolt + 1]
    output["Jrad"] = dd.get("jrad")[1 : dd.nvolt + 1]
    output["Jsrh"] = dd.get("jsrh")[1 : dd.nvolt + 1]
    output["Jaug"] = dd.get("jaug")[1 : dd.nvolt + 1]
    output["Jsur"] = dd.get("jsur")[1 : dd.nvolt + 1]

    if IV_info:
        # We calculate the solar cell parameters
        output["Jsc"] = -np.interp(0, output["V"], output["J"])  # dd.get('isc')[0]
        output["Voc"] = np.interp(
            0, output["J"][output["V"] > 0], output["V"][output["V"] > 0]
        )  # dd.get('voc')[0]
        print(output["Voc"])

        maxPP = np.argmin(output["V"] * output["J"])
        Vmax = output["V"][maxPP - 3 : maxPP + 3]
        Imax = output["J"][maxPP - 3 : maxPP + 3]
        Pmax = Vmax * Imax

        poly = np.polyfit(Vmax, Pmax, 2)
        output["Vmpp"] = -poly[1] / (2 * poly[0])
        output["Jmpp"] = -output["J"][np.argmin(np.abs(output["V"] - output["Vmpp"]))]
        output["FF"] = output["Vmpp"] * output["Jmpp"] / (output["Voc"] * output["Jsc"])

        print("Jsc  = %5.3f mA/cm2" % (output["Jsc"] / 10))
        print("Voc  = %4.3f  V" % (output["Voc"]))
        print("FF   = %3.3f " % (output["FF"]))
        print("Jmpp = %5.3f mA/cm2" % (output["Jmpp"] / 10))
        print("Vmpp = %4.3f  V" % (output["Vmpp"]))
        print("Power= %5.3f mW/cm2" % (output["Jmpp"] * output["Vmpp"] / 10))

    # If NP, V and J should be in the 3rd quadrant
    output["V"] = (-1) ** s * output["V"]
    output["J"] = (-1) ** s * output["J"]
    output["Jrad"] = (-1) ** s * output["Jrad"]
    output["Jsrh"] = (-1) ** s * output["Jsrh"]
    output["Jaug"] = (-1) ** s * output["Jaug"]
    output["Jsur"] = (-1) ** s * output["Jsur"]

    return output


def DumpQE():
    output = {}
    numwl = dd.numwl + 1
    output["EQE"] = dd.get("iqe")[0:numwl]
    output["EQEsrh"] = dd.get("iqesrh")[0:numwl]
    output["EQErad"] = dd.get("iqerad")[0:numwl]
    output["EQEaug"] = dd.get("iqeaug")[0:numwl]
    output["EQEsurf"] = dd.get("iqesurf")[0:numwl]
    output["EQEsurb"] = dd.get("iqesurb")[0:numwl]

    return output


# ----
# Functions for setting the parameters controling the recombination, meshing and the
# numerial algorithm
def SetMeshParameters(options):
    dd.set("coarse", options.coarse)
    dd.set("fine", options.fine)
    dd.set("ultrafine", options.ultrafine)
    dd.set("growth", options.growth_rate)


def SetRecombinationParameters(options):
    dd.srh = options.srh
    dd.rad = options.rad
    dd.aug = options.aug
    dd.sur = options.sur
    dd.gen = options.gen


def SetConvergenceParameters(options):
    dd.nitermax = options.nitermax
    dd.set("clamp", options.clamp)
    dd.set("atol", options.ATol)
    dd.set("rtol", options.RTol)


def find_minimum_bandgap(junction: Junction) -> float:
    """Finds the minimum bandgap of the stack of layers of the junction.

    Args:
        junction: The junction to find the minimum bandgap for.

    Returns:
        The minimum banfdgap
    """
    return min([layer.material.band_gap / q for layer in junction])


def find_voltage_limits(
    bandgap: float, vmax: float, vmin: float, p_on_n: bool, T: float
) -> Tuple[float, float]:
    """Finds the voltage limits based on a given bandgap.

    This is needed to ensure that the IV calculation doest not go to a too high
    voltage which could drive the material into the high injection regime, something
    the PDD solver cannot handle. So, the maxmum/minimum voltages allowed are ±3kbT
    with respect to the given bandgap or vmax/vmin if they are within those limits.

    FIXME: It seems that the current implementation actually goes above the high
    injection limit. The signs should be reversed such that the voltage is never higher,
    in absolute terms, than the bandgap.

    Args:
        bandgap: The reference bandgap to use when considering high injection.
        vmax: The desired maximum voltage.
        vmin: The desired minimum voltage.
        p_on_n: The polarity of the junction. For P-on-N junctions, the maximum
        voltage is the one to be limitted. For N-on-P junctions, it is the minimum
        one.
        T: Temperature of the cell.

    Returns:
        A tuple with the adapted vmax and vmin for the calculation.
    """

    if p_on_n:
        vmax = min(bandgap + 3 * kb * T / q, vmax)
    else:
        vmin = max(-bandgap - 3 * kb * T / q, vmin)

    return vmax, vmin


def consolidate_iv(
    positive: Dict[str, NDArray], negative: Dict[str, NDArray]
) -> Dict[str, NDArray]:
    """Consolidate the positive and negative IV curves in a single array.

    If there is no positive range, the negative is returned; if there is no
    negative, the positive is return; and if both are present, the arrays are
    concatenated.

    Args:
        positive: The positive range of the IV curves.
        negative: The negative range of the IV curves.

    Returns:
        A dictionary with the same keys than the inputs but with the ranges
        concatenaded, if relevant.
    """
    if not positive:
        return negative
    elif not negative:
        return positive
    else:
        return {
            mag: np.concatenate((negative[mag][:0:-1], positive[mag]))
            for mag in negative.keys()
        }
