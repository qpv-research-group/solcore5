import numpy as np

from solcore.constants import kb, q, hbar, c
from solcore.structure import Junction
from scipy.optimize import root

from .detailed_balance import iv_detailed_balance


def iv_2diode(junction, options):
    """ Calculates the IV curve of a junction object using the 2-diode equation. All parameters needed for the calculation need to be included in the junction object. Series resistance is included at solar cell level, not at junction level. The junction is then updated with an "iv" function that calculates the IV curve at any voltage.

    :param junction: A junction object.
    :param options: Solver options.
    :return: None.
    """

    T = options.T
    light = options.light_iv
    junction.voltage = options.internal_voltages
    wl = options.wavelength

    # We get some of the minimum parameters
    R_shunt = min(junction.R_shunt, 1e14) if hasattr(junction, 'R_shunt') else 1e14
    n1 = junction.n1 if hasattr(junction, 'n1') else 1
    n2 = junction.n2 if hasattr(junction, 'n2') else 2

    try:
        # We check if we are using a radiative recombination and reference current
        # If that is the case, we solve the properties first using the DB model with the Boltzmann aproximation
        if hasattr(junction, 'reff') and hasattr(junction, 'jref'):
            options.db_mode = 'boltzmann'
            iv_detailed_balance(junction, options)

            reff = max(junction.reff, 1e-10)
            jref = junction.jref
            j01 = junction.j01

            # The reference voltage becomes
            Vref = n1 * kb * T / q * np.log(reff * jref / j01 + 1)

            # Now we can calculate the j02
            junction.j02 = ((1 - reff) * jref - Vref / R_shunt) / (np.exp(q * Vref / (n2 * kb * T)) - 1)
            j02 = junction.j02

        # If not, we are in the normal 2D equation case
        else:
            j01 = junction.j01
            j02 = junction.j02 if hasattr(junction, 'j02') else 0

            # If the saturation currents correspond to a different temperature, we update them for the current temperature.
            if hasattr(junction, 'Tref') and (T != junction.Tref):
                assert hasattr(junction,
                               'Eg'), 'ERROR: The bandgap for each junction (Eg) must be provided if the working ' \
                                      'temperature (T) is different that the reference temperature (Tref). '

                Eg = junction.Eg
                Tref = junction.Tref
                kB = kb / q
                j01 = j01 * (T / Tref) ** 3 * np.exp(-Eg / (n1 * kB) * (1 / T - 1 / Tref))
                j02 = j02 * (T / Tref) ** (5. / 3.) * np.exp(-Eg / (n2 * kB) * (1 / T - 1 / Tref))

        if not light:
            jsc = 0
        else:
            if hasattr(junction, 'jsc'):
                jsc = junction.jsc
            elif hasattr(junction, 'eqe'):
                eqe = junction.eqe
                wl, ph = options.light_source.spectrum(x=wl, output_units='photon_flux_per_m')
                jsc = q * np.trapz(eqe(wl) * ph, wl)
            else:
                jsc = 0

    except AttributeError as err:
        raise AttributeError('ERROR in 2-diode equation. Junction is missing one essential argument. {}'.format(err))

    def iv(v):
        out = j01 * (np.exp(q * v / (n1 * kb * T)) - 1) + j02 * (np.exp(q * v / (n2 * kb * T)) - 1) + v / R_shunt - jsc
        return np.minimum(out, 1e8)

    junction.jsc = jsc
    junction.current = iv(junction.voltage)
    junction.iv = iv


##########################
## Old (v4) MJ calculator.
#
# def multijunctionIV(solar_cell, V=None, photon_recycling=False, n1=1, n2=2, fraction_coupled=None,
#                     coupling_iterations=10):
#     """ Calculates the overall IV characteristics of any number of junctions numerically at the requested voltage
#     points using the 2-diodes model. If photocurrent_np is not provided, the resultant IV characteristics are purely
#     recombination currents, otherwise light IVs are returned.
#
#     If the state_object includes a reference temperature "Tref", then the reverse saturation currents are assumed to be
#     at that temperature and are therefore re-calculated at the target "T". This is only valid for small temperature
#     differences.
#
#     Conventions followed:
#         - SI Units
#         - Photocurrents: Positive.
#         - Dark Currents: Negative
#
#     :param solar_cell: Structure object containing
#
#         - T: overall temperature
#         - Tref: (optional, in case the j01 and j02 correspond to a different temperature that the target temperature)
#         - One or more Junction objects, each with (at least) the following
#         parameters in SI units:
#
#             - j01
#             - j02
#             - n1
#             - n2
#             - R_series
#             - R_shunt
#             - photocurrent (optional, is assumed = 0 if not present)
#             - Eg (optional, bandgap, only used if Tref is given)
#
#     :param photon_recycling: If photon recycling calculation must be included (NOT IMPLEMENTED YET)
#     :param V: Voltages at which to have the resulting IV curve
#     :param n1: Default n1 ideality factor (=1) if not included in the Junction objects
#     :param n2: Default n2 ideality factor (=2) if not included in the Junction objects
#     :param fraction_coupled: fraction of the emitted light coupled from one junction to the next (NOT IMPLEMENTED YET)
#     :param coupling_iterations: the maximum number of iterations during the coupling calculation (NOT IMPLEMENTED YET)
#     :return: Dictionary including:
#
#         - "IV": (V, I) Calculated IV characteristics
#         - "junction IV": [(V junc 1, I junc 1), (V junc 2, I junc 2), ...]
#         - "Rseries IV": (V, I) Calculated IV characteristics of the series resistance
#         - "V": Device voltages
#         - "I": Device currents (the same for all junctions).
#         - "Isc", Voc", "P" and "FF": In case of calculation under illumination.
#     """
#
#     # First we have to prepare the arrays necessary for the calculation
#     junctions = [obj for obj in solar_cell if type(obj) == Junction]
#
#     if hasattr(solar_cell, 'T'):
#         T = solar_cell.T
#     else:
#         T = 298
#
#     # If the saturation currents correspond to a different temperature, we update them for the current temperature.
#     if hasattr(solar_cell, 'Tref'):
#         junctions = update_j0(junctions, T, solar_cell.Tref)
#
#     J01_array = [junc.j01 for junc in junctions]
#     J02_array = [junc.j02 for junc in junctions]
#     n1_array = [junc.n1 if "n1" in junc.__dict__.keys() else n1 for junc in junctions]
#     n2_array = [junc.n2 if "n2" in junc.__dict__.keys() else n2 for junc in junctions]
#     R_shunt_array = [junc.R_shunt for junc in junctions]
#     R_series = sum([junc.R_series for junc in junctions])
#
#     # Preparation in case there is photon recycling
#     try:
#         assert photon_recycling
#         materials = [junc.material for junc in junctions]
#         qes = [junc.qe for junc in junctions]
#         e = solar_cell.energies
#     except AssertionError:
#         materials = None
#         qes = None
#         e = None
#         if photon_recycling:
#             print('WARNING: The material and QE for each junction must be provided in order to include photon '
#                   'recycling effects. Switching-off the photon recycling calculation.')
#             photon_recycling = False
#
#     # Preparation of the photocurrent of each junctions, if pressent
#     try:
#         photocurrent_array = [junc.photocurrent for junc in junctions]
#     except:
#         photocurrent_array = np.zeros_like(J01_array)  # photocurrent=0, ie dark current
#
#     # If there is no input voltage, we create a sensible set of voltages
#     if V is None:
#         output_V = np.linspace(-2, 5, 1000)
#     else:
#         output_V = V
#
#     # And now we perform the actual calculation
#     result = _multijunctionIV(J01_array, J02_array, R_shunt_array, R_series, n1_array, n2_array, T,
#                               photocurrent_array, output_V, photon_recycling, materials, qes, e, fraction_coupled,
#                               coupling_iterations)
#
#     return result
#
#
# def _multijunctionIV(J01_array, J02_array, R_shunt_array, R_series, n1_array, n2_array, T, photocurrent_array,
#                      output_V, photon_recycling=False, materials=None, qes=None, energ=None, fraction_coupled=None,
#                      coupling_iterations=10):
#     """ Calculates the overall IV characteristics of any number of junctions numerically at the requested voltage
#     points. If photocurrent is not provided, the resultant IV characteristics are purely recombination currents,
#     otherwise light IVs are returned.
#
#     Conventions followed:
#         - SI Units
#         - Photocurrents: Positive.
#         - Dark Currents: Negative
#
#     :param J01_array: Array of J01 saturation currents
#     :param J02_array: Array of J02 saturation currents
#     :param R_shunt_array: Array of shunt resistances
#     :param R_series: Combined series resistance of all juncitons
#     :param n1_array: Array of n1 ideality factors
#     :param n2_array: Array of n2 ideality factors
#     :param T: Cell temperature
#     :param photocurrent_array: Array of photocurents
#     :param photon_recycling: (False) if calculation should include photon recycling
#     :param output_V: Array of voltages in which to calculate the data
#     :param materials: Array of materials of all the junctions
#     :param qes: Array of the quantum efficiencies of all the junctions
#     :param energ: Energies of the qes
#     :param fraction_coupled: Array with the fraction of the emitted light coupled from one junction to the next.
#     :param coupling_iterations: the maximum number of iterations during the coupling calculation
#     :return: dictionary with several entries:
#         "IV": (V, I) Calculated IV characteristics
#         "junction IV": [(V junc 1, I junc 1), (V junc 2, I junc 2), ...]
#         "Rseries IV": (V, I) Calculated IV characteristics of the series resistance
#         "V": Device voltages
#         "I": Device currents (the same for all junctions).
#         "Isc", Voc", "P" and "FF": In case of calculation under illumination.
#     """
#
#     # If there is no series resistance, we make it a really small number, rather than zero, for simplicity.
#     Rs = max(np.sqrt(np.finfo(float).eps), R_series)
#
#     # The current and junction voltage arrays
#     num_jun = len(J01_array)
#     V_junction_array = np.zeros((len(output_V), len(J01_array)))
#     output_J = np.zeros_like(output_V)
#
#     # The initial guess for the numerical calculation
#     vals = np.arange(num_jun, 0.5, -1)
#     total = np.sum(vals)
#     guess = output_V[0] * vals / total
#
#     # This is the function we want to minimize: the node equations
#     def F(Vext, voltages):
#
#         # Create the output
#         output = np.zeros_like(voltages)
#
#         # First equation
#         I1 = _current(voltages[0], J01_array[0], J02_array[0], R_shunt_array[0], n1_array[0], n2_array[0],
#                       photocurrent_array[0], T)
#         output[0] = (Vext - np.sum(voltages)) / Rs - I1
#
#         # The rest of the equations
#         for j in range(1, num_jun):
#             I2 = _current(voltages[j], J01_array[j], J02_array[j], R_shunt_array[j], n1_array[j], n2_array[j],
#                           photocurrent_array[j], T)
#             output[j] = I1 - I2
#
#             I1 = I2
#
#         return output
#
#     # Now we loop over all the voltages
#     for i, VV in enumerate(output_V):
#         # We particularise the function to minimize for the current output voltage
#         def FF(voltages):
#             return F(VV, voltages)
#
#         # Solve the problem
#         result = root(FF, guess, method='lm', tol=1e-10)
#         V_junction_array[i] = result['x']
#
#         # Calculate the current and create the next initial guess
#         output_J[i] = _current(V_junction_array[i, 0], J01_array[0], J02_array[0], R_shunt_array[0], n1_array[0],
#                                n2_array[0], photocurrent_array[0], T)
#         guess = V_junction_array[i]
#
#     # Finally, we calculate the solar cell parameters
#     Isc = None
#     Voc = None
#     Pmpp = None
#     Vmpp = None
#     Impp = None
#     FF = None
#
#     # If we are calculating the light IV, we also calculate the main parameters: Jsc, Voc, FF, MPP...
#     if min(photocurrent_array) > 0:
#         Isc = -np.interp([0], output_V, output_J)[0]
#         Voc = np.interp([0], output_J, output_V)[0]
#         Power = abs(output_J[output_V < Voc] * output_V[output_V < Voc])
#         Pmpp = np.max(Power)
#         idx = np.where(Power == Pmpp)
#         Vmpp = output_V[output_V < Voc][idx]
#         Impp = output_J[output_V < Voc][idx]
#         FF = Pmpp / (Isc * Voc)
#
#     return {
#         "IV": (output_V, -output_J),
#         "V": output_V,
#         "I": -output_J,
#         "junction IV": [(V_junction_array[:, i], -output_J) for i in range(num_jun)],
#         "Rseries IV": (output_J * R_series, -output_J),
#         "Isc": Isc,
#         "Voc": Voc,
#         "FF": FF,
#         "Pmpp": Pmpp,
#         "Vmpp": Vmpp,
#         "Impp": Impp
#     }
#
#
# def MJ_IV_calculator(solar_cell, output_V, mpp_parameters=False, photon_recycling=False, materials=None, qes=None,
#                      energy=None, fraction_coupled=None, coupling_iterations=10):
#     """ Calculates the overall IV characteristics of any number of junctions numerically at the requested voltage
#     points. If photocurrent is not provided, the resultant IV characteristics are purely recombination currents,
#     otherwise light IVs are returned.
#
#     Conventions followed:
#         - SI Units
#         - Photocurrents: Positive.
#         - Dark Currents: Negative
#
#     :param solar_cell: A solar cell object with one or more junctions. The IV of the individual junctions must have been calculated.
#     :param mpp_parameters: If Isc, Voc, FF, Vmpp, Impp and Pmpp must be calculated.
#     :param photon_recycling: (False) if calculation should include photon recycling
#     :param output_V: Array of voltages in which to calculate the data
#     :param materials: Array of materials of all the junctions
#     :param qes: Array of the quantum efficiencies of all the junctions
#     :param energy: Energies of the qes
#     :param fraction_coupled: Array with the fraction of the emitted light coupled from one junction to the next.
#     :param coupling_iterations: the maximum number of iterations during the coupling calculation
#     :return: dictionary with several entries:
#         "IV": (V, I) Calculated IV characteristics
#         "junction IV": [(V junc 1, I junc 1), (V junc 2, I junc 2), ...]
#         "Rseries IV": (V, I) Calculated IV characteristics of the series resistance
#         "V": Device voltages
#         "I": Device currents (the same for all junctions).
#         "Isc", Voc", "P" and "FF": In case of calculation under illumination.
#     """
#
#     # If there is no series resistance, we make it a really small number, rather than zero, for simplicity.
#     Rs = max(np.sqrt(np.finfo(float).eps), solar_cell.R_series)
#
#     # The current and junction voltage arrays
#     num_jun = solar_cell.junctions
#     V_junction_array = np.zeros((len(output_V), num_jun))
#     output_J = np.zeros_like(output_V)
#
#     # The initial guess for the numerical calculation
#     vals = np.arange(num_jun, 0.5, -1)
#     total = np.sum(vals)
#     guess = output_V[0] * vals / total
#
#     # This is the function we want to minimize: the node equations
#     def F(voltages, Vext):
#
#         # Create the output
#         output = np.zeros_like(voltages)
#
#         # First equation
#         I1 = solar_cell(0).iv(voltages[0])
#         output[0] = (Vext - np.sum(voltages)) / Rs - I1
#
#         # The rest of the equations
#         for j in range(1, num_jun):
#             I2 = solar_cell(j).iv(voltages[j])
#             output[j] = I1 - I2
#
#             I1 = I2
#
#         return output
#
#     # Now we loop over all the voltages
#     for i, VV in enumerate(output_V):
#         # Solve the problem
#         result = root(F, guess, args=(VV), method='lm', tol=1e-10)
#         V_junction_array[i] = result['x']
#
#         # Calculate the current and create the next initial guess
#         output_J[i] = solar_cell(0).iv(V_junction_array[i, 0])
#         guess = V_junction_array[i]
#
#     # Finally, we calculate the solar cell parameters
#     Isc = None
#     Voc = None
#     Pmpp = None
#     Vmpp = None
#     Impp = None
#     FF = None
#
#     # If we are calculating the light IV, we also calculate the main parameters: Jsc, Voc, FF, MPP...
#     if mpp_parameters:
#         try:
#             Isc = -np.interp([0], output_V, output_J)[0]
#             Voc = np.interp([0], output_J, output_V)[0]
#             Power = abs(output_J[output_V < Voc] * output_V[output_V < Voc])
#             Pmpp = np.max(Power)
#             idx = np.where(Power == Pmpp)
#             Vmpp = output_V[output_V < Voc][idx]
#             Impp = output_J[output_V < Voc][idx]
#             FF = Pmpp / (Isc * Voc)
#         except Exception as err:
#             print('Error calculating the MPP parameters: {}'.format(err))
#
#     return {
#         "IV": (output_V, -output_J),
#         "V": output_V,
#         "I": -output_J,
#         "junction IV": [(V_junction_array[:, i], -output_J) for i in range(num_jun)],
#         "Rseries IV": (output_J * Rs, -output_J),
#         "Isc": Isc,
#         "Voc": Voc,
#         "FF": FF,
#         "Pmpp": Pmpp,
#         "Vmpp": Vmpp,
#         "Impp": Impp
#     }
#
#
# def simple_gp(e, n, T, V, absorption):
#     result = n ** 2 * absorption * e ** 2 / (c ** 2 * np.exp((e - q * V) / (kb * T)) - 1) / hbar ** 2
#     return result
#
#
# def _current(V_junction, J01, J02, Rshunt, n1, n2, junction_photocurrent, Temp):
#     """ Calculates the current flowing through the solar cell using the 2-diode equation. Series resistance is included
#     externally, in the nodes equation.
#
#     :param V_junction: The voltage at which to calculate the current
#     :param J01: Reverse saturation current J01, typically the radiative component
#     :param J02: Reverse saturation current J02, tipically the non-radiative component
#     :param Rshunt: Shunt resistance
#     :param n1: Ideallity factor associated to J01, typically ~1
#     :param n2: Ideallity factor associated to J02, typically ~2
#     :param junction_photocurrent: Photocurrent
#     :param Temp: Junction temperature
#     :return: The current in the above conditions
#     """
#     output = J01 * (np.exp(q * V_junction / (n1 * kb * Temp)) - 1) + \
#              J02 * (np.exp(q * V_junction / (n2 * kb * Temp)) - 1) + \
#              V_junction / Rshunt - junction_photocurrent
#     return output
#

def calculate_J01(Eg_in_eV, T, n):
    """ Calculate the reverse saturation current J01, assumed radiative, considering an absorption equal to 1 above the
    bandgap. Light trapping is included by considering the refractive index of the material:

    .. math:: J_{01} = \\frac {q n^2 k_b T} {2 \\pi ^2 c^2 \\hbar ^3} e^{\\frac{-E_g}{k_b T}} (E_g^2 + 2 k_b T E_g + 2 k_b^2 T^2)


    :param Eg_in_eV: Bandgap  in eV
    :param T: Cell temperature
    :param n: Refractive index of the material
    :return: The reverse saturation current J01
    """
    Eg = Eg_in_eV * q
    Term1 = 2 * np.pi * n ** 2 * q / (4 * np.pi ** 3 * hbar ** 3 * c ** 2)
    Term2 = kb * T * np.exp(-Eg / (kb * T))
    Term3 = Eg ** 2 + (2 * Eg * kb * T) + (2 * kb ** 2 * T ** 2)

    J01 = Term1 * Term2 * Term3
    return J01


def calculate_J02_from_Voc(J01, Jsc, Voc, T, R_shunt=1e15):
    """ Calculates J02 based on the J01, Jsc and the Voc. It is just the result of solving the 2-diode equation for J02.
    Ideality factors n1 and n2 are assumed to be equal to 1 and 2, respectively.

    :param J01: Reverse saturation current J01, typically the radiative component
    :param Jsc: Short circuit current (=photocurrent)
    :param Voc: Open circuit voltage
    :param T: Temperature
    :param R_shunt: Shunt resistance (default = 1e15)
    :return: The reverse saturation current J02
    """
    Term1 = Jsc - J01 * (np.exp(q * Voc / (kb * T)) - 1) - Voc / R_shunt
    Term2 = np.exp(q * Voc / (2 * kb * T)) - 1

    J02 = Term1 / Term2
    return J02


def calculate_J02_from_rad_eff(J01, radiative_efficiency, V, T, R_shunt=1e15):
    """ Calculates J02 based on J01 and a radiative efficiency at a given voltage and temperature. Ideality factors n1
    and n2 are assumed to be equal to 1 and 2, respectively.

    :param J01: Reverse saturation current J01, typically the radiative component
    :param radiative_efficiency: Fraction of the dark current that is radiative
    :param V: Operating voltage
    :param T: Temperature
    :param R_shunt: Shunt resistance (default = 1e15)
    :return:
    """

    Term1 = J01 * (np.exp(q * V / (kb * T)) - 1)
    Term2 = 1 / radiative_efficiency - 1
    Term3 = np.exp(q * V / (2 * kb * T)) - 1

    J02 = (Term1 * Term2 + V / R_shunt) / Term3

    return J02


def calculate_j02_from_J01_Jsc_reference_radiative_efficiency(J01, Jsc_ref, radiative_efficiency, T):
    Voc_reference = (kb * T / q) * np.log((Jsc_ref * radiative_efficiency / J01) + 1)
    J02 = (Jsc_ref - (radiative_efficiency * Jsc_ref)) / (np.exp(q * Voc_reference / (2 * kb * T)) - 1)

    return J02


def update_j0(junctions, T, Tref):
    """ Updates the reverse saturation currents for the target temperature knowing their values at a reference
    temperature.

    :param junctions: List of junctions
    :param T: Target temperature
    :param Tref: Working temperature
    :return: List of junctions with their saturation currents updated for the new temperatures
    """
    kB = kb / q

    for junc in junctions:
        assert hasattr(junc, 'Eg'), 'ERROR: The bandgap for each junction (Eg) must be provided if the working ' \
                                    'temperature (T) is different that the reference temperature (Tref). '
        junc.j01 = junc.j01 * (T / Tref) ** 3 * np.exp(-junc.Eg / kB * (1 / T - 1 / Tref))
        junc.j02 = junc.j02 * (T / Tref) ** (5. / 3.) * np.exp(-0.5 * junc.Eg / kB * (1 / T - 1 / Tref))

    return junctions
