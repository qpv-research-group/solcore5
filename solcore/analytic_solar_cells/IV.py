import numpy as np
import math
from scipy.optimize import root

from solcore.state import State
from solcore.constants import kb, q, c, h


def iv_multijunction(solar_cell, options):
    """ Calculates the overall IV characteristics of any number of junctions numerically at the requested voltage
    points. If photocurrent is not provided, the resultant IV characteristics are purely recombination currents,
    otherwise light IVs are returned.

    In the end, the soalr_cell object is updated with an "iv" attribute containing a dictionary with:
        "IV": (V, I) Calculated IV characteristics
        "junction IV": [(V junc 1, I junc 1), (V junc 2, I junc 2), ...]
        "Rseries IV": (V, I) Calculated IV characteristics of the series resistance
        "coupled IV": For each junction but for the first one, a list with the coupled current coming form the upper junctions.
        "Isc", Voc", "P", "FF" and "Eta": In case of mpp = True and light IV.

    The sign convention is:
        - Photocurrents: Positive.
        - Dark Currents: Negative

    :param solar_cell: A solar cell object with one or more junctions. The IV of the individual junctions must have been calculated already.
    :param kwargs: A dictionary containing, at least, the following elements:
        - mpp: (Boolean) If Isc, Voc, FF, Vmpp, Impp and Pmpp must be calculated.
        - voltages: Array of voltages in which to calculate the data
        - light_iv: (Boolean) if light IV is being calculated
    :return: None

    """
    output_V = options.voltages
    mpp_parameters = options.mpp
    radiative_coupling = options.radiative_coupling
    light_iv = options.light_iv

    Rs = max(solar_cell.R_series, 1e-14)

    # The current and junction voltage arrays
    num_jun = solar_cell.junctions
    tunnel_jun = solar_cell.tunnel_indices
    V_junction_array = np.zeros((len(output_V), num_jun))

    # The following assumes that all junctions have the currents defined at the same voltages
    minimum_J = solar_cell(0).current
    for j in range(1, num_jun):
        minimum_J = np.where(np.absolute(minimum_J) < np.absolute(solar_cell(j).current), minimum_J,
                             solar_cell(j).current)

    minimum_J = np.sort(minimum_J)

    # For each junction, we find the voltage corresponding to each of these currents
    temp_V_junction_array = np.zeros((len(minimum_J), num_jun))
    for j in range(num_jun):
        temp_V_junction_array[:, j] = np.interp(minimum_J, solar_cell(j).current, solar_cell(j).voltage)

    # We calculate the total voltage related to the...
    # ... series resistance
    temp_V_total = Rs * minimum_J

    # ... tunnel junctions
    for j in tunnel_jun:
        temp_V_total -= solar_cell[j].vi(-minimum_J)

    # ... and the normal junctions
    for j in range(num_jun):
        temp_V_total += temp_V_junction_array[:, j]

    # Finally, we calculate the current at the requested voltages and the voltages for each junction
    output_J = np.interp(output_V, temp_V_total, minimum_J)
    for j in range(num_jun):
        V_junction_array[:, j] = np.interp(output_V, temp_V_total, temp_V_junction_array[:, j])

    coupled_J = []
    if radiative_coupling:
        # First of all, we check if all the junctions are DB or 2D
        ok = True
        for i in range(num_jun):
            ok = ok and solar_cell(i).kind in ['DB', '2D']

        if ok:
            output_J, V_junction_array, coupled_J = solve_radiative_coupling(solar_cell, options, V_junction_array)
        else:
            print(
                'WARNING: Only "DB" and "2D" junctions (using the DB solver) can use radiative coupling.\nSkipping calculation...')

    # Finally, we calculate the solar cell parameters
    Isc = None
    Voc = None
    Pmpp = None
    Vmpp = None
    Impp = None
    FF = None
    Eta = None

    # If we are calculating the light IV, we also calculate the main parameters: Jsc, Voc, FF, MPP...
    if mpp_parameters and light_iv:
        VV = output_V
        II = -output_J
        PP = VV * II
        power = options.light_source.power_density
        try:
            Isc = np.interp([0], VV, II)[0]
            Voc = np.interp([0], -II, VV)[0]
            Pmpp = np.max(PP)
            idx = np.argmax(PP)
            Vmpp = VV[idx]
            Impp = II[idx]
            FF = Pmpp / (Isc * Voc)
            Eta = Pmpp/power

        except Exception as err:
            print('Error calculating the MPP parameters: {}'.format(err))

    solar_cell.iv = State({
        "IV": (output_V, -output_J),
        "junction IV": [(V_junction_array[:, i], -output_J) for i in range(num_jun)],
        "coupled IV": coupled_J,
        "Rseries IV": (output_J * Rs, -output_J),
        "Isc": Isc,
        "Voc": Voc,
        "FF": FF,
        "Pmpp": Pmpp,
        "Vmpp": Vmpp,
        "Impp": Impp,
        "Eta": Eta
    })


def solve_radiative_coupling(solar_cell, options, V_junction_array):
    """ Calculates the radiative IV curve of a MJ solar cell in the presence of radiative coupling between subcells.

    WARNING: Tunnel junctions are not implemented in the radiative coupling mode, yet.

    :param solar_cell: The MJ solar cell structure
    :param options: General options of the solver
    :param V_junction_array: Array with all the junction voltages without coupling. Only the first element is used, actually.
    :return: A tupple with J, V_junction_array in the presence of coupling and the coupled current
    """
    print('Calculating IV curve with radiative coupling...')
    print('WARNING: Tunnel junctions are not implemented in the radiative coupling calculator, yet.')

    output_V = options.voltages
    wl = options.wavelength
    T = options.T
    num_jun = solar_cell.junctions

    # First we calculate the front surface reflection and effective refractive index, as we will do normally for a
    # detailed balanced calculation
    Rn = np.minimum(0.999999, solar_cell(0).reflected(wl))

    # Out of the normal incidence reflection we calculate an effective refractive index
    neff = (1 + np.sqrt(Rn)) / (1 - np.sqrt(Rn))

    # And the corresponding critical angle
    cos_tc = np.sqrt(1 - 1 / neff ** 2)
    steps = int(math.ceil(20 / (1 - min(cos_tc))))
    cos_t = np.linspace(0.001, 1, steps)

    # And now, we back calculate the reflection at any internal angle
    R = np.ones((len(wl), steps))

    for i in range(len(wl)):
        idx = np.where(cos_t > cos_tc[i])

        Rs = np.ones(steps)
        Rp = np.ones(steps)

        cos_i = np.sqrt(1 - neff[i] ** 2 * (1 - cos_t[idx] ** 2))

        Rs[idx] = ((cos_i - neff[i] * cos_t[idx]) / (cos_i + neff[i] * cos_t[idx])) ** 2
        Rp[idx] = ((cos_t[idx] - neff[i] * cos_i) / (cos_t[idx] + neff[i] * cos_i)) ** 2

        R[i, :] = (Rs + Rp) / 2

    # Now, we loop over the junctions, considering all possible couplings
    coupling_currents = np.zeros((num_jun - 1, num_jun - 1))

    # The emitting cell
    for i in range(num_jun - 1):
        # B = exp(-alpha*w) but it is not necessary to calculate alpha and w explicitly
        Bem = 1 - solar_cell(i).absorptance(wl)

        n = solar_cell(i).n
        j01 = 2 * q * n ** 2 * c

        # The absorbing cell:
        Btrans = 1
        for j in range(i + 1, num_jun):
            Bab = 1 - solar_cell(j).absorptance(wl)

            A = 0.5 * (1 + R[:, -1] * Bem) * (1 - Bem) * (1 - Bab) * Btrans

            for k in range(1, steps - 1):
                t = cos_t[k]

                # The integral of the absorpton/emission across the back surface
                A = A + (1 + R[:, k] * Bem ** (1 / t)) * (
                    1 - Bem ** (1 / t)) * (1 - Bab ** (1 / t)) * Btrans ** (1 / t) * t

            A = 2 * np.pi * A / (steps + 1)
            boltz = j01 * (A / wl ** 4) * np.exp(- (h * c / wl) / (kb * T))

            # Now the coupling factors are defined, we need to calculate the coupling current, integrating over wavlenegth
            # We use the boltzmann approximation to make our life easier (and faster).
            coupling_currents[i, j - 1] = np.trapz(boltz, wl)

            # And we reduce the fraction of light reaching the next junction
            Btrans = Btrans * Bab

    # With this, we can calculate the IV curves numerically
    R_series = max(solar_cell.R_series, 1e-14)

    # This is the function we want to minimize: the node equations
    def F(Vext, voltages):

        # Create the output
        output = np.zeros_like(voltages)

        # First equation
        I1 = solar_cell(0).iv(voltages[0])
        output[0] = (Vext - np.sum(voltages)) / R_series - I1

        # The rest of the equations
        for j in range(1, num_jun):
            I2 = solar_cell(j).iv(voltages[j])

            # Here is where the coupling comes into play
            for i in range(j):
                I2 = I2 - coupling_currents[i, j - 1] * np.exp((q * voltages[i]) / (kb * T))

            output[j] = I1 - I2

            I1 = I2

        return output

    # Now we loop over all the voltages
    # We get the initial condition from the solution without coupling
    guess = V_junction_array[0]
    output_J = np.zeros_like(output_V)

    for i, VV in enumerate(output_V):
        # We particularise the function to minimize for the current output voltage
        def FF(voltages):
            return F(VV, voltages)

        # Solve the problem
        result = root(FF, guess, method='lm', tol=1e-10)
        V_junction_array[i] = result['x']

        # Finally, we calculate the current, based on the current of the first junction
        output_J[i] = solar_cell(0).iv(result['x'][0])

        guess = V_junction_array[i]

    # We calculate the coupling current and add it to the output
    coupled_J = []
    for j in range(1, num_jun):
        coupled = []
        for i in range(j):
             coupled.append((V_junction_array[:, i], coupling_currents[i, j - 1] * np.exp((q * V_junction_array[:, i]) / (kb * T))))

        coupled_J.append(coupled)

    return output_J, V_junction_array, coupled_J
