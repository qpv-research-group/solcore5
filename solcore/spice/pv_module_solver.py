import numpy as np
from solcore.solar_cell_solver import solar_cell_solver
from solcore.spice import solve_circuit


def solve_pv_module(solar_cell, options, totalcells=25, bias_start=0, bias_end=75, bias_step=0.1, jscSigma=2e-4,
                    shading=None):
    """ Calculate the IV curve of a PV module made of a certain number solar_cells connected in series in a single string. A certain dispersion to the distribution of photocurrents among cells and shadowing losses can be added to ahoeve more realistic results.

    :param solar_cell: A solar cell object containing all the junctions
    :param options: A state object with all the options needed to solve the solar cell IV
    :param totalcells: Total number of cells in a string
    :param bias_start: Initial voltage of the IV curve of the module
    :param bias_end: Final voltage of the IV curve of the module
    :param bias_step: Step for the bias
    :param jscSigma: Width of the distribution of short circuit currents around the ideal (eg: 0.02 = 2%). Default is set to 0.02% to aid convergence
    :param shading: Array containing the shading loses for each cell. It can be a 2D array (0 = no shading, 1 = full shading)
    :return: A tuple with the array of voltages, currents, a list with the Isc of each junction in each cell and the raw output from SPICE.
    """
    # We first start by the solar cell as if it were a normal, isolated cell
    solar_cell_solver(solar_cell, 'iv', user_options=options)

    # We don't care about this IV curve, in principle, but we care about some of the parameters calculated, like jsc,
    # j01 or j02 if calculated from detailed balance. We extract those parameters from the cell
    totaljuncs = solar_cell.junctions

    area = solar_cell.area if hasattr(solar_cell, 'area') else 1
    Rs = max(solar_cell.R_series / area, 1e-14)
    cell_temp = solar_cell.T - 273

    Isc_array = np.zeros(totaljuncs)
    I01_array = np.zeros(totaljuncs)
    n1_array = np.zeros(totaljuncs)
    I02_array = np.zeros(totaljuncs)
    n2_array = np.zeros(totaljuncs)
    Eg_array = np.zeros(totaljuncs)
    rsh_array = np.zeros(totaljuncs)

    for i in range(totaljuncs):

        try:
            Isc_array[i] = solar_cell(i).jsc * area
            I01_array[i] = solar_cell(i).j01 * area
            n1_array[i] = solar_cell(i).n1 if hasattr(solar_cell(i), 'n1') else 1
            I02_array[i] = solar_cell(i).j02 * area
            n2_array[i] = solar_cell(i).n2 if hasattr(solar_cell(i), 'n2') else 2
            Eg_array[i] = solar_cell(i).Eg
            rsh_array[i] = min(solar_cell(i).R_shunt / area, 1e14) if hasattr(solar_cell(i), 'R_shunt') else 1e14
        except AttributeError as err:
            raise AttributeError('ERROR in PV module: Junction is missing one essential argument. {}'.format(err))

    if shading is not None:
        shade = shading.flatten()
        assert len(shade) == totalcells
    else:
        shade = np.zeros(totalcells)

    # Calculate the probability distribution for the cells
    # Since all cells will not generate the same current, some random scatter can be added to the Isc.
    # note that the same random value is applied to all sub-cells in a particular junction.
    jsc_random_normal = np.random.normal(0.0, jscSigma, size=totalcells)

    # Shift the distribution so it is always <=1 so Isc never exceeds the ideal value,
    # and also take into account shading
    jsc_random_normal = abs(1 - max(jsc_random_normal) * jscSigma + jscSigma * jsc_random_normal) * (1 - shade)

    # We choose to keep track of the sub-cell Isc values so that they can be plotted at the end.
    all_Isc_values = []

    # Now we start to write the spice file
    spice_file = """ * SolCore spice Simulation \n.model bypassdiode d \n"""

    nodeCounter = 0
    junctionCounter = 0

    # We loop over all solar cells and all junctions within each solar cell
    for n in range(1, totalcells + 1):
        # Defines the lower connection for the bypass diode.
        lowerBypassConnection = nodeCounter
        temp_isc = []

        # Loop through junctions writing the spice code for each junction in turn
        for j in range(totaljuncs):
            junctionCounter += 1

            isc = Isc_array[j] * jsc_random_normal[n - 1]

            # Write the spice code for the junction to file
            spiceout = spice_junction(junctionCounter, nodeCounter, isc, I01_array[j], I02_array[j], n1_array[j],
                                      n2_array[j], Eg_array[j], rsh_array[j])

            spice_file = spice_file + spiceout
            nodeCounter = nodeCounter + 1
            temp_isc.append(isc)  # Save the isc value for later plotting

        # Entire cell has been built
        all_Isc_values.append(temp_isc)  # Add isc values to the main list.

        # Add the series resistance
        spiceout = 'r{0} {1} {2} {3}\n'.format(2 * junctionCounter - 1, nodeCounter+1, nodeCounter, Rs)
        nodeCounter = nodeCounter + 1
        spice_file = spice_file + spiceout

        # Connect bypass diode.  Connections are: lowerBypassConnection & nodeCounter
        spiceout = 'd{0} {1} {2} bypassdiode\n'.format(2 * junctionCounter + 1, lowerBypassConnection, nodeCounter)
        spice_file = spice_file + spiceout

        junctionCounter += 1

    spice_file = spice_file + ".OPTIONS TNOM=25 TEMP=" + str(cell_temp) + "\n"
    spice_file = spice_file + "vbias " + str(nodeCounter) + " 0\n"
    spice_file = spice_file + ".dc vbias " + str(bias_start) + " " + str(bias_end) + " " + str(bias_step) + "\n"
    spice_file = spice_file + ".plot dc i(vbias)\n"
    spice_file = spice_file + ".end \n\n"

    # And the definition of the spice file is finished. Now it's time to solve the problem
    raw_data = solve_circuit(spice_file)

    voltage = []
    current = []

    raw_data_in_lines = raw_data.split('\n')

    for line in raw_data_in_lines:
        if len(line) > 2 and line[1] in '0123456789':
            splitstring = line.split()
            voltage.append(float(splitstring[0]))
            current.append(float(splitstring[1]))

    return voltage, current, all_Isc_values, raw_data


def spice_junction(jc, nc, isc, j01, j02, n1, n2, Eg, rsh):
    """ Creates the string representation in SPICE of the junciton defined by the input values.

    :param jc: junction counter
    :param nc: node counter
    :param isc: Short circuit current
    :param j01: Reverse saturation current corresponding to ideality factor n1
    :param j02: Reverse saturation current corresponding to ideality factor n2
    :param n1: Ideality factor, typically = 1
    :param n2: Ideality factor, typically = 2
    :param Eg: Bandgap of the material the junction is made off
    :param rsh: Shunt resistance in the junction
    :return: A string representation of the junciton in SPICE
    """

    isource = 'i{0} {1} {2} dc {3}\n'.format(jc, nc, nc + 1, isc)
    d1 = 'd{0} {1} {2} diode{3} OFF\n'.format(2 * jc - 1, nc + 1, nc, 2 * jc - 1)
    d1deff = '.model diode{0} d(is={1},n={2},eg={3})\n'.format(2 * jc - 1, j01, n1, Eg)
    d2 = 'd{0} {1} {2} diode{3} OFF\n'.format(2 * jc, nc + 1, nc, 2 * jc)
    d2deff = '.model diode{0} d(is={1},n={2},eg={3})\n'.format(2 * jc, j02, n2, Eg)
    rshunt = 'r{0} {1} {2} {3}\n'.format(2 * jc, nc + 1, nc, rsh)

    junction = isource + d1 + d1deff + d2 + d2deff + rshunt


    return junction


if __name__ == '__main__':
    out = spice_junction(1, 1, 300, 1e-25, 1e-16, 1, 2, 1.4, 0.0001, 1e10)
    print(out)
