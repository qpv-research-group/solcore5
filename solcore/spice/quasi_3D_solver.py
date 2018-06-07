import numpy as np
from solcore.solar_cell_solver import solar_cell_solver
from solcore.spice import solve_circuit


def solve_quasi_3D(solar_cell, injection, contacts, options=None, Lx=10e-6, Ly=10e-6, h=2e-6, R_back=1e-16,
                   R_contact=1e-16, R_line=1e-16, bias_start=0, bias_end=1.8, bias_step=0.01):
    """ Entry function for the quasi-3D solver

    :param solar_cell: A solar cell object
    :param injection: 2D array indicating the (optical) injection mask
    :param contacts: 2D array indicating the electrical contacts
    :param options: Options for the 1D solar cell solver
    :param Lx: Pixel size in the X direction
    :param Ly: Pixel size in the Y direction
    :param h: Height of the metal fingers
    :param R_back: Resistance back contact
    :param R_contact: Contact resistance
    :param R_line: Resistivity metal fingers
    :param bias_start: Initial voltage (V)
    :param bias_end: Final voltage (V)
    :param bias_step: Voltage step (V)
    :return: A tuple with:

        - V [steps + 1] : 1D Array with the external voltages
        - I [steps + 1] : 1D Array with the current at all external V
        - Vall [xnodes, ynodes, 2 * junctions, steps + 1] : 4D Array with the voltages in all nodes, at all external V
        - Vmet [xnodes, ynodes, steps + 1] : 3D Array with the voltages in the metal nodes, at all external V
    """
    # We first start by the solar cell as if it were a normal, isolated cell
    print("Solving 1D Solar Cell...")
    solar_cell_solver(solar_cell, 'iv', user_options=options)
    print("... Done!\n")

    # We don't care about this IV curve, in principle, but we care about some of the parameters calculated, like jsc,
    # j01 or j02 if calculated from detailed balance. We extract those parameters from the cell
    totaljuncs = solar_cell.junctions

    Isc_array = np.zeros(totaljuncs)
    I01_array = np.zeros(totaljuncs)
    n1_array = np.zeros(totaljuncs)
    I02_array = np.zeros(totaljuncs)
    n2_array = np.zeros(totaljuncs)
    Eg_array = np.zeros(totaljuncs)
    rsh_array = np.zeros(totaljuncs)
    rshTop_array = np.zeros(totaljuncs)
    rshBot_array = np.zeros(totaljuncs)
    rseries_array = np.ones(totaljuncs) * 1e-16

    for i in range(totaljuncs):

        n1_array[i] = solar_cell(i).n1 if hasattr(solar_cell(i), 'n1') else 1
        n2_array[i] = solar_cell(i).n2 if hasattr(solar_cell(i), 'n2') else 2
        rsh_array[i] = min(solar_cell(i).R_shunt, 1e16) if hasattr(solar_cell(i), 'R_shunt') else 1e16
        rshTop_array[i] = max(solar_cell(i).R_sheet_top, 1e-16) if hasattr(solar_cell(i), 'R_sheet_top') else 1e-16
        rshBot_array[i] = max(solar_cell(i).R_sheet_bot, 1e-16) if hasattr(solar_cell(i), 'R_sheet_bot') else 1e-16

        try:
            Isc_array[i] = solar_cell(i).jsc
            I01_array[i] = solar_cell(i).j01
            I02_array[i] = solar_cell(i).j02
            Eg_array[i] = solar_cell(i).Eg
        except AttributeError as err:
            raise AttributeError('ERROR in quasi-3D solver: Junction is missing one essential argument. {}'.format(err))

    j = 0
    for i in solar_cell.tunnel_indices:
        rseries_array[j] = max(solar_cell[i].R_series, 1e-16) if hasattr(solar_cell[i], 'R_series') else 1e-16
        j += 1

    rseries_array[-1] = max(R_back, 1e-16)

    print("Solving quasi-3D Solar Cell...")
    V, I, Vall, Vmet = solve_circuit_quasi3D(bias_start, bias_end, bias_step, Isc_array, I01_array, I02_array, n1_array,
                                             n2_array, Eg_array, rsh_array, rseries_array, injection, contacts,
                                             rshTop_array, rshBot_array, R_line / h, R_contact, Lx, Ly)
    print("... Done!!")
    return V, I, Vall, Vmet


def create_node(type, idx, idy, Lx, Ly, Isc, topLCL, botLCL, rshunt, rseries, xMetalTop, yMetalTop, contact):
    """ Creates a node of the solar cell, meaning all the circuit elements at an XY location in the plane. This includes all the diodes, resistances and current sources for all the junctions at that location.

    :param type: The type of the node, 'Normal', 'Finger' or 'Bus'
    :param idx: Index with the location in the X direction
    :param idy: Index with the location in the Y direction
    :param Lx: Pixel size in the X direction
    :param Ly: Pixel size in the Y direction
    :param Isc: Array of Isc for each of the junctions
    :param topLCL: Array of resistances of the top lateral conductive layer
    :param botLCL: Array of resistances of the bottom lateral conductive layers
    :param rshunt: Array of Rshunt for each of the junctions
    :param rseries: Array of Rseries for each of the junctions
    :param xMetalTop: Resistance of the metal in the X direction
    :param yMetalTop: Resistance of the metal in the Y direction
    :param contact: Contact resistance
    :return: The node define in SPICE file as a string.
    """
    node = ''
    for j in range(len(Isc)):
        loc = str(j) + "_" + str(idx).zfill(3) + "_" + str(idy).zfill(3)
        locXR = str(j) + "_" + str(idx + 1).zfill(3) + "_" + str(idy).zfill(3)
        locYR = str(j) + "_" + str(idx).zfill(3) + "_" + str(idy + 1).zfill(3)

        if j + 1 == len(Isc):
            locLow = 0
        else:
            locLow = "t_" + str(j + 1) + "_" + str(idx).zfill(3) + "_" + str(idy).zfill(3)

        s = Ly / Lx

        # We add the diodes
        diode1 = "d1_{0} t_{0} b_{0} diode1_{1}\n".format(loc, j)
        diode2 = "d2_{0} t_{0} b_{0} diode2_{1}\n".format(loc, j)

        # Now the shunt resistance
        rshuntJ = "Rshunt_{0} t_{0} b_{0} {1}\n".format(loc, rshunt[j])

        # And add the source
        source = 'i{0} b_{0} t_{0} {1}\n'.format(loc, Isc[j])

        # Now we add the sheet resistances
        rbotLCLX = "RbX{0}to{1} b_{0} b_{1} {2}\n".format(loc, locXR, botLCL[j] / s)
        rbotLCLY = "RbY{0}to{1} b_{0} b_{1} {2}\n".format(loc, locYR, botLCL[j] * s)
        rtopLCLX = "RtX{0}to{1} t_{0} t_{1} {2}\n".format(loc, locXR, topLCL[j] / s)
        rtopLCLY = "RtY{0}to{1} t_{0} t_{1} {2}\n".format(loc, locYR, topLCL[j] * s)

        # Now the series resistance with the back of the junction
        rseriesJ = "Rseries{0}to{1} b_{0} {1} {2}\n".format(loc, locLow, rseries[j])

        if j == 0 and type == 'Finger':
            rcontact = "Rcontact{0} t_{0} m_{0} {1}\n".format(loc, contact)
            rmetalX = "RbusX{0}to{1} m_{0} m_{1} {2}\n".format(loc, locXR, xMetalTop / s)
            rmetalY = "RbusY{0}to{1} m_{0} m_{1} {2}\n".format(loc, locYR, yMetalTop * s)
            rext = ""

        elif j == 0 and type == 'Bus':
            rcontact = "Rcontact{0} t_{0} m_{0} {1}\n".format(loc, contact)
            rmetalX = "RbusX{0}to{1} m_{0} m_{1} {2}\n".format(loc, locXR, 1e-16)
            rmetalY = "RbusY{0}to{1} m_{0} m_{1} {2}\n".format(loc, locYR, 1e-16)

            # This is the connection to the external voltage
            rext = "Rext{0} in m_{0} {1}\n".format(loc, 1e-16)

        else:
            rcontact = ""
            rmetalX = ""
            rmetalY = ""
            rext = ""

        # Finally, we create the output statement for this node
        if j == 0 and type in ['Finger', 'Bus']:
            output = ".PRINT DC v(t_{0}) v(b_{0}) v(m_{0})\n\n".format(loc)
        else:
            output = ".PRINT DC v(t_{0}) v(b_{0})\n\n".format(loc)

        # and put all the instructtions together
        node = node + diode1 + diode2 + source + rshuntJ + rtopLCLX + rtopLCLY + rbotLCLX + rbotLCLY + rseriesJ + \
               rcontact + rmetalX + rmetalY + rext + output

    return node


def create_header(I01, I02, n1, n2, Eg, T=20):
    """ Creates the header of the SPICE file, where the diode models, the temperature and the independent voltage source are defined.

    :param I01: Array of I01 for each of the junctions
    :param I02: Array of I02 for each of the junctions
    :param n1: Array of n1 for each of the junctions
    :param n2: Array of n2 for each of the junctions
    :param Eg: Array of Eg for each of the junctions
    :param T: Temperature of the device
    :return: The header of the SPICE file as a string.
    """
    title = "*** A SPICE simulation with python\n\n"

    diodes = ""
    for j in range(len(I01)):
        modelDiode1 = ".model diode1_{0} d(is={1},n={2},eg={3})\n".format(j, I01[j], n1[j], Eg[j])
        modelDiode2 = ".model diode2_{0} d(is={1},n={2},eg={3})\n".format(j, I02[j], n2[j], Eg[j])

        diodes = diodes + modelDiode1 + modelDiode2

    options = ".OPTIONS TNOM=20 TEMP={0}\n\n".format(T)
    independent_source = """ 
    vdep in 0 DC 0
    """

    SPICEheader = title + diodes + options + independent_source

    return SPICEheader


def solve_circuit_quasi3D(vini, vfin, step, Isc, I01, I02, n1, n2, Eg, Rshunt, Rseries, injection, contacts, RsTop,
                          RsBot, Rline, Rcontact, Lx, Ly):
    """ This is the function that actually dumps all the information to the Spice engine, runs the calculation, and retrieves the datafrom the calculator.

    :param vini: Initial voltage (V)
    :param vfin: Final voltage (V)
    :param step: Voltage step (V)
    :param Isc: Array of Isc for each of the junctions
    :param I01: Array of I01 for each of the junctions
    :param I02: Array of I02 for each of the junctions
    :param n1: Array of n1 for each of the junctions
    :param n2: Array of n2 for each of the junctions
    :param Eg: Array of Eg for each of the junctions
    :param Rshunt: Array of Rshunt for each of the junctions
    :param Rseries: Array of Rseries for each of the junctions
    :param injection: 2D array indicating the (optical) injection mask
    :param contacts: 2D array indicating the electrical contacts
    :param RsTop: Array of sheet resistance on the top for each of the junctions
    :param RsBot: Array of sheet resistance on the bottom for each of the junctions
    :param Rline: Resistance of the metal fingers
    :param Rcontact: Contact resistance
    :param Lx: Pixel size in the X direction
    :param Ly: Pixel size in the Y direction
    :return: A tuple with:

        - V [steps + 1] : 1D Array with the external voltages
        - I [steps + 1] : 1D Array with the current at all external V
        - Vall [xnodes, ynodes, 2 * junctions, steps + 1] : 4D Array with the voltages in all nodes, at all external V
        - Vmet [xnodes, ynodes, steps + 1] : 3D Array with the voltages in the metal nodes, at all external V
    """

    # Scaling factor to bring the magnitudes to a regime where the solver is comfortable
    gn = np.sqrt(1.0 / I01[0])

    areaPerPixel = Lx * Ly

    isc = Isc * areaPerPixel * gn
    i01 = I01 * areaPerPixel * gn
    i02 = I02 * areaPerPixel * gn
    rsTop = RsTop / gn
    rsBot = RsBot / gn
    rseries = Rseries / areaPerPixel / gn
    rshunt = Rshunt / areaPerPixel / gn
    contact = Rcontact / areaPerPixel / gn
    metal = Rline / gn

    illumination = injection / 255
    pads = np.where(contacts > 200, 1, 0)
    shadow = np.where(contacts > 55, 1, 0)

    xnodes, ynodes = injection.shape
    junctions = len(I01)
    steps = round((vfin - vini) / step)

    SPICEheader = create_header(i01, i02, n1, n2, Eg)
    SPICEfooter = ".end"

    # We prepare the SPICE execution
    SPICEexec = ".PRINT DC i(vdep)\n.DC vdep {0} {1} {2}\n".format(vini, vfin, step)

    # Now we prepare the subcircuits on each node
    SPICEbody = """"""

    x = np.arange(xnodes)
    y = np.arange(ynodes)
    Vall = np.zeros((xnodes, ynodes, 2 * junctions, steps + 1))
    Vmet = np.zeros((xnodes, ynodes, steps + 1))
    I = np.zeros(steps + 1)
    V = np.zeros(steps + 1)

    # We create the solar cell units and the series resistances at each node
    for xx in x:
        for yy in y:

            # This calculates the metal resistance depending on having metal on top or not. It leaves some dangling
            # resistors in the circuit, but it shouldn't be a problem.
            metalR = max(metal * shadow[xx, yy], 1e-16) + 1e-16 * (1 - shadow[xx, yy])

            if not shadow[xx, yy]:
                # we create a normal node
                SPICEbody = SPICEbody + create_node('Normal', xx, yy, Lx, Ly, Isc=illumination[xx, yy] * isc,
                                                    topLCL=rsTop, botLCL=rsBot, rshunt=rshunt, rseries=rseries,
                                                    xMetalTop=metalR, yMetalTop=metalR, contact=contact)
            elif pads[xx, yy]:
                # we create at bus node, with no resistance in the metal and direct electrical injection
                SPICEbody = SPICEbody + create_node('Bus', xx, yy, Lx, Ly, Isc=0 * isc, topLCL=rsTop,
                                                    botLCL=rsBot, rshunt=rshunt, rseries=rseries, xMetalTop=metalR,
                                                    yMetalTop=metalR, contact=contact)
            else:
                # We create a finger node, with resistance in the metal and not direct injection
                SPICEbody = SPICEbody + create_node('Finger', xx, yy, Lx, Ly, Isc=0 * isc, topLCL=rsTop,
                                                    botLCL=rsBot, rshunt=rshunt, rseries=rseries, xMetalTop=metalR,
                                                    yMetalTop=metalR, contact=contact)

    # We combine the different bits to create the SPICE input file
    SPICEcommand = SPICEheader + SPICEbody + SPICEexec + SPICEfooter
    raw_results = solve_circuit(SPICEcommand)

    # The raw results are are a very long chunk of text. We have to clean it and pick just the info we want,
    # that is the voltages at certain nodes.
    lines = raw_results.split("\n")
    for line in lines:
        if len(line) == 0:
            continue

        if line[:5] == 'Index':
            headers = line.split()

            if len(headers) >= 4:
                indices = headers[2][4:13].split('_')
                junc = int(indices[0])
                x = int(indices[1])
                y = int(indices[2])

        if line[0] not in '01234567890.':
            continue

        i, *rest = line.split()
        i = int(i)

        if len(rest) == 3:
            Vall[x, y, 2 * junc + 1, i] = float(rest[2])
            Vall[x, y, 2 * junc, i] = float(rest[1])
        elif len(rest) == 4:
            Vall[x, y, 2 * junc + 1, i] = float(rest[2])
            Vall[x, y, 2 * junc, i] = float(rest[1])
            Vmet[x, y, i] = float(rest[3])
        else:
            V[i] = float(rest[0])
            I[i] = -float(rest[1])

    # Finally, we un-do the scaling
    I = I / gn

    return V, I, Vall, Vmet
