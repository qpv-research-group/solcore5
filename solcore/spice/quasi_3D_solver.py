import numpy as np
from solcore.solar_cell_solver import solar_cell_solver
from solcore.spice import solve_circuit


def solve_quasi_3D(solar_cell, options, injection, contacts, Lx=10e-6, Ly=10e-6, h=2e-6, R_back=1e-16,
                   R_contact=1e-16, R_line=1e-16,
                   bias_start=0, bias_end=1.8, bias_step=0.01):
    """ Calculate the IV curve of a PV module made of a certain number solar_cells connected in series in a single string. A certain dispersion to the distribution of photocurrents among cells and shadowing losses can be added to ahoeve more realistic results.

    :param solar_cell: A solar cell object containing all the junctions
    :param options: A state object with all the options needed to solve the solar cell IV
    :param bias_start: Initial voltage of the IV curve of the module
    :param bias_end: Final voltage of the IV curve of the module
    :param bias_step: Step for the bias
    :return: A tuple with the array of voltages, currents, a list with the Isc of each junction in each cell and the raw output from SPICE.
    """
    # We first start by the solar cell as if it were a normal, isolated cell
    print("Solving 1D Solar Cell...")
    solar_cell_solver(solar_cell, 'iv', user_options=options)
    print("... Done!\n")

    # We don't care about this IV curve, in principle, but we care about some of the parameters calculated, like jsc,
    # j01 or j02 if calculated from detailed balance. We extract those parameters from the cell
    totaljuncs = solar_cell.junctions

    cell_temp = solar_cell.T - 273

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
            raise AttributeError('ERROR in PV module: Junction is missing one essential argument. {}'.format(err))

    j = 0
    for i in solar_cell.tunel_indices:
        rseries_array[j] = max(solar_cell[i].R_series, 1e-16) if hasattr(solar_cell[i], 'R_series') else 1e-16
        j += 1

    rseries_array[-1] = max(R_back, 1e-16)

    print("Solving quasi-3D Solar Cell...")
    V, I, Vall, Vmet = solve_circuit_quasi3D(bias_start, bias_end, bias_step, Isc_array, I01_array, I02_array, n1_array,
                                             n2_array, Eg_array, rsh_array, rseries_array, injection, contacts,
                                             rshTop_array, rshBot_array, R_line / h, R_contact, Lx, Ly)
    print("... Done!!")
    return V, I, Vall, Vmet


def create_node(type, idx, idy, Lx, Ly, Isc, topLCL, botLCL, rshunt, rseries, xMetalTop, yMetalTop,
                contact):
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


def solve_circuit_quasi3D(vini, vfin, step, Isc, I01, I02, n1, n2, Eg, Rshunt, Rseries, injection, contacts,
                          RsTop, RsBot, Rline, Rcontact, Lx, Ly):

    gn = np.sqrt(1.0 / I01[0])  # Scaling factor to bring the magnitudes to a regime where the solver is comfortable

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

    illumination = injection/255
    pads = np.where(contacts == 255, 1, 0)
    shadow = np.where(contacts > 0, 1, 0)

    xnodes = injection.shape[0]
    ynodes = injection.shape[1]
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

            metalX = 1e12
            metalY = 1e12
            if xx < xnodes - 1:
                isMetalX = shadow[xx, yy] and shadow[xx - 1, yy]
                metalX = max(metal * isMetalX, 1e-12) + metalX * (1 - isMetalX)
            if yy < ynodes - 1:
                isMetalY = shadow[xx, yy] and shadow[xx, yy - 1]
                metalY = max(metal * isMetalY, 1e-12) + metalY * (1 - isMetalY)

            if not shadow[xx, yy]:
                # we create a normal node
                SPICEbody = SPICEbody + create_node('Normal', xx, yy, Lx, Ly, Isc=illumination[xx, yy] * isc,
                                                    topLCL=rsTop, botLCL=rsBot, rshunt=rshunt, rseries=rseries,
                                                    xMetalTop=metalX, yMetalTop=metalY, contact=contact)
            elif pads[xx, yy]:
                # we create at bus node, with no resistance in the metal and direct electrical injection
                SPICEbody = SPICEbody + create_node('Bus', xx, yy, Lx, Ly, Isc=0 * isc, topLCL=rsTop,
                                                    botLCL=rsBot, rshunt=rshunt, rseries=rseries, xMetalTop=metalX,
                                                    yMetalTop=metalY, contact=contact)
            else:
                # We create a finger node, with resistance in the metal and not direct injection
                SPICEbody = SPICEbody + create_node('Finger', xx, yy, Lx, Ly, Isc=0 * isc, topLCL=rsTop,
                                                    botLCL=rsBot, rshunt=rshunt, rseries=rseries, xMetalTop=metalX,
                                                    yMetalTop=metalY, contact=contact)


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

    I = I / gn

    return V, I, Vall, Vmet
