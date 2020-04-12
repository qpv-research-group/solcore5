import numpy as np
from scipy.interpolate import interp1d

from solcore import asUnit, constants
from solcore.state import State
from .ddModel import driftdiffusion as dd
from .DeviceStructure import LoadAbsorption, CreateDeviceStructure, CalculateAbsorptionProfile

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


# Functions for creating the streucture in the fortran variables and executing the DD solver
def ProcessStructure(device, meshpoints, wavelengths=None):
    """ This function reads a dictionary containing all the device structure, extract the electrical and optical
    properties of the materials, and loads all that information into the Fortran variables. Finally, it initialises the
    device (in Fortran) calculating an initial mesh and all the properties as a function of the position.

    :param device: A dictionary containing the device structure. See PDD.DeviceStructure
    :param wavelengths: (Optional) Wavelengths at which to calculate the optical properties.
    :return: Dictionary containing the device structure properties as a function of the position.
    """
    print('Processing structure...')
    # First, we clean any previous data from the Fortran code
    dd.reset()
    output = {}

    if wavelengths is not None:
        output['Optics'] = {}
        output['Optics']['wavelengths'] = wavelengths

    # We dump the structure information to the Fotran module and initialise the structure
    i = 0
    while i < device['numlayers']:
        layer = device['layers'][i]['properties']
        args_list = [layer['width'],
                     asUnit(layer['band_gap'], 'eV'),
                     asUnit(layer['electron_affinity'], 'eV'),
                     layer['electron_mobility'],
                     layer['hole_mobility'],
                     layer['Nc'],
                     layer['Nv'],
                     layer['electron_minority_lifetime'],
                     layer['hole_minority_lifetime'],
                     layer['relative_permittivity'] * Epsi0,
                     layer['radiative_recombination'],
                     layer['electron_auger_recombination'],
                     layer['hole_auger_recombination'],
                     layer['Na'],
                     layer['Nd']]
        dd.addlayer(args=args_list)
        i = i + device['layers'][i]['numlayers']

    # We set the surface recombination velocities. This needs to be improved at some point
    # to consider other boundary conditions
    dd.frontboundary("ohmic", device['sn'], device['sp'], 0)
    dd.backboundary("ohmic", device['sn'], device['sp'], 0)

    dd.initdevice(meshpoints)
    print('...done!\n')

    output['Properties'] = DumpInputProperties()
    return output


def equilibrium_pdd(junction, options):
    """ Solves the PDD equations under equilibrium: in the dark with no external current and zero applied voltage.

    :param junction: A junction object
    :param options: Options to be passed to the solver
    :return: None
    """
    T = options.T
    wl = options.wavelength
    output_info = options.output_equilibrium

    SetConvergenceParameters(options)
    SetMeshParameters(options)
    SetRecombinationParameters(options)

    device = CreateDeviceStructure('Junction', T=T, layers=junction)

    print('Solving equilibrium...')
    output = ProcessStructure(device, options.meshpoints, wavelengths=wl)
    dd.gen = 0

    dd.equilibrium(output_info)
    print('...done!\n')

    output['Bandstructure'] = DumpBandStructure()
    junction.equilibrium_data = State(**output)


def short_circuit_pdd(junction, options):
    """ Solves the devices electronic properties at short circuit. Internally, it calls Equilibrium.

    :param junction: A junction object
    :param options: Options to be passed to the solver
    :return: None
    """

    # We run equilibrium
    equilibrium_pdd(junction, options)

    wl = options.wavelength
    wl, ph = options.light_source.spectrum(x=wl, output_units='photon_flux_per_m')

    z = junction.equilibrium_data['Bandstructure']['x']
    absorption = junction.absorbed if hasattr(junction, 'absorbed') else lambda x: 0
    abs_profile = CalculateAbsorptionProfile(z, wl, absorption)

    dd.set_generation(abs_profile)
    dd.gen = 1

    dd.illumination(ph)

    dd.lightsc(options.output_sc, 1)
    print('...done!\n')

    output = {'Bandstructure': DumpBandStructure()}
    junction.short_circuit_data = State(**output)


def iv_pdd(junction, options):
    """ Calculates the IV curve of the device between 0 V and a given voltage. Depending on the options, the IV will be calculated in the dark (calling the equilibrium_pdd function) or under illumination (calling the short_circuit_pdd function). If the voltage range has possitive and negative values, the problem is solved twice: from 0 V to the maximu positive and from 0 V to the maximum negative, concatenating the results afterwards.

    :param junction: A junction object
    :param options: Options to be passed to the solver
    :return: None
    """
    print('Solving IV...')
    light = options.light_iv
    output_info = options.output_iv
    T = options.T

    junction.voltage = options.internal_voltages

    Eg = 10
    for layer in junction:
        Eg = min([Eg, layer.material.band_gap / q])

    s = True if junction[0].material.Na >= junction[0].material.Nd else False

    if s:
        vmax = min(Eg + 3 * kb * T / q, max(junction.voltage))
        vmin = min(junction.voltage)
    else:
        vmax = max(junction.voltage)
        vmin = max(-Eg - 3 * kb * T / q, min(junction.voltage))

    vstep = junction.voltage[1] - junction.voltage[0]

    # POSITIVE RANGE
    output_pos = None
    if vmax > 0:
        if light:
            short_circuit_pdd(junction, options)
        else:
            equilibrium_pdd(junction, options)

        dd.runiv(vmax, vstep, output_info, 0)
        output_pos = {'Bandstructure': DumpBandStructure(), 'IV': DumpIV()}

    # NEGATIVE RANGE
    output_neg = None
    if vmin < 0:
        if light:
            short_circuit_pdd(junction, options)
        else:
            equilibrium_pdd(junction, options)

        dd.runiv(vmin, (-1 * vstep), output_info, 0)
        output_neg = {'Bandstructure': DumpBandStructure(), 'IV': DumpIV()}

    print('...done!\n')

    # Now we need to put together the data for the possitive and negative regions.
    junction.pdd_data = State({'possitive_V': output_pos, 'negative_V': output_neg})

    if output_pos is None:
        V = output_neg['IV']['V']
        J = output_neg['IV']['J']
        Jrad = output_neg['IV']['Jrad']
        Jsrh = output_neg['IV']['Jsrh']
        Jaug = output_neg['IV']['Jaug']
        Jsur = output_neg['IV']['Jsur']
    elif output_neg is None:
        V = output_pos['IV']['V']
        J = output_pos['IV']['J']
        Jrad = output_pos['IV']['Jrad']
        Jsrh = output_pos['IV']['Jsrh']
        Jaug = output_pos['IV']['Jaug']
        Jsur = output_pos['IV']['Jsur']
    else:
        V = np.concatenate((output_neg['IV']['V'][:0:-1], output_pos['IV']['V']))
        J = np.concatenate((output_neg['IV']['J'][:0:-1], output_pos['IV']['J']))
        Jrad = np.concatenate((output_neg['IV']['Jrad'][:0:-1], output_pos['IV']['Jrad']))
        Jsrh = np.concatenate((output_neg['IV']['Jsrh'][:0:-1], output_pos['IV']['Jsrh']))
        Jaug = np.concatenate((output_neg['IV']['Jaug'][:0:-1], output_pos['IV']['Jaug']))
        Jsur = np.concatenate((output_neg['IV']['Jsur'][:0:-1], output_pos['IV']['Jsur']))

    # Finally, we calculate the currents at the desired voltages
    R_shunt = min(junction.R_shunt, 1e14) if hasattr(junction, 'R_shunt') else 1e14

    junction.current = np.interp(junction.voltage, V, J) + junction.voltage / R_shunt
    Jrad = np.interp(junction.voltage, V, Jrad)
    Jsrh = np.interp(junction.voltage, V, Jsrh)
    Jaug = np.interp(junction.voltage, V, Jaug)
    Jsur = np.interp(junction.voltage, V, Jsur)

    junction.iv = interp1d(junction.voltage, junction.current, kind='linear', bounds_error=False, assume_sorted=True,
                           fill_value=(junction.current[0], junction.current[-1]))
    junction.recombination_currents = State({"Jrad": Jrad, "Jsrh": Jsrh, "Jaug": Jaug, "Jsur": Jsur})


def qe_pdd(junction, options):
    """ Calculates the quantum efficiency of the device at short circuit. Internally it calls ShortCircuit

    :param junction: A junction object
    :param options: Options to be passed to the solver
    :return: None
    """
    print('Solving quantum efficiency...')
    output_info = options.output_qe

    short_circuit_pdd(junction, options)

    dd.runiqe(output_info)
    print('...done!\n')

    output = {}
    output['QE'] = DumpQE()
    output['QE']['WL'] = options.wavelength

    z = junction.short_circuit_data['Bandstructure']['x']
    absorbed_per_wl = np.trapz(junction.absorbed(z), z, axis=0)

    # This is redundant but useful to keep the same format than the other solvers
    junction.qe = State(**output['QE'])
    junction.qe['IQE'] = junction.qe['EQE']/np.maximum(absorbed_per_wl, 0.00001)

    # The EQE is actually the IQE inside the fortran solver due to an error in the naming --> to be changed
    junction.eqe = interp1d(options.wavelength, output['QE']['EQE'], kind='linear', bounds_error=False,
                            assume_sorted=True, fill_value=(output['QE']['EQE'][0], output['QE']['EQE'][-1]))


# ----
# Functions for dumping data from the fortran variables
def DumpInputProperties():
    output = {}
    output['x'] = dd.get('x')[0:dd.m + 1]
    output['Xi'] = dd.get('xi')[0:dd.m + 1]
    output['Eg'] = dd.get('eg')[0:dd.m + 1]
    output['Nc'] = dd.get('nc')[0:dd.m + 1]
    output['Nv'] = dd.get('nv')[0:dd.m + 1]
    output['Nd'] = dd.get('nd')[0:dd.m + 1]
    output['Na'] = dd.get('na')[0:dd.m + 1]

    return output


def DumpBandStructure():
    output = {}
    output['x'] = dd.get('x')[0:dd.m + 1]
    output['n'] = dd.get('n')[0:dd.m + 1]
    output['p'] = dd.get('p')[0:dd.m + 1]
    output['ni'] = dd.get('ni')[0:dd.m + 1]
    output['Rho'] = dd.get('rho')[0:dd.m + 1]
    output['Efe'] = dd.get('efe')[0:dd.m + 1]
    output['Efh'] = dd.get('efh')[0:dd.m + 1]
    output['potential'] = dd.get('psi')[0:dd.m + 1]
    output['Ec'] = dd.get('ec')[0:dd.m + 1]
    output['Ev'] = dd.get('ev')[0:dd.m + 1]
    output['GR'] = dd.get('gr')[0:dd.m + 1]
    output['G'] = dd.get('g')[0:dd.m + 1]
    output['Rrad'] = dd.get('rrad')[0:dd.m + 1]
    output['Rsrh'] = dd.get('rsrh')[0:dd.m + 1]
    output['Raug'] = dd.get('raug')[0:dd.m + 1]

    return output


def DumpIV(IV_info=False):
    # Depending of having PN or NP the calculation of the MPP is a bit different. We move everithing to the 1st quadrant and then send it back to normal
    Nd = dd.get('nd')[0:dd.m + 1][0]
    Na = dd.get('na')[0:dd.m + 1][0]
    s = Nd > Na

    output = {}
    output['V'] = (-1) ** s * dd.get('volt')[1:dd.nvolt + 1]
    output['J'] = dd.get('jtot')[1:dd.nvolt + 1]
    output['Jrad'] = dd.get('jrad')[1:dd.nvolt + 1]
    output['Jsrh'] = dd.get('jsrh')[1:dd.nvolt + 1]
    output['Jaug'] = dd.get('jaug')[1:dd.nvolt + 1]
    output['Jsur'] = dd.get('jsur')[1:dd.nvolt + 1]

    if IV_info:
        # We calculate the solar cell parameters
        output['Jsc'] = -np.interp(0, output['V'], output['J'])  # dd.get('isc')[0]
        output['Voc'] = np.interp(0, output['J'][output['V'] > 0], output['V'][output['V'] > 0])  # dd.get('voc')[0]
        print(output['Voc'])

        maxPP = np.argmin(output['V'] * output['J'])
        Vmax = output['V'][maxPP - 3:maxPP + 3]
        Imax = output['J'][maxPP - 3:maxPP + 3]
        Pmax = Vmax * Imax

        poly = np.polyfit(Vmax, Pmax, 2)
        output['Vmpp'] = -poly[1] / (2 * poly[0])
        output['Jmpp'] = -output['J'][np.argmin(np.abs(output['V'] - output['Vmpp']))]
        output['FF'] = output['Vmpp'] * output['Jmpp'] / (output['Voc'] * output['Jsc'])

        print("Jsc  = %5.3f mA/cm2" % (output['Jsc'] / 10))
        print("Voc  = %4.3f  V" % (output['Voc']))
        print("FF   = %3.3f " % (output['FF']))
        print("Jmpp = %5.3f mA/cm2" % (output['Jmpp'] / 10))
        print("Vmpp = %4.3f  V" % (output['Vmpp']))
        print("Power= %5.3f mW/cm2" % (output['Jmpp'] * output['Vmpp'] / 10))

    # If NP, V and J should be in the 3rd quadrant
    output['V'] = (-1) ** s * output['V']
    output['J'] = (-1) ** s * output['J']
    output['Jrad'] = (-1) ** s * output['Jrad']
    output['Jsrh'] = (-1) ** s * output['Jsrh']
    output['Jaug'] = (-1) ** s * output['Jaug']
    output['Jsur'] = (-1) ** s * output['Jsur']

    return output


def DumpQE():
    output = {}
    numwl = dd.numwl + 1
    output['EQE'] = dd.get('iqe')[0:numwl]
    output['EQEsrh'] = dd.get('iqesrh')[0:numwl]
    output['EQErad'] = dd.get('iqerad')[0:numwl]
    output['EQEaug'] = dd.get('iqeaug')[0:numwl]
    output['EQEsurf'] = dd.get('iqesurf')[0:numwl]
    output['EQEsurb'] = dd.get('iqesurb')[0:numwl]

    return output


# ----
# Functions for setting the parameters controling the recombination, meshing and the numerial algorithm
def SetMeshParameters(options):
    dd.set('coarse', options.coarse)
    dd.set('fine', options.fine)
    dd.set('ultrafine', options.ultrafine)
    dd.set('growth', options.growth_rate)


def SetRecombinationParameters(options):
    dd.srh = options.srh
    dd.rad = options.rad
    dd.aug = options.aug
    dd.sur = options.sur
    dd.gen = options.gen


def SetConvergenceParameters(options):
    dd.nitermax = options.nitermax
    dd.set('clamp', options.clamp)
    dd.set('atol', options.ATol)
    dd.set('rtol', options.RTol)
