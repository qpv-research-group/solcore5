import numpy as np
from scipy.interpolate import interp1d

from solcore import siUnits, si, asUnit, constants, material
from solcore.structure import Layer, Junction
from solcore.state import State
from .ddModel import driftdiffusion as dd
from solcore.light_source import LightSource
from .DeviceStructure import LoadAbsorption, calculate_reflection, calculate_optics, CreateDeviceStructure

Epsi0 = constants.vacuum_permittivity
q = constants.q
pi = constants.pi
h = constants.h
kb = constants.kb
m0 = constants.electron_mass

log = dd.log_file

mesh_control = {'meshpoints': -400, 'growth_rate': 0.7, 'coarse': 20e-9, 'fine': 1e-9, 'ultrafine': 0.2e-9}
convergence_control = {'clamp': 20, 'nitermax': 100, 'ATol': 1e14, 'RTol': 1e-6}
recombination_control = {'srh': 1, 'rad': 1, 'aug': 0, 'sur': 1, 'gen': 0}

default_wavelengths = np.linspace(300e-9, 1100e-9, 200)
default_photon_flux = LightSource(source_type='standard', version='AM1.5g', x=default_wavelengths,
                                  output_units='photon_flux_per_m').spectrum()[1]


def iv_pdd(junction, options):
    T = options.T
    light = options.light_iv
    mpp = options.mpp

    junc = CreateDeviceStructure('Junction', T=T, layers=junction)

    # We set the surface recombination of the first and last layer to those indicated in the Junction
    if hasattr(junction, 'sn'):
        junc['layers'][0]['properties']['sn'] = junction.sn
        junc['layers'][-1]['properties']['sn'] = junction.sn

    if hasattr(junction, 'sp'):
        junc['layers'][0]['properties']['sp'] = junction.sp
        junc['layers'][-1]['properties']['sp'] = junction.sp

    # We define a sensible voltage range in which to calculate the results to avoid going above the bandgap
    Eg = 10
    for lay in junc['layers']:
        Eg = min(Eg, lay['properties']['band_gap'] / q)

    s = junc['layers'][0]['properties']['Nd'] > junc['layers'][0]['properties']['Na']
    junction.voltage = options.internal_voltages
    if not s:
        volt = np.where(junction.voltage < Eg - 3 * kb * T / q, junction.voltage, Eg - 3 * kb * T / q)
    else:
        volt = np.where(junction.voltage > - Eg + 3 * kb * T / q, junction.voltage, -Eg + 3 * kb * T / q)

    vmin = min(volt)
    vmax = max(volt)
    vstep = junction.voltage[1] - junction.voltage[0]

    # Now it is time to perform the calculation, separating the positive and negative ranges
    if not light:
        output_pos = IV(junc, vfin=vmax, vstep=vstep, IV_info=False) if vmax > 0 else None
        output_neg = IV(junc, vfin=vmin, vstep=(-1 * vstep), IV_info=False) if vmin < 0 else None

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


def solve_pdd(solar_cell, task, wavelength=None, incident_light=None, output_info=1, rs=0, vfin=1, vstep=0.01,
              light=False, IV_info=True, escape=True):
    """ Solves the chosen task of the Poisson-drift-difussion solver taking as input a standard SolarCell object with
    any number of junctions. The only necessary inputs are the solar_cell and the task to be done. All other parameters
    have default values.

    :param solar_cell:
    :param task:
    :param wavelength:
    :param incident_light:
    :param output_info:
    :param rs:
    :param vfin:
    :param vstep:
    :param light:
    :param IV_info:
    :param escape:
    :return:
    """

    # Get the energy range and incident spectrum. If they are not inputs, we create some sensible values.
    if wavelength is None:
        if incident_light is not None:
            wavelength = incident_light[0]
            bs = np.copy(incident_light[1])
        else:
            wavelength = default_wavelengths
            bs = default_photon_flux
            print('Using default light source, the AM1.5d solar spectrum between 300 and 1100 nm.')
    else:
        if incident_light is not None:
            bs = np.interp(wavelength, incident_light[0], incident_light[1])
        else:
            bs = np.interp(wavelength, default_wavelengths, default_photon_flux)
            print('Using default light source, the AM1.5d solar spectrum between 300 and 1100 nm.')

    bs_initial = np.copy(bs)

    # We include the shadowing losses
    if hasattr(solar_cell, 'shading'):
        bs *= (1 - solar_cell.shading)

    # And the reflexion losses
    if hasattr(solar_cell, 'reflectivity') and solar_cell.reflectivity is not None:
        ref = solar_cell.reflectivity(wavelength)
        bs *= (1 - ref)
        reflected = ref * bs_initial
    else:
        reflected = np.zeros_like(bs)

    R = reflected / bs_initial

    # And now we perform the calculation, each junction at a time
    passive_loss = np.ones_like(bs)

    for i, layer_object in enumerate(solar_cell):

        # Attenuation due to absorption in the AR coatings or any layer in the front that is not part of the
        # junction
        if type(layer_object) is Layer:
            bs = bs * np.exp(-layer_object.material.alpha(wavelength) * layer_object.width)
            passive_loss *= np.exp(-layer_object.material.alpha(wavelength) * layer_object.width)

        # For each junction, we calculate the chosen task
        elif type(layer_object) is Junction:

            junction = CreateDeviceStructure('Junction', T=solar_cell.T, layers=layer_object,
                                             substrate=solar_cell.substrate, reflection=solar_cell.reflectivity)

            # We set the surface recombination of the first and last layer to those indicated in the Junction
            if hasattr(layer_object, 'sn'):
                junction['layers'][0]['properties']['sn'] = layer_object.sn
                junction['layers'][-1]['properties']['sn'] = layer_object.sn

            if hasattr(layer_object, 'sp'):
                junction['layers'][0]['properties']['sp'] = layer_object.sp
                junction['layers'][-1]['properties']['sp'] = layer_object.sp

            if task == 'Equilibrium':
                output = Equilibrium(junction, output_info=output_info)
                solar_cell[i].equilibrium = output
            elif task == 'ShortCircuit':
                output = ShortCircuit(junction, wavelengths=wavelength, photon_flux=bs, output_info=output_info, rs=rs)
                solar_cell[i].short_circuit = output
            elif task == 'IV':
                output = IV(junction, vfin=vfin, vstep=vstep, output_info=output_info, IV_info=IV_info, rs=rs,
                            escape=escape, light=light, wavelengths=wavelength, photon_flux=bs)
                solar_cell[i].iv = output
            elif task == 'QE':
                output = QE(junction, wavelengths=wavelength, photon_flux=bs, rs=rs, output_info=output_info)
                T = bs / bs_initial
                output['QE']['EQE'] = output['QE']['IQE'] * (1 - R) * T
                solar_cell[i].qe = output
            else:
                raise RuntimeError(
                    'ERROR in the PDD solver: Valid tasks are: "Equilibrium", "ShortCircuit", "IV" and "QE".')

            # And we reduce the amount of light reaching the next junction
            for junction_layer_object in layer_object:
                bs *= np.exp(-junction_layer_object.material.alpha(wavelength) * junction_layer_object.width)

        else:
            raise ValueError("Strange layer-like object discovered in structure stack: {}".format(type(layer_object)))

    solar_cell.R = R
    solar_cell.passive_loss = 1 - passive_loss
    solar_cell.T = bs / bs_initial
    return


# Functions for creating the streucture in the fortran variables and executing the DD solver
def ProcessStructure(device, wavelengths=None, use_Adachi=False):
    """ This function reads a dictionary containing all the device structure, extract the electrical and optical
    properties of the materials, and loads all that information into the Fortran variables. Finally, it initiallise the
    device (in fortran) calculating an initial mesh and all the properties as a function of the possition.

    :param device: A dictionary containing the device structure. See PDD.DeviceStructure
    :param wavelengths: (Optional) Wavelengths at which to calculate the optical properties.
    :param use_Adachi: (Optional) If Adachi model should be use to calculate the dielectric constant of the material.
    :return: Dictionary containing the device structure properties as a function of the position.
    """
    print('Processing structure...')
    # First, we clean any previous data from the Fortran code
    dd.reset()
    output = {}

    if wavelengths is not None:
        output['Optics'] = {}
        output['Optics']['wavelengths'] = wavelengths
        calculate_absorption = True
    else:
        calculate_absorption = False

    # We dump the structure information to the Fotran module and initialise the structure
    i = 0
    first = 0
    last = -1
    while i < device['numlayers']:

        if device['layers'][i]['label'] in ['optics', 'Optics', 'metal', 'Metal']:
            # Optics or metal layers. They are not included in the PDD solver and are just used for optics
            if i == first:
                first = i + 1
                i = i + device['layers'][i]['numlayers']
                continue
            else:
                last = i - device['numlayers'] - 1
                break

        if device['layers'][i]['group'] is None:
            # We have a normal layer
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
                         layer['permittivity'] * Epsi0,
                         layer['radiative_recombination'],
                         layer['electron_auger_recombination'],
                         layer['hole_auger_recombination'],
                         layer['Na'],
                         layer['Nd']]
            dd.addlayer(args=args_list)

            # We load the absorption coeficients if necessary
            if calculate_absorption:
                if 'absorption' in layer.keys():
                    layer['absorption'][1] = np.interp(wavelengths, layer['absorption'][0],
                                                       layer['absorption'][1]).tolist()
                    layer['absorption'][0] = wavelengths.tolist()
                else:
                    layer['absorption'] = LoadAbsorption(device['layers'][i], device['T'], wavelengths,
                                                         use_Adachi=use_Adachi)
                dd.addabsorption(layer['absorption'][1], wavelengths)

        else:
            # We have a group of several layers, usually a QW with 'numlayers' repeated 'repeat' times. 
            for j in range(device['layers'][i]['repeat']):
                for k in range(device['layers'][i]['numlayers']):
                    layer = device['layers'][i + k]['properties']
                    args_list = [layer['width'],
                                 asUnit(layer['band_gap'], 'eV'),
                                 asUnit(layer['electron_affinity'], 'eV'),
                                 layer['electron_mobility'],
                                 layer['hole_mobility'],
                                 layer['Nc'],
                                 layer['Nv'],
                                 layer['electron_minority_lifetime'],
                                 layer['hole_minority_lifetime'],
                                 layer['permittivity'] * Epsi0,
                                 layer['radiative_recombination'],
                                 layer['electron_auger_recombination'],
                                 layer['hole_auger_recombination'],
                                 layer['Na'],
                                 layer['Nd']]
                    dd.addlayer(args=args_list)

                    # We load the absorption coeficients if necessary
                    if calculate_absorption:
                        if 'absorption' in layer.keys():
                            layer['absorption'][1] = np.interp(wavelengths, layer['absorption'][0],
                                                               layer['absorption'][1]).tolist()
                            layer['absorption'][0] = wavelengths.tolist()
                        else:
                            layer['absorption'] = LoadAbsorption(device['layers'][i + k], device['T'], wavelengths,
                                                                 use_Adachi=use_Adachi)
                        dd.addabsorption(layer['absorption'][1], wavelengths)

        i = i + device['layers'][i]['numlayers']

    # We set the surface recombination velocities. This needs to be improved at some to consider other boundary conditions
    dd.frontboundary("ohmic", device['layers'][first]['properties']['sn'], device['layers'][first]['properties']['sp'],
                     0)
    dd.backboundary("ohmic", device['layers'][last]['properties']['sn'], device['layers'][last]['properties']['sp'], 0)

    SetMeshParameters()

    dd.initdevice(mesh_control['meshpoints'])
    print('...done!\n')

    output['Properties'] = DumpInputProperties()
    return output


def Equilibrium(device, output_info=2, wavelengths=None):
    """ Solves the Poisson-DD equations under equilibrium: in the dark with no external current and zero applied voltage. Internally, it calls *ProcessStructure*. Absorption coeficients are not calculated unless *wavelengths* is given as input.

    :param device: A dictionary containing the device structure. See PDD.DeviceStructure
    :param output_info: Indicates how much information must be printed by the fortran solver (1=less, 2=more)
    :param wavelengths: (Optional) Wavelengths at which to calculate the optical properties.
    :return: Dictionary containing the device properties as a function of the position at equilibrium.
    """
    print('Solving equilibrium...')
    output = ProcessStructure(device, wavelengths=wavelengths)
    SetRecombinationParameters(gen=0)
    SetConvergenceParameters()

    dd.equilibrium(output_info)
    print('...done!\n')

    output['Bandstructure'] = DumpBandStructure()
    return output


def ShortCircuit(device, wavelengths, photon_flux, rs=0, output_info=1):
    """ Solves the devices electronic properties at short circuit. Internally, it calls Equilibrium.

    :param device: A dictionary containing the device structure. See PDD.DeviceStructure
    :param output_info: Indicates how much information must be printed by the fortran solver (0=min, 2=max)
    :param rs: Series resistance. Default=0
    :param wavelengths: Array with the wavelengths (in m)
    :param photon_flux: Array with the photon_flux (in photons/m2/m) corresponding to the above wavelengths
    :return: A dictionary containing the device properties as a function of the position at short circuit.
    """

    # We run equilibrium
    output = Equilibrium(device, output_info=output_info, wavelengths=wavelengths)

    SetRecombinationParameters(gen=1)

    dd.illumination(photon_flux)
    dd.tmm = 0

    dd.set('rs', rs)
    dd.lightsc(output_info, 1)
    print('...done!\n')

    output['Bandstructure'] = DumpBandStructure()

    return output


def IV(device, vfin, vstep, output_info=1, IV_info=True, rs=0, escape=1, light=False, wavelengths=None,
       photon_flux=None):
    """ Calculates the IV curve of the device between 0 V and a given voltage. Depending if the "sol" parameter is set
    or not, the IV will be calculated in the dark (calling the Equilibrium function) or under illumination (calling
    the ShortCircuit function).

    :param device: A dictionary containing the device structure. See PDD.DeviceStructure
    :param vfin: Final voltage. If it is negative, vstep must also be negative.
    :param vstep: Maximum step size for the IV curve. This is adapted dynamically to ensure that the shape is reproduced correctly.
    :param output_info: Indicates how much information must be printed by the fortran solver (0=min, 2=max)
    :param IV_info: If information about the Voc, Isc and FF should be provided after the calculation. Default=True
    :param rs: Series resistance. Default=0
    :param escape: Indicates if the calculation should stop when Voc is reached (0=False, 1=True). Default=1
    :param light: Indicates if the light IV curve should be calculated
    :param wavelengths: Array with the wavelengths (in m)
    :param photon_flux: Array with the photon_flux (in photons/m2/m) corresponding to the above wavelengths
    :return: A dictionary containing the IV curves, the different components and also the output of Equilibrium or ShortCircuit.
    """
    print('Solving IV...')
    if light:
        output = ShortCircuit(device, wavelengths=wavelengths, photon_flux=photon_flux, rs=rs, output_info=output_info)
        escape = escape
        IV_info = IV_info
    else:
        output = Equilibrium(device, output_info=output_info)
        escape = 0
        IV_info = False

    dd.set('rs', rs)

    dd.runiv(vfin, vstep, output_info, escape)
    print('...done!\n')

    # This is the bandstructure at the last V point, which might or might not be useful
    output['Bandstructure'] = DumpBandStructure()
    output['IV'] = DumpIV(IV_info=IV_info)

    return output


def qe_pdd(solar_cell, args):
    pass


def QE(device, wavelengths, photon_flux, rs=0, output_info=1):
    """ Calculates the quantum efficiency of the device at short circuit. Internally it calls ShortCircuit

    :param device: A dictionary containing the device structure. See PDD.DeviceStructure
    :param wavelengths: Array with the wavelengths (in m)
    :param photon_flux: Array with the photon_flux (in photons/m2/m) corresponding to the above wavelengths
    :param output_info: Indicates how much information must be printed by the fortran solver (0=min, 2=max)
    :param rs: Series resistance. Default=0
    :return: The internal and external quantum efficiencies, in adition to the output of ShortCircuit.
    """
    print('Solving quantum efficiency...')
    output = ShortCircuit(device, wavelengths=wavelengths, photon_flux=photon_flux, rs=rs, output_info=output_info)

    dd.runiqe(output_info)
    print('...done!\n')

    output['QE'] = DumpQE()
    output['QE']['wavelengths'] = output['Optics']['wavelengths']

    return output


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
    output['IQE'] = dd.get('iqe')[0:numwl]
    output['IQEsrh'] = dd.get('iqesrh')[0:numwl]
    output['IQErad'] = dd.get('iqerad')[0:numwl]
    output['IQEaug'] = dd.get('iqeaug')[0:numwl]
    output['IQEsurf'] = dd.get('iqesurf')[0:numwl]
    output['IQEsurb'] = dd.get('iqesurb')[0:numwl]

    return output


# ----
# Functions for setting the parameters controling the recombination, meshing and the numerial algorithm
def SetMeshParameters(**kwargs):
    global mesh_control
    # We override the default mesh control, if necessary    
    for key in kwargs.keys():
        if key in mesh_control.keys():
            mesh_control[key] = kwargs[key]

    dd.set('coarse', mesh_control['coarse'])
    dd.set('fine', mesh_control['fine'])
    dd.set('ultrafine', mesh_control['ultrafine'])
    dd.set('growth', mesh_control['growth_rate'])


def SetRecombinationParameters(**kwargs):
    global recombination_control
    # We override the default recombination control, if necessary 
    for key in kwargs.keys():
        if key in recombination_control.keys():
            recombination_control[key] = kwargs[key]

    dd.srh = recombination_control['srh']
    dd.rad = recombination_control['rad']
    dd.aug = recombination_control['aug']
    dd.sur = recombination_control['sur']
    dd.gen = recombination_control['gen']


def SetConvergenceParameters(**kwargs):
    global convergence_control
    # We override the default convergence control, if necessary 
    for key in kwargs.keys():
        if key in convergence_control.keys():
            convergence_control[key] = kwargs[key]

    dd.nitermax = convergence_control['nitermax']
    dd.set('clamp', convergence_control['clamp'])
    dd.set('atol', convergence_control['ATol'])
    dd.set('rtol', convergence_control['RTol'])
