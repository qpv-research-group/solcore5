import copy
import os
import sys

import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quad

from solcore import constants, material
from solcore.absorption_calculator import adachi_alpha, calculate_absorption_profile, calculate_rat
from solcore.material_data import calculate_mobility
from solcore.structure import *
from .QWunit import QWunit

Epsi0 = constants.vacuum_permittivity
q = constants.q
pi = constants.pi
h = constants.h
kb = constants.kb
m0 = constants.electron_mass
vacuum_permittivity = constants.vacuum_permittivity
c = constants.c

# The default material is a set of minimum properties that are used for the layers unless otherwise stated
DefaultMaterial = material("GaAs")(T=293)
DefaultProperties = {'band_gap': DefaultMaterial.band_gap,  # J
                     'electron_affinity': DefaultMaterial.electron_affinity,  # J
                     'eff_mass_electron_Gamma': DefaultMaterial.eff_mass_electron_Gamma,  # relative to m0
                     'eff_mass_hh_z': DefaultMaterial.eff_mass_hh_z,  # relative to m0
                     'eff_mass_lh_z': DefaultMaterial.eff_mass_lh_z,  # relative to m0
                     'electron_mobility': calculate_mobility("GaAs", 0, 1),  # m2 V-1 s-1
                     'hole_mobility': calculate_mobility("GaAs", 1, 1),  # m2 V-1 s-1
                     'ni': DefaultMaterial.ni,  # m-3
                     'Nc': DefaultMaterial.Nc,  # m-3
                     'Nv': DefaultMaterial.Nv,  # m-3
                     'electron_minority_lifetime': 3e-6,  # s
                     'hole_minority_lifetime': 2.5e-7,  # s
                     'permittivity': 12.9,  # relative to epsilon0
                     'electron_auger_recombination': 1e-42,  # m6 s-1,
                     'hole_auger_recombination': 1e-42,  # m6 s-1
                     # 'radiative_recombination': 7.2e-16,  # m3 s-1
                     'radiative_recombination': DefaultMaterial.radiative_recombination,  # m3 s-1
                     'Nd': 1,  # m-3
                     'Na': 1,  # m3
                     'sn': 1e6,  # m s-1
                     'sp': 1e6}  # m s-1


def CreateDeviceStructure(name, role='device', T=293, layers=None, comments='', repeat=1, substrate=DefaultMaterial,
                          reflection=None):
    """ Creates the structure of a device to be used in drift diffusion calculations.
    
    This class defines an object that can be used in drift diffusion calculations. 
    It has the advantage over the old Sample class (based on the Structure class) of being organised in a format that than be easily stored
    using Json or managed by the upcoming front end, Sunglasses.
    
    With minimum coding effort, it should be possible to use it with the Schrodinger solver, analythical methods, etc.

    TO BE DONE: This should be replaced by a state_object, the same way it is done in the analytical solar cell solvers

    :param name: Name of the structure
    :param role: Role of the structure (eg. "device", "MQW", etc)
    :param T: Temperature
    :param layers: List containing all the layers of the structure
    :param comments: A string with comments
    :param repeat: If the structure is actually something that is repeated many times, for example a quantum well
    :param substrate: A solcore material defining the substrate of the structure.
    :param reflection: The reflexion of the structure
    :return: A dictionary-like object with all the information related to the device
    
    """
    output = {}
    output['name'] = name
    output['role'] = role
    output['T'] = T
    output['numlayers'] = 0
    output['comments'] = comments
    output['repeat'] = repeat
    output['layers'] = []
    output['substrate'] = SolcoreMaterialToStr(substrate)
    output['reflection'] = reflection
    output['absorption'] = layers.absorbed if hasattr(layers, 'absorbed') else lambda x: 0
    output['sn'] = layers.sn if hasattr(layers, 'sn') else 1e6
    output['sp'] = layers.sp if hasattr(layers, 'sp') else 1e6

    AddLayers(output, layers)

    return output


def AddLayers(device, layers):
    """ Add layers to the structure

    :param device: The device structure
    :param layers: A list with the layers to add
    :return: None
    """
    for layer in layers:
        NewLayer = {}

        NewLayer['label'] = layer.role
        NewLayer['class'] = 'Layer'
        NewLayer['group'] = None
        NewLayer['numlayers'] = 1
        NewLayer['repeat'] = 1
        NewLayer['properties'] = GetLayerProperties(layer, device['T'])

        device['layers'].append(NewLayer)
        device['numlayers'] = device['numlayers'] + 1


def RemoveLayer(device, i):
    """ Remove a layer from the device structure

    :param device: the device structure
    :param i: the index of the layer to remove
    :return: None
    """
    del device['layers'][i]
    device['numlayers'] = device['numlayers'] - 1


def GetLayerProperties(layer, T):
    """ Given a layer, get all the properties of that layer and the material it is made of at the given temperature

    :param layer: The layer of interest
    :param T: The temperature
    :return: Dictionary with all the properties
    """
    # Get all the properties of the new layer. If they dont exist in the layer definition, the default values are used. 
    NewProperties = {}
    NewProperties['composition'] = SolcoreMaterialToStr(layer.material)
    NewProperties['width'] = layer.width

    for key in DefaultProperties.keys():
        try:
            NewProperties[key] = getattr(layer.material, key)
        except ValueError:
            NewProperties[key] = DefaultProperties[key]

    return NewProperties


def LoadAbsorption(layer, T, wavelengths, use_Adachi=False):
    """ If there is a file containing the absorption in this layer, we try to load it. It must have at least two
    columns: wavelengths and absorption coefficient. Both must be in SI units. The result is interpolated to the current
    working wavelengths.

    The file can not be loaded: Maybe it does not exists or is not in the correct format. We move forward to
    calculating the absorption. Later, the calculated absorption can be saved with this filename.

    :param layer: The layer we want to get the absorption coefficient from.
    :param T: The temperature
    :param wavelengths: The wavelength we want the absorption at
    :param use_Adachi: If we should use the Adachi method to get the absorption coeficient.
    :return: A list with two columns, the wavelnegths and the absorption coeficient.
    """
    if 'absorption_file' in layer.keys():

        try:
            data = np.loadtxt(layer['absorption_file'])
            absorption = np.interp(wavelengths, data[:, 0], data[:, 1])
            return [wavelengths.tolist(), absorption.tolist()]
        except:
            print('Error loading the absorption file or layer {}'.format(layer))
            sys.exit()

    if use_Adachi:
        try:
            # 0 = Energy, 1 = n, 2 = k.txt, 3 = Absorption
            absorption = adachi_alpha.create_adachi_alpha(InLineComposition(layer), T=T, wl=wavelengths)[3]
        except:
            print("Warning: Using experimental data to estimate the absorption coefficient of material: ",
                  InLineComposition(layer))
            absorption = ToSolcoreMaterial(layer['properties']['composition'], T, execute=True).alpha(wavelengths)
            if layer['properties']['composition']['material'] == 'InGaAs':
                print(
                    "Warning: Extrapolation of experimental absorption data for InGaAs is not reliable at longer wavelengths.")
                print("    >>>: We truncate the absorption at the bandgap wavelength.")
                edge = 1240e-9 / (layer['properties']['band_gap'] / q)
                edgeidx = np.abs(wavelengths - edge).argmin()
                absorption[edgeidx:] = 0

    else:
        print(layer['properties']['composition'])
        absorption = ToSolcoreMaterial(layer['properties']['composition'], T, execute=True).alpha(wavelengths)
        try:

            if layer['properties']['composition']['material'] == 'InGaAs':
                print(
                    "Warning: Extrapolation of experimental absorption data for InGaAs is not reliable at longer wavelengths.")
                print("    >>>: We truncate the absorption at the bulk bandgap wavelength.")
                edge = 1240e-9 / (layer['properties']['band_gap'] / q)
                edgeidx = np.abs(wavelengths - edge).argmin()
                absorption[edgeidx:] = 0
        except Exception as err:
            print("Warning: Using Adachi calculation to estimate the absorption coefficient of material: ",
                  InLineComposition(layer))
            # 0 = Energy, 1 = n, 2 = k.txt, 3 = Absorption
            try:
                absorption = adachi_alpha.create_adachi_alpha(InLineComposition(layer), T=T, wl=wavelengths)[3]
            except:
                print("Warning: No absorption information found for material {}. Setting it equal to zero.".format(
                    InLineComposition(layer)))
                absorption = 0 * wavelengths

    return [wavelengths.tolist(), absorption.tolist()]


def SolveQWproperties(device, calculate_absorption=True, WLsteps=(300e-9, 1100e-9, 201), wavelengths=None,
                      periodic=True, filter_strength=0.0, blur=None, blurmode="left", mode='kp8x8_bulk',
                      use_Adachi=False,
                      alpha_params=None):
    """ Considers the device as a QW and solves its properties, including the modification of the bandeges due to strain, the efective mases and the absorption coefficient. Without calling this function, the structure is just a colection of layers with bulk-like properties.

    :param device: The device structure
    :param calculate_absorption: If absorption must be calculated
    :param WLsteps: wavelengths in which to calculate the absorption (input for np.linspace function)
    :param wavelengths: An array with the waveengths
    :param periodic: If it has to be assumed that the structure is perdiodic
    :param filter_strength:
    :param blur:
    :param blurmode:
    :param mode:
    :param use_Adachi:
    :param alpha_params:
    :return: A dictionary with the output of the Schrodinger solver.
    """
    print('Solving QW properties...')
    T = device['T']

    QW = QWunit(ToStructure(device), substrate=ToSolcoreMaterial(device['substrate'], device['T'], execute=True))
    output = QW.solve(calculate_absorption=calculate_absorption, WLsteps=WLsteps, wavelengths=wavelengths,
                      T=device['T'], periodic=periodic, filter_strength=filter_strength, blur=blur, blurmode=blurmode,
                      mode=mode, use_Adachi=use_Adachi, alpha_params=alpha_params)

    for i in range(len(QW)):
        device['layers'][i]['properties']['band_gap'] = QW[i].eff_band_gap
        device['layers'][i]['properties']['electron_affinity'] = QW[i].eff_electron_affinity
        device['layers'][i]['properties']['eff_mass_electron_Gamma'] = QW[i].material.eff_mass_electron_Gamma
        device['layers'][i]['properties']['eff_mass_hh_z'] = QW[i].material.eff_mass_hh_z
        device['layers'][i]['properties']['eff_mass_lh_z'] = QW[i].material.eff_mass_lh_z
        device['layers'][i]['properties']['Nc'] = QW[i].material.Nc
        device['layers'][i]['properties']['Nv'] = QW[i].material.Nv
        device['layers'][i]['properties']['ni'] = np.sqrt(
            QW[i].material.Nc * QW[i].material.Nv * np.exp(-QW[i].eff_band_gap / (kb * T)))

        if calculate_absorption:
            device['layers'][i]['properties']['absorption'] = [QW.wl.tolist(), QW[i].material.absorption.tolist()]

    # Finally, we re-build a list of layers with the effective properties
    N = device['repeat']
    new_QW = []

    for i in range(len(QW)):
        # First, we create a dictionary with all the updated parameters
        param = dict(device['layers'][i]['properties'])
        del param['absorption']
        del param['composition']
        del param['width']

        # We recover the composition and thickness
        mat = device['layers'][i]['properties']['composition']
        width = device['layers'][i]['properties']['width']

        # Create the material with the updated properties
        layer_mat = ToSolcoreMaterial(mat, T, execute=True, **param)

        # In the end, we convert the absorption coeficient in extinction coefficient
        kk = QW[i].material.absorption * QW.wl / 4 / np.pi
        layer_mat.k = interp1d(QW.wl, kk, bounds_error=False, fill_value=(0, 0))
        # layer_mat.alpha = interp1d(QW.wl, QW[i].material.absorption, bounds_error=False, fill_value=(0, 0))

        # And the radiative recombination parameter
        inter = lambda E: layer_mat.n(E) ** 2 * layer_mat.alphaE(E) * np.exp(-E / (kb * T)) * E ** 2
        upper = layer_mat.band_gap + 10 * kb * T
        Br = 1.0 / layer_mat.ni ** 2 * 2 * pi / (h ** 3 * c ** 2) * quad(inter, 0, upper)[0]
        layer_mat.radiative_recombination = Br

        # And add the layer to the list of layers
        new_QW.append(Layer(width, layer_mat))

    # As the QW might be actually a MQW, we repeat this as many times as needed
    new_QW = N * new_QW

    return new_QW


def CalculateAbsorptionProfile(z, wl, absorption):
    out = np.array(wl)
    out = np.vstack((out, absorption(z)))

    return out

    #
    # def calculate_reflection(device, wavelengths):
    #     """ Calculates the reflexion of a device structure
    #
    #     :param device: The device structure
    #     :param wavelengths: The wavelengths
    #     :return: The reflexion at the given wavelengths
    #     """
    #     if device['reflection'] is not None:
    #         try:
    #             # We asumme that we have a file with the reflection
    #             print('Loading reflection from file %s ...' % (device['reflection']))
    #             data = np.loadtxt(device['reflection'])
    #             R = np.interp(wavelengths, data[:, 0], data[:, 1])
    #             print('...Sucess!!')
    #         except:
    #             # If it fails, we calculate it based on the refractive index of the first layer
    #             comp = device['layers'][0]['properties']['composition']
    #             n = ToSolcoreMaterial(comp, device['T'], execute=True).n(wavelengths)
    #             R = ((1 - n) / (1 + n)) ** 2
    #             print(
    #                 '... Error!! Device surface reflection calculated from the refractive index of the first layer, instead. ')
    #     else:
    #         # Otherwise, we calculate it based on the refractive index of the first layer
    #         comp = device['layers'][0]['properties']['composition']
    #         n = ToSolcoreMaterial(comp, device['T'], execute=True).n(wavelengths)
    #         k = ToSolcoreMaterial(comp, device['T'], execute=True).k(wavelengths)
    #         nc = n + k * 1.0j
    #         R = np.abs((1 - nc) / (1 + nc)) ** 2
    #         print('Device surface reflection calculated from the refractive index of the first layer. ')
    #
    #     return R
    #
    # def calculate_optics(device, wavelengths, dist=None):
    #     """ Uses the transfer matrix solver to calculate the optical properties of the structure: that is, the reflection
    #     and the absorption as a function of the position.
    #
    #     :param device: A device structure
    #     :param wavelengths: The wavelengths at which to calculate the optical information (in m)
    #     :param dist: The positions at which to calculate the absorption (in m). If None, it is calculated internally.
    #     :return: A dictionary with the reflection, the position, the wavelengths and the absorption as a function of
    #     the wavelength and position.
    #     """
    #
    #     output = {}
    #     output['wavelengths'] = wavelengths
    #     wl = wavelengths * 1e9  # Input is in meters but the calculators use nm
    #     if dist is None:
    #         d = dist
    #     else:
    #         d = dist * 1e9  # Input is in meters but the calculators use nm
    #
    #     rat = calculate_rat(device, wl)
    #     output['R'] = rat['R']
    #
    #     absorption = calculate_absorption_profile(device, wl, dist=d)
    #
    #     output['absorption'] = absorption['absorption'] * 1e9
    #     output['position'] = absorption['position'] * 1e-9
    #
    #     optics_thickness = 0
    #     for layer in device['layers']:
    #         if layer['label'] in ['optics', 'Optics']:
    #             optics_thickness += layer['properties']['width']
    #         else:
    #             break
    #
    #     output['position'] -= optics_thickness
    #
    #     return output


    #
    #
    # def Load(filename, yaml=False):
    #     """ Loads a device structure stored in a file. By default, the file must be JSON format
    #
    #     :param filename: The filename
    #     :param yaml: If the format is YALM rather than JSON
    #     :return: A device structure
    #     """
    #     if yaml:
    #         try:
    #             # If chosen, we try to use yaml. If anything fails, we try with json.
    #             import yaml
    #             stream = open(filename + '.yaml', 'r')
    #             output = yaml.load(stream)
    #         except:
    #             import json
    #             stream = open(filename + '.json', 'r')
    #             output = json.load(stream)
    #     else:
    #         import json
    #         stream = open(filename + '.json', 'r')
    #         output = json.load(stream)
    #
    #     stream.close()
    #
    #     for i in range(output['numlayers']):
    #         if 'absorption_file' in output['layers'][i].keys():
    #             # If there is a file containing the absorption in this layer, we try to load it. It must have at least to columns: wavelengths and absorption coefficient. Both must be in SI units. The result is interpolated to the current working wavelengths.
    #             try:
    #                 data = np.loadtxt(output['layers'][i]['absorption_file'])
    #                 output['layers'][i]['properties']['absorption'] = [data[:, 0].tolist(), data[:, 1].tolist()]
    #             except:
    #                 # The file can not be loaded. Maybe it does not exists or is not in the correct format. We move forward to calculating the absorption. Later, the calculated absorption can be saved with this filename.
    #                 print('Error loading absorption file %s. The absorption will be calculated later if necesary. ' % (
    #                     output['layers'][i]['absorption_file']))
    #
    #     return output
    #
    #
    # def Save(device, filename, save_absorptions_individually=False, remove_absorption_from_json=False,
    #          override_absorption=False, directory='default', yaml=False):
    #     """ Save the device structure to a file.
    #
    #     :param device: The device structure to save
    #     :param filename: The filename
    #     :param save_absorptions_individually: If the absorption of the materials must be saved in individual files
    #     :param remove_absorption_from_json: If the absorption must be removed from the JSON file to make it more readable
    #     :param override_absorption: If the external absrption files must override the default absorption
    #     :param directory: Directory in which to save the absorption data
    #     :param yaml: If YALM should be used rather than JSON
    #     :return: None
    #     """
    #     if save_absorptions_individually:
    #         # The absorption coefficient of each layer is saved in an individual file. Optionally, it removes them from the structure so they are not saved also in in the json file
    #         for i in range(device['numlayers']):
    #             if 'absorption' not in device['layers'][i]['properties'].keys(): continue
    #             # if ('absorption_file' in device['layers'][i].keys()) and override_absorption:
    #             #     abs_filename = device['layers'][i]['absorption_file']
    #             # else:
    #             #     abs_filename = '%s_inputs/%s_%s_%s_%s.dat'   %(filename.split('.')[0], device['name'], InLineComposition(device['layers'][i]), device['layers'][i]['label'], i)
    #
    #             if directory == 'default': directory = '%s_inputs' % (filename)
    #
    #             abs_filename = '%s/%s_%s_%s_%s.dat' % (
    #                 directory, device['name'], i, InLineComposition(device['layers'][i]), device['layers'][i]['label'])
    #
    #             os.makedirs(directory, exist_ok=True)
    #             device['layers'][i]['absorption_file'] = abs_filename
    #             np.savetxt(abs_filename, np.transpose(device['layers'][i]['properties']['absorption']))
    #
    #     if remove_absorption_from_json:
    #         for i in range(device['numlayers']):
    #             if 'absorption' not in device['layers'][i]['properties'].keys(): continue
    #             del device['layers'][i]['properties']['absorption']
    #
    #     if yaml:
    #         try:
    #             # If chosen, we try to use yaml. If anything fails, we try with json.
    #             import yaml
    #             stream = open(filename + '.yaml', 'w')
    #             yaml.dump(device, stream, default_flow_style=False)
    #         except:
    #             import json
    #             stream = open(filename + '.json', 'w')
    #             json.dump(device, stream, sort_keys=True, indent=2)
    #     else:
    #         import json
    #         stream = open(filename + '.json', 'w')
    #         json.dump(device, stream, sort_keys=True, indent=2)
    #
    #     stream.close()
    #
    #
    #
    # def OldAddLayers(device, layers):
    #     """ Add layers to the structure
    #
    #     :param device: The device structure
    #     :param layers: A list with the layers to add
    #     :return: None
    #     """
    #     for layer in layers:
    #         NewLayer = {}
    #
    #         if type(layer) is Layer:
    #             NewLayer['label'] = layer.role
    #             NewLayer['class'] = 'Layer'
    #             NewLayer['group'] = None
    #             NewLayer['numlayers'] = 1
    #             NewLayer['repeat'] = 1
    #             NewLayer['properties'] = GetLayerProperties(layer, device['T'])
    #             if 'absorption_file' in layer.material.__dict__.keys():
    #                 NewLayer['absorption_file'] = layer.material.absorption_file
    #
    #             device['layers'].append(NewLayer)
    #             device['numlayers'] = device['numlayers'] + 1
    #
    #         elif type(layer) is dict:
    #             for sublayer in layer['layers']:
    #                 device['layers'].append(copy.deepcopy(sublayer))
    #                 device['numlayers'] = device['numlayers'] + 1
    #
    #                 device['layers'][-1]['group'] = layer['name']
    #                 device['layers'][-1]['numlayers'] = layer['numlayers']
    #                 device['layers'][-1]['repeat'] = layer['repeat']
