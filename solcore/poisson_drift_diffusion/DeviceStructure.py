import sys

import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quad

from solcore import constants, material
from solcore.absorption_calculator import adachi_alpha
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
                     'relative_permittivity': DefaultMaterial.relative_permittivity,  # relative to epsilon0
                     'electron_auger_recombination': 1e-42,  # m6 s-1,
                     'hole_auger_recombination': 1e-42,  # m6 s-1
                     # 'radiative_recombination': 7.2e-16,  # m3 s-1
                     'radiative_recombination': DefaultMaterial.radiative_recombination,  # m3 s-1
                     'Nd': 1,  # m-3
                     'Na': 1}  # m3


def CreateDeviceStructure(name=None, role='device', T=293, layers=None, comments='', repeat=1, substrate=DefaultMaterial,
                          reflection=None):
    """ Creates the structure of a device to be used in drift diffusion calculations.

    :param name: Name of the structure
    :param role: Role of the structure (eg. "device", "MQW", etc)
    :param T: Temperature
    :param layers: List containing all the layers of the structure. It can also be a Structure or a Junction.
    :param comments: A string with comments
    :param repeat: If the structure is actually something that is repeated many times, for example a quantum well
    :param substrate: A Solcore material defining the substrate of the structure.
    :param reflection: The reflection of the structure -> Not use from Solcore v4
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
    """ Loads the absorption coefficient for the material of a given layer. It is only useful when solving QWs, since it is necessary to calculate an effective absorption coefficient before solving the opticalproperties of the structure using the BL, TMM or RCWA solvers.

    :param layer: The layer we want to get the absorption coefficient from.
    :param T: The temperature
    :param wavelengths: The wavelength we want the absorption at
    :param use_Adachi: If we should use the Adachi method to get the absorption coeficient.
    :return: A list with two columns, the wavelengths and the absorption coeficient.
    """

    if use_Adachi:
        try:
            # 0 = Energy, 1 = n, 2 = k, 3 = Absorption
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
            # 0 = Energy, 1 = n, 2 = k, 3 = Absorption
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

        # In the end, we convert the absorption coefficient in extinction coefficient
        kk = QW[i].material.absorption * QW.wl / 4 / np.pi
        layer_mat.k = interp1d(QW.wl, kk, bounds_error=False, fill_value=(0, 0))
        layer_mat.n = QW[i].material.n
        # layer_mat.alpha = interp1d(QW.wl, QW[i].material.absorption, bounds_error=False, fill_value=(0, 0))
        import matplotlib.pyplot as plt
        plt.semilogy(QW.wl*1e9, layer_mat.n(QW.wl), label=i)
        plt.semilogy(QW.wl*1e9, layer_mat.k(QW.wl))

        # And the radiative recombination parameter
        inter = lambda E: layer_mat.n(E) ** 2 * layer_mat.alphaE(E) * np.exp(-E / (kb * T)) * E ** 2
        upper = layer_mat.band_gap + 10 * kb * T
        Br = 1.0 / layer_mat.ni ** 2 * 2 * pi / (h ** 3 * c ** 2) * quad(inter, 0, upper)[0]
        layer_mat.radiative_recombination = Br

        # And add the layer to the list of layers
        new_QW.append(Layer(width, layer_mat))

    # As the QW might be actually a MQW, we repeat this as many times as needed
    new_QW = N * new_QW

    plt.ylim(1e-5, 10)
    plt.legend()
    plt.show()
    import sys
    sys.exit()

    return new_QW


def CalculateAbsorptionProfile(z, wl, absorption):
    out = np.array(wl)
    out = np.vstack((out, absorption(z)))

    return out
