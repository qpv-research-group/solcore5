import copy
import numpy as np

from solcore import material
from solcore.solar_cell import SolarCell
from solcore.structure import Layer, Junction, TunnelJunction
from solcore.solar_cell_solver import solar_cell_solver

from .properties_and_settings import *


def create_reflectivity(filename):
    """

    :param filename:
    :return:
    """
    return None


def run_solar_cell_model(task, model):
    """

    :param model:
    :return: a Solar cell object
    """

    # First, we take care of the options
    # ----------------------------------
    options = copy.copy(model['Global'])
    options.update(model['Electrical'])
    options.update(model['Optical'])
    T = float(options['T'])

    # We have to validate the options and check they have the correct format and type
    array_based_options = ['voltages', 'internal_voltages', 'position', 'wavelength']
    float_options = ['T', 'T_ambient']
    bool_options = ['mpp', 'light_iv', 'radiative_coupling']

    try:
        for p in options:
            # For the options that should be arrays, we create the arrays. The wavelength need ot be converted to meters
            if p in array_based_options:
                c = options[p].split(', ')
                ini = float(c[0])
                end = float(c[1])
                num = int(c[2])
                options[p] = np.linspace(ini, end, num) * 1e-9 ** (p == 'wavelength')

            # For the options that need to be floats
            elif p in float_options:
                options[p] = float(options[p])

            elif p in bool_options:
                options[p] = options[p].lower() in ['true', 'yes', 't', 1]

    except Exception as err:
        print('ERROR parsing the solver options: Option format not recognised')
        print(err)
        raise

    # Now we create the layers and junctions to the solar cell object
    # ---------------------------------------------------------------
    sc_data = model['Solar cell']
    all_layers = []
    all_materials = []

    for i in sc_data['structure']:
        current_element = sc_data[i[0]]

        # First the individual layers
        if 'Layer' in current_element['type']:
            layer_properties = {}
            width = current_element['width'] * 1e-9

            # Set the composition and get the properties, converting them to the correct units
            for key in current_element['options']:
                if key in ['element', 'x']:
                    layer_properties[current_element['options']['element']] = current_element['options']['x']
                else:
                    layer_properties[key] = current_element['options'][key] * conversion[key]

            all_materials.append(material(current_element['material'])(T=T, **layer_properties))

            all_layers.append(Layer(width, all_materials[-1], role=current_element['name']))
            continue

        # Unless it is an individual layer, we have junctions
        properties = {}
        for p in properties_junctions[current_element['type']]:
            properties[p] = default_junction_properties[p]

        properties.update(**current_element['options'])
        kind = current_element['type'].split('-')[-1]
        if 'TJ' in current_element['type']:
            new_junction = TunnelJunction(name=current_element['name'], kind=kind, T=T, **properties)
        else:
            new_junction = Junction(name=current_element['name'], kind=kind, T=T, **properties)

        # Now we add the layers, if any
        for l in i[1:]:
            current_child = sc_data[l]

            width = current_child['width'] * 1e-9
            layer_properties = {}

            # Set the composition and get the properties, converting them to the correct units
            for key in current_child['options']:
                if key in ['element', 'x']:
                    layer_properties[current_child['options']['element']] = current_child['options']['x']
                else:
                    layer_properties[key] = current_child['options'][key] * conversion[key]

            all_materials.append(material(current_child['material'])(T=T, **layer_properties))
            new_junction.append(Layer(width, all_materials[-1], role=current_child['name']))

        all_layers.append(new_junction)

    # Now we have to create the solar cell
    # ------------------------------------
    reflectivity = create_reflectivity(sc_data['reflectivity']) if sc_data['reflectivity'] != '' else None

    try:
        substrate = material(sc_data['substrate'])(T=T, **{sc_data['element']: sc_data['composition']})
    except KeyError:
        substrate = material(sc_data['substrate'])(T=T)

    sc = SolarCell(layers=all_layers, name=sc_data['name'], substrate=substrate, T=T, cell_area=sc_data['size'],
                   shading=sc_data['shading'], R_series=sc_data['r_series'], reflectivity=reflectivity)

    # With all said and done, we can run the solver
    # ---------------------------------------------
    solar_cell_solver(sc, task, user_options=options)
    print('Done!')

    return sc
