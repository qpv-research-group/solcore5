from collections import defaultdict
from solcore import ParameterSystem
import solcore


class Structure(list):
    """ Subclass of lists that stores information of a 'sample' consisting of several 'layers'."""

    def __init__(self, *args, **kwargs):
        super(Structure, self).__init__(*args)
        self.__dict__.update(kwargs)
        self.labels = [None] * len(self)

    def append(self, new_layer, layer_label=None, repeats=1):
        # Pass the arguments to the superclass for extending
        for i in range(repeats):
            # Extend the structure labels
            self.labels.append(layer_label)
            super(Structure, self).append(new_layer)

    def append_multiple(self, layers, layer_labels=None, repeats=1):

        assert type(layers) == type([]), "`append_multiple` only accepts lists for the first argument."

        if layer_labels is not None:
            assert len(layers) == len(
                layer_labels), "When using `layer_labels` keyword a label must be specified for each layer added i.e. layers and layer_labels must have the same number of elements.  Either fix this or simply do not assign any labels (i.e. layer_labels=None)."

        for i in range(repeats):
            # Extend the structure by appending layers
            self.extend(layers)

            # Extend the structure labels by appending an equal number of None values
            # or by appending the actual labels.
            if layer_labels is None:
                self.labels.extend([None] * len(layers))
            else:
                self.labels.extend(layer_labels)

    def __str__(self):

        layer_info = ["  {} {}".format(
            layer,
            self.labels[i] if self.labels[i] is not None else "",
        ) for i, (layer, label) in enumerate(zip(self, self.labels))]

        return "<Structure object\n{}\n{}>".format(str(self.__dict__), "\n".join(layer_info))

    def width(self):
        return sum([layer.width for layer in self])

    def relative_widths(self):
        total = 0
        aggregate_widths = defaultdict(float)
        for layer, comment in zip(self, self.labels):
            aggregate_widths[comment] += layer.width
            total += layer.width

        for layername in aggregate_widths.keys():
            aggregate_widths[layername] /= total

        return aggregate_widths


class Layer:
    """ Class that stores the information about layers of materials, such as thickness and composition.
    It is the building block of the 'Structures' """

    def __init__(self, width, material, role=None, geometry=None, **kwargs):
        """ Layer class constructor.

        :param width: Width of the layer, in SI units.
        :param material: Solcore material
        :param role: Role of the layer.
        :param kwargs: Any other keyword parameter which will become part of the layer attributes
        """
        self.width = width
        self.material = material
        self.role = role
        self.geometry = geometry
        self.__dict__.update(**kwargs)

    def __str__(self):
        """ Representation of the Layer object
        :return: A string with a summary of the layer properties
        """
        widthstring = "{:.3}nm".format(self.width * 1e9)
        return "<{}layer {} {}>".format(
            self.role + " " if self.role != None else "",
            widthstring,
            self.material,
        )


class Junction(list):
    """ Class that groups together the layers that make a junction. Esentially, it is just a list with attributes you can update.
    """

    def __init__(self, *args, **kwargs):
        self.__dict__.update(kwargs)
        list.__init__(self, *args)

    def __str__(self):
        layer_info = ["{}".format(layer) for layer in self]
        return "<Junction object \n\t{}\n\t{}>".format(str(self.__dict__), "\n\t".join(layer_info))


class TunnelJunction(Junction):
    """ Class that contains the minimum definitions for a tunnel junction, ie. a series resistance in our case.
    """
    def __init__(self, *args, **kwargs):
        Junction.__init__(self, *args, **kwargs)

        self.R = kwargs['R'] if 'R' in kwargs.keys() else 1e-16
        self.pn = kwargs['pn'] if 'pn' in kwargs.keys() else True


# CONVERSION UTILITIES
# The following functions are used to move back and forth between the Solcore structures and the Device structures used
# in the PDD solver
def InLineComposition(layer):
    """ Hack to use the Adachi-alfa calculator, provinding the composition as a single string

    :param layer: A layer as defined in the Device structures of the PDD solver
    :return: A mterial string
    """
    comp = layer['properties']['composition']
    if 'element' in comp.keys():
        return comp['material'].replace(comp['element'], comp['element'] + str(comp['fraction']))
    else:
        return comp['material']


def SolcoreMaterialToStr(material_input):
    """ Translate a solcore material into something than can be easily stored in a file and read

    :param material_input: A solcore material
    :return: A dictionary with the name, consituents and composition of the material
    """
    material_string = material_input.__str__().strip('<>').split(" ")
    material_name = material_string[0].strip("'")
    composition = {'material': material_name}

    alloy = True if len(material_input.composition) > 0 else False

    if alloy:
        material_composition = material_string[2].split("=")
        for i, comp in enumerate(material_composition):
            if comp in material_name:
                composition['element'] = material_composition[i]
                composition['fraction'] = float(material_composition[i + 1])

    return composition


def ToSolcoreMaterial(comp, T, execute=False, **kwargs):
    """ It provides a solcore material out of its string composition. The output can be a string with the command or
    a solcore material itself.

    :param comp: A Solcore material as a string
    :param T: The temperature
    :param execute: If a Solcore material must be provided or just a string that would do that.
    :return: A Solcore material or a string with the command to calculate the solcore material
    """
    # It provides a solcore material out of its string composition. The output can be a string with the comand or a solcore material itself.
    if 'element' in comp.keys():
        # A ternary material
        out = 'solcore.material("%s")(T=%s, %s=%s' % (comp['material'], T, comp['element'], comp['fraction'])
    else:
        # A binary material
        out = 'solcore.material("%s")(T=%s' % (comp['material'], T)

    for key in kwargs.keys():
        out = out + ', {}={}'.format(key, kwargs[key])

    out = out + ') '

    if execute:
        return eval(out)
    else:
        return out


def ToLayer(width, material, role):
    """ Creates a Layer based on strings containing the width, the material and the role

    :param width: Width of the layer
    :param material: Material of the layer, as a string
    :param role: The role of the layer
    :return: A Solcore Layer
    """
    return eval('Layer( width=%s, material=%s, role="%s"  )' % (width, material, role))


def ToStructure(device):
    """ Translate the device object (as used by the PDD solver) into a list of solcore Structure that can be used in
    other plugings. Only translate the composition, for now. It works only on non-nested device structures
    (QW, for example, but not in devices with QW)

    :param device: A solcore device structure as used by the PDD solver
    :return: A Solcore Structure
    """

    LayersList = []
    MatList = []

    for i in range(device['numlayers']):
        layer = device['layers'][i]
        MatList.append(ToSolcoreMaterial(layer['properties']['composition'], device['T']))
        LayersList.append(ToLayer(layer['properties']['width'], MatList[i], layer['label']))
        LayersList[-1].material.strained = 'True'

    LayersList = Structure(LayersList)
    LayersList.substrate = ToSolcoreMaterial(device['substrate'], device['T'], execute=True)

    return LayersList
