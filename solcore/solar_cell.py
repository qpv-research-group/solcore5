from solcore.structure import Layer, Junction, Structure, TunnelJunction
from solcore import material, si

import numpy as np


def default_GaAs(T):
    # We create the other materials we need for the device
    window = material('AlGaAs')(T=T, Na=5e24, Al=0.8)
    p_GaAs = material('GaAs')(T=T, Na=1e24)
    n_GaAs = material('GaAs')(T=T, Nd=8e22)
    bsf = material('GaAs')(T=T, Nd=2e24)

    output = Junction([Layer(width=si('30nm'), material=window, role="Window"),
                       Layer(width=si('150nm'), material=p_GaAs, role="Emitter"),
                       Layer(width=si('3000nm'), material=n_GaAs, role="Base"),
                       Layer(width=si('200nm'), material=bsf, role="BSF")], sn=1e6, sp=1e6, T=T, kind='PDD')

    return output


class SolarCell(Structure):
    """ This class is almost identical to the basic Structure class in Solcore (it is a subclass of it, actually) but implementing some default parameter values and control about the types of layers. It should work anywhere where a Structure object works.
    """

    def __init__(self, layers=None, T=298, cell_area=1, reflectivity=None, shading=0, substrate=None, incidence=None,
                 R_series=0, **kwargs):
        """ Constructor of the class.

        :param layers: A list with the layers to add. The layers might include individual Layer objects of whole Junction objects.
        :param T: Temperature.
        :param cell_area: The area of the cell.
        :param reflectivity: Function that calculates the reflectivity as a function of energy.
        :param shading: Shading losses due to the front metal contacts.
        :param substrate: Substrate of the solar cell.
        :param incidence: Material above the solar cell. Nominally air.
        :param R_series: Series resistance of the structure
        :param kwargs: Other possible attributes.
        """

        assert isinstance(layers, list), "Layers must be provided inside a list, even if it is just one layer."

        super(SolarCell, self).__init__(layers)
        self.__dict__.update(kwargs)

        self.T = T
        self.cell_area = cell_area
        self.shading = shading
        self.reflectivity = reflectivity
        self.junctions = 0
        self.junction_indices = []
        self.tunnel_indices = []
        self.substrate = substrate
        self.incidence = incidence
        self.R_series = R_series

        for i, element in enumerate(layers):
            self.sort_layer_type(element, i)
            self.R_series += element.R_series if hasattr(element, 'R_series') else 0

    def sort_layer_type(self, layer, i):
        """ Sorts the layer in different categories, depending on its type, and keeps record on the indices of that type of layer.

        :param layer: The layer to check
        :param i: The index of that layer
        :return: None
        """

        if type(layer) == Junction:
            self.junction_indices.append(i)
            self.junctions += 1
        if type(layer) == TunnelJunction:
            self.tunnel_indices.append(i)

    def append(self, new_layer, layer_label=None, repeats=1):
        """ Appends a layer to the structure a certain number of times to the structure.

        :param new_layer: The layer to append.
        :param layer_label: An optional label for that layer.
        :param repeats: Number of times to repeat add the layer.
        :return: None
        """

        for i in range(repeats):
            self.sort_layer_type(new_layer, len(self))
            super(SolarCell, self).append(new_layer, layer_label)

    def append_multiple(self, layers, layer_labels=None, repeats=1):
        """ Appends multiple layers a certain umber of times to the structure.

        :param layers: A list with the layers to append.
        :param layer_labels: An optional list with the labels. If present, it must have the same length that the layers list.
        :param repeats: Number of times to add this set of layers.
        :return: None
        """

        assert isinstance(layers, list), "'append_multiple' only accepts lists for the first argument."

        if layer_labels is not None:
            assert len(layers) == len(layer_labels), "When using 'layer_labels' keyword a label must be specified for " \
                                                     "each layer added i.e. layers and layer_labels must have the same " \
                                                     "number of elements. Either fix this or simply do not assign any " \
                                                     "labels (i.e. layer_labels=None)."
        else:
            layer_labels = [None] * len(layers)

        for i in range(repeats):
            for j, element in enumerate(layers):
                self.append(element, layer_labels[j])

    def update_junction(self, junction, **kwargs):
        """ Adds or updates the attributes - not the layers - of a junction.

        :param junction: The junction to update.
        :param kwargs: The attributes to update.
        :return: None
        """

        try:
            num = self.junction_indices[junction]
            self[num].__dict__.update(kwargs)
        except IndexError:
            print('ERROR updating junction: The junction index must be {} or less.'.format(len(self.junction_indices)))

    def __call__(self, i):
        return self[self.junction_indices[i]]


if __name__ == '__main__':
    window = material('AlGaAs')(T=298, Na=1e24, Al=0.8)
    stack = [Layer(width=si("50nm"), material=window),
             default_GaAs(298)]

    # stack = [default_GaAs(298)]

    my_cell = SolarCell(layers=stack)
    #
    # my_cell.append_multiple([default_GaAs(298), default_GaAs(298), default_GaAs(298)])
    # print(my_cell)

    from solcore.poisson_drift_diffusion.DriftDiffusionUtilities import solve_pdd, default_photon_flux, \
        default_wavelengths
    import matplotlib.pyplot as plt

    solve_pdd(my_cell, 'QE', vfin=1.2, vstep=0.05, light=True)
    QE = my_cell(0).qe

    # Finally, we plot the internal and external quantum efficiencies using the information stored in the output dictionaries
    plt.plot(QE['QE']['wavelengths'] / 1e-9, QE['QE']['IQE'] * 100, label='IQE')
    plt.plot(QE['QE']['wavelengths'] / 1e-9, QE['QE']['EQE'] * 100, label='EQE')
    plt.plot(QE['QE']['wavelengths'] / 1e-9, my_cell.T * 100, label='T')
    plt.ylim(0, 100)
    plt.legend(loc='lower left')
    plt.ylabel('QE (%)')
    plt.xlabel('Wavelength (nm)')

    plt.show()
