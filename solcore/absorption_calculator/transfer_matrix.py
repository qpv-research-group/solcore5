""" This module serves as interface between solcore structures (junctions, layers, materials...) and the
transfer matrix package developed by Steven Byrnes and included in the PyPi repository.

"""
import numpy as np
import solcore
from solcore.interpolate import interp1d
from solcore.structure import ToStructure
from solcore.absorption_calculator import tmm_core_vec as tmm

degree = np.pi / 180


class OptiStack(object):
    """ Class that contains an optical structure: a sequence of layers with a thickness and a complex refractive index.

    It serves as an intermediate step between solcore layers and materials and the stack of thicknesses and
    and n and k.txt values necessary to run calculations involving TMM. When creating an OptiStack object, the thicknesses
    of all the layers forming the Solcore structure and the optical data of the materials of the layers are extracted
    and arranged in such a way they can be easily and fastly read by the TMM functions.

    In addition to a solcore structure with Layers, it can also take a list where each element represent a layer
    written as a list and contains the layer thickness and the dielectrical model, the raw n and k data as a function
    of wavelengths, or a whole Device structure as the type used in the PDD model.

    In summary, this class acepts:

        - A solcore structure with layers
        - A list where each element is [thickness, DielectricModel]
        - A list where each element is [thickness, wavelength, n, k]
        - A list mixing the above:
            [ [thickness, DielectricModel],
              [thickness, wavelength, n, k],
              solcore.Layer,
              solcore.Layer ]

    This allows for maximum flexibility when creating the optical model, allowing to construct the stack with
    experimental data, modelled data and known material properties from the database.

    Yet anther way of defining the layers mixes experimental data with a DielectricModel within the same layer but in
    spectrally distinct regions. The syntaxis for the layer is:

    layer = [thickness, wavelength, n, k, DielectricModel, mixing]

    where mixing is a list containing three elements: [the mixing point (nm), the mixing width (nm),  zero or one]
    depending if the mixing function should be increasing with the wavelength or decreasing. If increasing (zero), the
    Dielectric model will be used at long wavelengths and the experimental data at short wavelengths. If decreasing
    (one) the oposite is done. The mixing point and mixing width control how smooth is the transition between one and
    the other type of data.

    Extra layers such as he semi-infinite, air-like first and last medium, and a back highly absorbing layer are
    included at runtime to fulfill the requirements of the TMM solver or to solve some of its limitations.
    """

    def __init__(self, structure=(), no_back_reflexion=False, substrate=None):
        """ Class constructor. It takes a Solcore structure and extract the thickness and optical data from the
        Layers and the materials. Option is given to indicate if the reflexion from the back of the structure must be
        supressed, usefull for ellipsometry calculations. This is done by creating an artificial highly absorbing but
        not reflecting layer just at the back.

        Alternativelly, it can take simply a list of [thickness, DielectricModel] or [thickness, wavelength, n, k.txt] for
        each layer accounting for the refractive index of the layers. The three options can be mixed for maximum
        flexibility.

        :param structure: A list with one or more layers.
        :param no_back_reflexion: If reflexion from the back must be supressed. Default=False.
        """

        self.widths = []
        self.n_data = []
        self.k_data = []
        self.models = []
        self.layers = []
        self.substrate = substrate

        self.num_layers = 0
        self.add_layers(structure)

        self.no_back_reflexion = no_back_reflexion

    def get_indices(self, wl):
        """ Returns the complex refractive index of the stack.

        :param wl: Wavelength of the light in nm.
        :return: A list with the complex refractive index of each layer, including the semi-infinite front and back
        layers and, opionally, the back absorbing layer used to suppress back surface relfexion.
        """

        out = []
        wl_m = solcore.si(wl, 'nm')

        if hasattr(self, 'n_sub'):
            n1 = self.n_sub(wl_m) + self.k_sub(wl_m)*1.0j
        else:
            if hasattr(wl, 'size'):
                n1 = np.ones_like(wl, dtype=complex)
            else:
                n1 = 1

        if hasattr(wl, 'size'):
            n0 = np.ones_like(wl, dtype=complex)
        else:
            n0 = 1

        for i in range(self.num_layers):
            out.append(self.n_data[i](wl_m) + self.k_data[i](wl_m) * 1.0j)

        # substrate irrelevant if no_back_reflexion = True
        if self.no_back_reflexion:
            return [n0] + out + [self.n_data[-1](wl_m) + self._k_absorbing(wl_m) * 1.0j, n1] # look at last entry in stack,
            # make high;y absorbing layer based on it.

        else:
            return [n0] + out + [n1]


    def get_widths(self):
        """ Returns the widths of the layers of the stack.

        :return: A list with the widths each layer, including the semi-infinite front and back layers and, opionally,
        the back absorbing layer used to suppress back surface relfexion, defined as 1 mm thick.
        """

        if self.no_back_reflexion:
            return [np.inf] + self.widths + [1e6, np.inf]
        else:
            return [np.inf] + self.widths + [np.inf]

    def _k_absorbing(self, wl):
        """ k.txt value of the back highly absorbing layer. It is the maximum between the bottom layer of the stack or a
        finite, small value that will absorb all light within the absorbing layer thickness.

        :param wl: Wavelength of the light in nm.
        :return: The k.txt value at each wavelength.
        """
        return np.maximum(wl / 1e-3, self.k_data[-1](wl))

    @staticmethod
    def _k_dummy(wl):
        """ Dummy k.txt value to be used with the dielectric model, which produces the refractive index as a complex
        number.

        :param wl: Wavelength of the light in nm.
        :return: The k.txt value at each wavelength... set to zero.
        """
        return 0.

    def add_layers(self, layers):
        """ Generic function to add layers to the OptiStack. Internally, it calls add_solcore_layer,
        add_modelled_layer or add_raw_nk_layer.

        :param layers: A list with the layers to add (even if it is just one layer) It can be one or more and it can
        mixed, Solcore-based and modelled layers.
        :return: None
        """
        try:
            # If the input is a whole device structure, we get just the layers information
            if type(layers) is dict:
                layers = ToStructure(layers)

            if self.substrate is not None:
                self.n_sub = self.substrate.n
                self.k_sub = self.substrate.k

            for layer in layers:
                self.layers.append(layer)
                if 'Layer' in str(type(layer)):
                    self._add_solcore_layer(layer)
                elif 'DielectricConstantModel' in str(type(layer[1])):
                    self._add_modelled_layer(layer)
                else:
                    self._add_raw_nk_layer(layer)

                self.num_layers += 1

        except Exception as err:
            print('Error when adding a new layer to the OptiStack. {}'.format(err))

    def remove_layer(self, idx):
        """ Removes layer with index idx from the OptiStack

        :param idx: Index of the layer to remove
        :return: None
        """
        if idx >= self.num_layers:
            print('Error when removing layers. idx must be: 0 <= idx <= {}.'.format(self.num_layers - 1))
            return

        self.widths.pop(idx)
        self.models.pop(idx)
        self.n_data.pop(idx)
        self.k_data.pop(idx)
        self.num_layers -= 1

    def swap_layers(self, idx1, idx2):
        """ Swaps two layers in the OptiStack.

        :param idx1: The index of one of the layers.
        :param idx2: The index of the other.
        :return: None
        """
        if idx1 >= self.num_layers or idx2 >= self.num_layers:
            print('Error when removing layers. idx must be: 0 <= idx1, idx2 <= {}.'.format(self.num_layers - 1))
            return

        self.widths[idx1], self.widths[idx2] = self.widths[idx2], self.widths[idx1]
        self.models[idx1], self.models[idx2] = self.models[idx2], self.models[idx1]
        self.n_data[idx1], self.n_data[idx2] = self.n_data[idx2], self.n_data[idx1]
        self.k_data[idx1], self.k_data[idx2] = self.k_data[idx2], self.k_data[idx1]

    def _add_solcore_layer(self, layer):
        """ Adds a Solcore layer to the end (bottom) of the stack, extracting its thickness and n and k.txt data.

        :param layer: The Solcore layer
        :return: None
        """
        self.widths.append(solcore.asUnit(layer.width, 'nm'))
        self.models.append([])
        self.n_data.append(layer.material.n)
        self.k_data.append(layer.material.k)

    def _add_modelled_layer(self, layer):
        """ Adds a layer to the end (bottom) of the stack. The layer must be defined as a list containing the layer
        thickness in nm and a dielectric model.

        :param layer: The new layer to add as [thickness, DielectricModel]
        :return: None
        """
        self.widths.append(layer[0])
        self.models.append(layer[1])
        self.n_data.append(self.models[-1].n_and_k)
        self.k_data.append(self._k_dummy)

    def _add_raw_nk_layer(self, layer):
        """ Adds a layer to the end (bottom) of the stack. The layer must be defined as a list containing the layer
        thickness in nm, the wavelength, the n and the k data as array-like objects.

        :param layer: The new layer to add as [thickness, wavelength, n, k]
        :return: None
        """
        # We make sure that the wavelengths are increasing, revering the arrays otherwise.
        if layer[1][0] > layer[1][-1]:
            layer[1] = layer[1][::-1]
            layer[2] = layer[2][::-1]
            layer[3] = layer[3][::-1]

        self.widths.append(layer[0])

        if len(layer) >= 5:
            self.models.append(layer[4])
            c = solcore.si(layer[5][0], 'nm')
            w = solcore.si(layer[5][1], 'nm')
            d = layer[5][2]  # = 0 for increasing, =1 for decreasing

            def mix(x):

                out = 1 + np.exp(-(x - c) / w)
                out = d + (-1) ** d * 1 / out

                return out

            n_data = lambda x: self.models[-1].n_and_k(x) * mix(x) + (1 - mix(x)) * interp1d(
                x=solcore.si(layer[1], 'nm'), y=layer[2], fill_value=layer[2][-1])(x)
            k_data = lambda x: interp1d(x=solcore.si(layer[1], 'nm'), y=layer[3], fill_value=layer[3][-1])(x)

            self.n_data.append(n_data)
            self.k_data.append(k_data)

        else:
            self.models.append([])
            self.n_data.append(interp1d(x=solcore.si(layer[1], 'nm'), y=layer[2], fill_value=layer[2][-1]))
            self.k_data.append(interp1d(x=solcore.si(layer[1], 'nm'), y=layer[3], fill_value=layer[3][-1]))


def calculate_rat(structure, wavelength, angle=0, pol='u', coherent=True, coherency_list=None, no_back_reflexion=True):
    """ Calculates the reflected, absorbed and transmitted intensity of the structure for the wavelengths and angles
    defined.

    :param structure: A solcore Structure object with layers and materials or a OptiStack object.
    :param wavelength: Wavelengths (in nm) in which calculate the data. An array.
    :param angle: Angle (in degrees) of the incident light. Default: 0 (normal incidence).
    :param pol: Polarisation of the light: 's', 'p' or 'u'. Default: 'u' (unpolarised).
    :param coherent: If the light is coeherent or not. If not, a coherency list must be added.
    :param coherency_list: A list indicating in which layers light should be treated as coeherent ('c') and in which
    incoherent ('i'). It needs as many elements as layers in the structure.
    :param no_back_reflexion: If reflexion from the back must be supressed. Default=True.
    :return: A dictionary with the R, A and T at the specified wavelengths and angle.
    """
    num_wl = len(wavelength)

    if 'OptiStack' in str(type(structure)):
        stack = structure
        stack.no_back_reflexion = no_back_reflexion
    else:
        stack = OptiStack(structure, no_back_reflexion=no_back_reflexion)

    if not coherent:
        if coherency_list is not None:
            assert len(coherency_list) == len(stack.get_widths()), \
                """Error: The coherency list must have as many elements (now {}) as the
                number of layers (now {}).""".format(
                        len(coherency_list)-len(stack.get_widths())+stack.num_layers,
                        stack.num_layers)

        else:
            raise Exception('Error: For incoherent or partly incoherent calculations you must supply the '
                            'coherency_list parameter with as many elements as the number of layers in the '
                            'structure')

    output = {'R': np.zeros(num_wl), 'A': np.zeros(num_wl), 'T': np.zeros(num_wl), 'all_p': [], 'all_s': []}

    if pol in 'sp':
        if coherent:
            out = tmm.coh_tmm(pol,
                              stack.get_indices(wavelength),
                              stack.get_widths(),
                              angle * degree,
                              wavelength)
            
            output['R'] = out['R']
            output['A'] = 1 - out['R'] - out['T']
            output['T'] = out['T']
        else:
            out = tmm.inc_tmm(pol,
                              stack.get_indices(wavelength),
                              stack.get_widths(),
                              coherency_list,
                              angle * degree,
                              wavelength)
            
            output['R'][i] = out['R']
            output['A'][i] = 1 - out['R'] - out['T']
            output['T'][i] = out['T']

    else:
        if coherent:
            out = tmm.unpolarized_RT(stack.get_indices(wavelength),
                                     stack.get_widths(),
                                     angle * degree,
                                     wavelength)
            output['R'] = out['R']
            output['A'] = 1 - out['R'] - out['T']
            output['T'] = out['T']
        else:
            out_p = tmm.inc_tmm('p',
                                stack.get_indices(wavelength),
                                stack.get_widths(),
                                coherency_list,
                                angle * degree,
                                wavelength)
            
            out_s = tmm.inc_tmm('s',
                                stack.get_indices(wavelength),
                                stack.get_widths(),
                                coherency_list,
                                angle * degree,
                                wavelength)

            output['R'] = 0.5 * (out_p['R'] + out_s['R'])
            output['T'] = 0.5 * (out_p['T'] + out_s['T'])
            output['A'] = 1 - output['R'] - output['T']
            output['all_p'] = out_p['power_entering_list']
            output['all_s'] = out_s['power_entering_list']

    return output


def calculate_ellipsometry(structure, wavelength, angle, no_back_reflexion=True):
    """ Calculates the ellipsometric parameters psi and delta. It can only deal with coherent light and the whole stack
    (including back surface) is considered, so caution must be taken when comparing the simulated results with
    experiments where the back surface is rough or layers are thick and coherent light propagation makes no sense.

    The optional argument no_back_reflexion can be included to add an extra layer on the back absorbing all light that
    reaches that position without any reflexion, to remove the reflexion from the back surface.

    :param structure: A solcore structure with layers and materials.
    :param wavelength: Wavelengths (in nm) in which calculate the data. An array.
    :param angle: A tupple or list with the angles (in degrees) in which to calculate the data.
    :param no_back_reflexion: If reflexion from the back must be supressed. Default=True.
    :return: A dictionary with psi and delta at the specified wavelengths and angles (2D arrays).
    """

    num_wl = len(wavelength)
    num_ang = len(angle)

    if 'OptiStack' in str(type(structure)):
        stack = structure
        stack.no_back_reflexion = no_back_reflexion
    else:
        stack = OptiStack(structure, no_back_reflexion=no_back_reflexion)

    output = {'psi': np.zeros((num_wl, num_ang)), 'Delta': np.zeros((num_wl, num_ang))}

    for i, ang in enumerate(angle):
        out = tmm.ellips(stack.get_indices(wavelength), stack.get_widths(), ang * degree, wavelength)
        output['psi'][:, i] = out['psi'] / degree

        # We revere the sign of Delta in order to use Woollam sign convention
        output['Delta'][:, i] = -out['Delta'] / degree

        output['Delta'][:, i] = np.where(output['Delta'][:, i] > 0, -output['Delta'][:, i], output['Delta'][:, i])
        output['Delta'][:, i] = np.where(output['Delta'][:, i] < 180, 180 - output['Delta'][:, i],
                                         output['Delta'][:, i])

    return output


# def calculate_absorption_profile(structure, wavelength, z_limit=None, steps_size=2, dist=None,
#                                  no_back_reflexion=True):
#     """ It calculates the absorbed energy density within the material. From the documentation:
#
#     'In principle this has units of [power]/[volume], but we can express it as a multiple of incoming light power
#     density on the material, which has units [power]/[area], so that absorbed energy density has units of 1/[length].'
#
#     Integrating this absorption profile in the whole stack gives the same result that the absorption obtained with
#     calculate_rat as long as the spacial mesh (controlled by steps_thinest_layer) is fine enough. If the structure is
#     very thick and the mesh not thin enough, the calculation might diverege at short wavelengths.
#
#     For now, it only works for normal incident, coherent light.
#
#     :param structure: A solcore structure with layers and materials.
#     :param wavelength: Wavelengths in which calculate the data (in nm). An array-like object.
#     :param z_limit: Maximum value in the z direction
#     :return: A dictionary containing the positions (in nm) and a 2D array with the absorption in the structure as a
#     function of the position and the wavelength.
#     """
#
#     num_wl = len(wavelength)
#
#     if 'OptiStack' in str(type(structure)):
#         stack = structure
#         stack.no_back_reflexion = no_back_reflexion
#     else:
#         stack = OptiStack(structure, no_back_reflexion=no_back_reflexion)
#
#     if dist is None:
#         if z_limit is None:
#             z_limit = np.sum(np.array(stack.widths))
#         dist = np.arange(0, z_limit, steps_size)
#
#     output = {'position': dist, 'absorption': np.zeros((num_wl, len(dist)))}
#
#     for i, wl in enumerate(wavelength):
#         out = tmm.coh_tmm('p', stack.get_indices(wl), stack.get_widths(), 0, wl)
#         for j, d in enumerate(dist):
#             layer, d_in_layer = tmm.find_in_structure_with_inf(stack.get_widths(), d)
#             data = tmm.position_resolved(layer, d_in_layer, out)
#             output['absorption'][i, j] = data['absor']
#
#     return output


def calculate_absorption_profile(structure, wavelength, z_limit=None, steps_size=2, dist=None,
                                   no_back_reflexion=True, angle=0):
    """ It calculates the absorbed energy density within the material. From the documentation:

    'In principle this has units of [power]/[volume], but we can express it as a multiple of incoming light power
    density on the material, which has units [power]/[area], so that absorbed energy density has units of 1/[length].'

    Integrating this absorption profile in the whole stack gives the same result that the absorption obtained with
    calculate_rat as long as the spacial mesh (controlled by steps_thinest_layer) is fine enough. If the structure is
    very thick and the mesh not thin enough, the calculation might diverege at short wavelengths.

    For now, it only works for normal incident, coherent light.

    :param structure: A solcore structure with layers and materials.
    :param wavelength: Wavelengths in which calculate the data (in nm). An array
    :param z_limit: Maximum value in the z direction
    :return: A dictionary containing the positions (in nm) and a 2D array with the absorption in the structure as a
    function of the position and the wavelength.
    """

    num_wl = len(wavelength)

    if 'OptiStack' in str(type(structure)):
        stack = structure
        stack.no_back_reflexion = no_back_reflexion
    else:
        stack = OptiStack(structure, no_back_reflexion=no_back_reflexion)

    if dist is None:
        if z_limit is None:
            z_limit = np.sum(np.array(stack.widths))
        dist = np.arange(0, z_limit, steps_size)

    output = {'position': dist, 'absorption': np.zeros((num_wl, len(dist)))}

    # print(stack.get_indices(wavelength).shape)
    out1 = tmm.coh_tmm('s', stack.get_indices(wavelength), stack.get_widths(), angle*degree, wavelength)
    layer, d_in_layer = tmm.find_in_structure_with_inf(stack.get_widths(), dist)
    data = tmm.position_resolved(layer, d_in_layer, out1)
    output['absorption'] = data['absor']

    return output

def calculate_inc_absorption_profile(structure, wavelength, coherency_list, angle=0, 
                                   no_back_reflexion=True):
    """ It calculates the absorbed energy density within the material for a partly incoherent stack. From the documentation:

    'In principle this has units of [power]/[volume], but we can express it as a multiple of incoming light power
    density on the material, which has units [power]/[area], so that absorbed energy density has units of 1/[length].'

    The function returns the total amount of light absorbed in each layer and a function for the position dependend absorption.
    Incoherent layers are treated as beer-lambert absorber.

    For now, it only works for unpolorized light.

    :param structure: A solcore structure with layers and materials.
    :param wavelength: Wavelengths in which calculate the data (in nm). An array
    :param coherency_list: A list indicating in which layers light should be treated as coeherent ('c') and in which
    incoherent ('i'). It needs as many elements as layers in the structure.
    :angle Angle (in degrees) of the incident light. Default: 0 (normal incidence).
    :param z_limit: Maximum value in the z direction
    :return: A tuple with the list of total absoprtion in each layer and a list of functions for the position resolved
    absorption.
    """

    if 'OptiStack' in str(type(structure)):
        stack = structure
        stack.no_back_reflexion = no_back_reflexion
    else:
        stack = OptiStack(structure, no_back_reflexion=no_back_reflexion)
        
    absorption_per_layer = []
    absorption_position_resolved = []
    stack_widths = stack.get_widths()
    indices = stack.get_indices(wavelength)

    out_p = tmm.inc_tmm('p', indices, stack_widths, coherency_list, angle * degree, wavelength)
    out_s = tmm.inc_tmm('s', indices, stack_widths, coherency_list, angle * degree, wavelength)
    
    abs_p = tmm.inc_absorp_in_each_layer(out_p)
    abs_s = tmm.inc_absorp_in_each_layer(out_s)
       
    for i in range(1, stack.num_layers+1):
        absorption_per_layer.append((abs_p[i]+abs_s[i])/2)
        if coherency_list[i] == 'c':
            abs_resolved_p = tmm.inc_find_absorp_analytic_fn(i, out_p)
            abs_resolved_s = tmm.inc_find_absorp_analytic_fn(i, out_s)
        elif coherency_list[i] == 'i':
            inc_index = out_p['all_from_inc'].index(i)
            abs_resolved_p = tmm.beer_lambert_derivative(out_p['VW_list'][:,inc_index],
                                                     out_p['th_list'][i],
                                                     indices[i],
                                                     wavelength)
            abs_resolved_s = tmm.beer_lambert_derivative(out_s['VW_list'][:,inc_index],
                                                     out_s['th_list'][i],
                                                     indices[i],
                                                     wavelength)
        absorption_position_resolved.append(
                lambda z: (abs_resolved_p(z)+abs_resolved_s(z))/2)
            
    return np.stack(absorption_per_layer), absorption_position_resolved

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from solcore import material, si
    from solcore.structure import Layer, Structure

    GaAs = material('GaAs')(T=300)
    InGaAs = material('InGaAs')(T=300, In=0.1)

    my_structure = Structure([
        Layer(si(3000, 'nm'), material=InGaAs),
        Layer(si(30, 'um'), material=GaAs),

    ])

    wavelength = np.linspace(450, 1100, 300)

    out = calculate_rat(my_structure, wavelength, coherent=True, no_back_reflexion=False)
    # #
    # plt.plot(wavelength, out['R'], 'b', label='Reflexion')
    # plt.plot(wavelength, out['A'], 'r', label='Absorption')
    # plt.plot(wavelength, out['T'], 'g', label='Transmission')
    # plt.legend()
    # plt.show()

    angles = [60, 65, 70]
    out = calculate_ellipsometry(my_structure, wavelength, angle=angles)
    #
    plt.plot(wavelength, out['psi'][:, 0], 'b', label='psi')
    plt.plot(wavelength, out['Delta'][:, 0], 'r', label='Delta')
    for i in range(1, len(angles)):
        plt.plot(wavelength, out['psi'][:, i], 'b')
        plt.plot(wavelength, out['Delta'][:, i], 'r')
    #
    # plt.legend()

    # out = calculate_absorption_profile_2(my_structure, wavelength, z_limit=3000)
    # # print(tuple(out['absorption'][0]))
    # # plt.plot(out['position'], out['absorption'][0])
    # A = np.zeros_like(wavelength)
    # #
    # for i, absorption in enumerate(out['absorption'][:]):
    #     A[i] = np.trapz(absorption, out['position'])
    # #
    # plt.plot(wavelength, A, 'k', label='Integrated Abs')

    #
    # plt.contourf(out['position'], wavelength, out['absorption'], 200)
    plt.legend()
    plt.show()
