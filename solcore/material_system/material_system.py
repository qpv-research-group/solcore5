import math  # hyperbolic functions etc in parameterisation
import os
import sys
from functools import lru_cache  # cache function calls to stop things taking forever / recalculating smae things

import numpy as np
from scipy.integrate import quad
import configparser

import solcore
from solcore.material_system import critical_point_interpolate
from solcore.parameter_system import ParameterSystem
from solcore.constants import h, c, q, kb, pi, electron_mass as m0, vacuum_permittivity
from solcore.singleton import Singleton
from solcore.absorption_calculator.sopra_db import sopra_database, compounds_info
from solcore.absorption_calculator.nk_db import nkdb_load_n, nkdb_load_k
from solcore.material_data import calculate_mobility


class MaterialSystem(metaclass=Singleton):
    """ The core class that manage the materials in solcore.
    """

    def __init__(self, sources=None):
        self.known_materials = {}
        self.sources = {}
        if sources is not None:
            self.sources = {k.lower(): sources(k) for k in sources()}

    def read(self, source, value):
        """ Reads a new source. """
        self.sources[source.lower()] = value

    def material(self, name, sopra=False, nk_db=False):
        """ This function checks if the requested material exists and creates a class
        that contains its properties, assuming that the material does not exists in
        the database, yet.

        Such class will serve as the base class for all the derived materials based on that SpecificMaterial.
        For example, N-type GaAs and P-type GaAs use the same SpecificMaterial, just with a different doping, and the
        happens with GaAs at 300K or 200K.

        The derived materials based on a SpecificMaterial are instances of the SpecificMaterial class.

        >>> GaAs = solcore.material('GaAs')      # The SpecificMaterial class
        >>> n_GaAs = GaAs(Nd=1e23)               # Instance of the class
        >>> p_GaAs = GaAs(Na=1e22)               # Another instance of GaAs with different doping

        >>> AlGaAs = solcore.material('AlGaAs')  # The SpecificMaterial class
        >>> AlGaAs_1 = AlGaAs(Al=0.3)            # Instance of the class. For compounds, the variable element MUST be present
        >>> AlGaAs_2 = AlGaAs(Al=0.7, T=290)     # Different composition and T (the default is T=300K)


        The material is created from the parameters in the parameter_system and the n and k data if available. If the
        n and k data does not exists - at all or for that composition - then n=1 and k=0 at all wavelengths. Keep in
        mind that the available n and k data is valid only at room temperature.

        :param name: Name of the material
        :param sopra: If a SOPRA material must be used, rather than the normal database material, in case both exist.
        :return: A class of that material
        """

        suffix = ''
        # First we check if the material exists. If not, some help is provided
        try:
            if sopra:
                sopra_database(Material=name)
                suffix = '_sopra'
            elif nk_db:
                suffix = '_nk'
            else:
                ParameterSystem().database.options(name)
        except configparser.NoSectionError:
            try:
                sopra_database(Material=name)
                sopra = True
                suffix = '_sopra'
            except:
                valid_materials = sorted(ParameterSystem().database.sections())
                valid_materials.remove('Immediate Calculables')
                valid_materials.remove('Final Calculables')
                print(
                    '\nMaterial ERROR: "{}" is not in the semiconductors database or in the SOPRA database. Valid semiconductor materials are: '.format(
                        name))
                for v in valid_materials:
                    if "x" in ParameterSystem().database.options(v):
                        x = ParameterSystem().database.get(v, 'x')
                        print('\t {}  \tx = {}'.format(v, x))
                    else:
                        print('\t {}'.format(v))
                print(
                    '\nIn compounds, check that the order of the elements is the correct one (eg. GaInSb is OK but InGaSb'
                    ' is not).')
                val = input('\nDo you want to see the list of available SOPRA materials (y/n)?')
                if val in 'Yy':
                    sopra_database.material_list()

                sys.exit()
        except solcore.absorption_calculator.sopra_db.SOPRAError:
            pass

        # Then we check if the material has already been created. If not, we create it.
        if name + suffix in self.known_materials:
            return self.known_materials[name + suffix]
        elif sopra:
            return self.sopra_material(name)
        elif nk_db:
            return self.nk_material(name)
        else:
            return self.parameterised_material(name)

    def parameterised_material(self, name):
        """ The function that actually creates the material class. """

        if "x" in ParameterSystem().database.options(name):
            self.composition = [ParameterSystem().database.get(name, "x")]
        else:
            self.composition = []

        class SpecificMaterial(BaseMaterial):
            material_string = name
            composition = self.composition

            def __init__(self, T=300, **kwargs):
                BaseMaterial.__init__(self, T=T, **kwargs)

        if name.lower() in self.sources.keys():
            SpecificMaterial.material_directory = self.sources[name.lower()]
            extension = ''
            if len(SpecificMaterial.composition) == 0:
                extension = '.txt'

            SpecificMaterial.n_path = os.path.join(SpecificMaterial.material_directory, 'n' + extension)
            SpecificMaterial.k_path = os.path.join(SpecificMaterial.material_directory, 'k' + extension)

        SpecificMaterial.__name__ = name

        self.known_materials[name] = SpecificMaterial
        return SpecificMaterial

    def sopra_material(self, name):
        """ Creates an optical material fromt he SOPRA database. """

        try:
            self.composition = [compounds_info.get(name, "x")]
        except:
            self.composition = []

        class SpecificMaterial(BaseMaterial):
            material_string = name
            composition = self.composition
            data = sopra_database(Material=name)

            def __init__(self, T=300, **kwargs):
                BaseMaterial.__init__(self, T=T, **kwargs)

            def __getattr__(self, attrname):  # only used for unknown attributes
                if attrname == "n":
                    return self.n_interpolated
                if attrname == "k":
                    return self.k_interpolated

                raise AttributeError(
                    'Parameter "{}" not available for this SOPRA material "{}".'.format(attrname, self.material_string))

            @lru_cache(maxsize=1)
            def load_n_data(self):
                if len(self.composition) == 0:
                    wl, n = self.data.load_n()
                    self.n_data = np.vstack((wl * 1e-9, n))

            @lru_cache(maxsize=1)
            def load_k_data(self):
                if len(self.composition) == 0:
                    wl, k = self.data.load_k()
                    self.k_data = np.vstack((wl * 1e-9, k))

            def n_interpolated(self, x):
                assert len(self.composition) <= 1, "Can't interpolate 2d spectra yet"

                try:
                    self.load_n_data()
                except:
                    print('Material "{}" does not have n-data defined. Returning "ones"'.format(self.material_string))
                    return np.ones_like(x)

                if len(self.composition) == 0:
                    y = np.interp(x, self.n_data[0], self.n_data[1])
                else:
                    wl, y, trash = self.data.load_composition(Lambda=x * 1e9,
                                                              **{self.composition[0]: self.main_fraction * 100})

                return y

            def k_interpolated(self, x):
                assert len(self.composition) <= 1, "Can't interpolate 2d spectra yet"

                try:
                    self.load_k_data()
                except:
                    print('Material "{}" does not have k-data defined. Returning "zeroes"'.format(self.material_string))
                    return np.zeros_like(x)

                if len(self.composition) == 0:
                    y = np.interp(x, self.k_data[0], self.k_data[1])
                else:
                    wl, trash, y = self.data.load_composition(Lambda=x * 1e9,
                                                              **{self.composition[0]: self.main_fraction * 100})

                return y

        SpecificMaterial.__name__ = name + '_sopra'

        self.known_materials[name + '_sopra'] = SpecificMaterial
        return SpecificMaterial

    def nk_material(self, name):
        """ Creates an optical material from the refractiveindex.info database. """

        class SpecificMaterial(BaseMaterial):

            def __init__(self, T=300, **kwargs):
                BaseMaterial.__init__(self, T=T, **kwargs)

            def __getattr__(self, attrname):  # only used for unknown attributes
                if attrname == "n":
                    return self.n_interpolated
                if attrname == "k":
                    return self.k_interpolated

                raise AttributeError(
                    'Parameter "{}" not available for this refractiveindex.info material "{}".'.format(attrname, self.material_string))

            @lru_cache(maxsize=1)
            def load_n_data(self):
                    wl, n = nkdb_load_n(name)
                    self.n_data = np.vstack((wl * 1e-6, n))  # wavelengths in DB are in microns

            @lru_cache(maxsize=1)
            def load_k_data(self):
                    wl, k = nkdb_load_k(name)
                    self.k_data = np.vstack((wl * 1e-6, k))

            def n_interpolated(self, x):
                assert len(self.composition) <= 1, "Can't interpolate 2d spectra yet"

                try:
                    self.load_n_data()
                except:
                    print('Material "{}" does not have n-data defined. Returning "ones"'.format(self.material_string))
                    return np.ones_like(x)

                y = np.interp(x, self.n_data[0], self.n_data[1])

                return y

            def k_interpolated(self, x):
                assert len(self.composition) <= 1, "Can't interpolate 2d spectra yet"

                try:
                    self.load_k_data()
                except:
                    print('Material "{}" does not have k-data defined. Returning "zeroes"'.format(self.material_string))
                    return np.zeros_like(x)

                y = np.interp(x, self.k_data[0], self.k_data[1])

                return y

        SpecificMaterial.__name__ = name + '_nk'

        self.known_materials[name + '_nk'] = SpecificMaterial
        return SpecificMaterial


class BaseMaterial:
    """ The solcore base material class
    """
    material_string = "Unnamed"
    composition = []
    main_fraction = 0
    material_directory = None
    k_path = None
    n_path = None
    strained = False

    def __init__(self, T, **kwargs):
        self.Na = kwargs["Na"] if "Na" in kwargs else 1
        self.Nd = kwargs["Nd"] if "Nd" in kwargs else 1
        self.__dict__.update(kwargs)
        self.T = T

        for element in self.composition:
            setattr(self, element, kwargs[element])

        if len(self.composition) != 0:
            self.main_fraction = kwargs[self.composition[0]]

        self.key_parameters = sorted(list(kwargs.keys()))

    def __getattr__(self, attrname):  # only used for unknown attributes.
        if attrname == "n":
            try:
                return self.n
            except:
                return self.n_interpolated
        if attrname == "k":
            try:
                return self.k
            except:
                return self.k_interpolated
        if attrname == "electron_affinity":
            try:
                return ParameterSystem().get_parameter(self.material_string, attrname)
            except:
                # from http://en.wikipedia.org/wiki/Anderson's_rule and GaAs values
                return (0.17 + 4.59) * q - self.valence_band_offset - self.band_gap
        if attrname == "electron_mobility":
            try:
                return ParameterSystem().get_parameter(self.material_string, attrname)
            except:
                return calculate_mobility(self.material_string, False, self.Nd, self.main_fraction)
        if attrname == "hole_mobility":
            try:
                return ParameterSystem().get_parameter(self.material_string, attrname)
            except:
                return calculate_mobility(self.material_string, True, self.Na, self.main_fraction)
        if attrname == "Nc":
            # return 2 * (2 * pi * self.eff_mass_electron_Gamma * m0 * kb * self.T / h ** 2) ** 1.5
            # Changed by Diego 22/10/18: some materials dont have _Gamma but all have the normal electron effective mass
            return 2 * (2 * pi * self.eff_mass_electron * m0 * kb * self.T / h ** 2) ** 1.5
        if attrname == "Nv":
            # Strictly speaking, this is valid only for III-V, zinc-blend semiconductors
            mhh = self.eff_mass_hh_z
            mlh = self.eff_mass_lh_z
            Nvhh = 2 * (2 * pi * mhh * m0 * kb * self.T / h ** 2) ** 1.5
            Nvlh = 2 * (2 * pi * mlh * m0 * kb * self.T / h ** 2) ** 1.5
            return Nvhh + Nvlh
        if attrname == "ni":
            return np.sqrt(self.Nc * self.Nv * np.exp(-self.band_gap / (kb * self.T)))
        if attrname == "radiative_recombination":
            inter = lambda E: self.n(E) ** 2 * self.alphaE(E) * np.exp(-E / (kb * self.T)) * E ** 2
            upper = self.band_gap + 10 * kb * self.T
            return 1.0 / self.ni ** 2 * 2 * pi / (h ** 3 * c ** 2) * quad(inter, 0, upper)[0]
        if attrname == "permittivity":
            return self.relative_permittivity * vacuum_permittivity

        kwargs = {element: getattr(self, element) for element in self.composition}
        kwargs["T"] = self.T
        return ParameterSystem().get_parameter(self.material_string, attrname, **kwargs)

    @lru_cache(maxsize=1)
    def load_n_data(self):
        if len(self.composition) == 0:
            self.n_data = np.loadtxt(self.n_path, unpack=True)
        else:
            self.n_data, self.n_critical_points = critical_point_interpolate.load_data_from_directory(self.n_path)

    @lru_cache(maxsize=1)
    def load_k_data(self):
        if len(self.composition) == 0:
            self.k_data = np.loadtxt(self.k_path, unpack=True)
        else:
            self.k_data, self.k_critical_points = critical_point_interpolate.load_data_from_directory(self.k_path)

    def n_interpolated(self, x):
        assert len(self.composition) <= 1, "Can't interpolate 2d spectra yet"

        try:
            self.load_n_data()
        except:
            print('Material "{}" does not have n-data defined. Returning "ones"'.format(self.material_string))
            return np.ones_like(x)

        if len(self.composition) == 0:
            y = np.interp(x, self.n_data[0], self.n_data[1])
        else:
            x, y, self.critical_points_n = critical_point_interpolate.critical_point_interpolate(
                self.n_data,
                self.n_critical_points,
                self.main_fraction,
                x
            )

        return y

    def k_interpolated(self, x):
        assert len(self.composition) <= 1, "Can't interpolate 2d spectra yet"

        try:
            self.load_k_data()
        except:
            print('Material "{}" does not have k-data defined. Returning "zeros"'.format(self.material_string))
            return np.zeros_like(x)

        if len(self.composition) == 0:
            y = np.interp(x, self.k_data[0], self.k_data[1])
        else:
            x, y, self.critical_points_k = critical_point_interpolate.critical_point_interpolate(
                self.k_data,
                self.k_critical_points,
                self.main_fraction,
                x
            )
        return y

    def get(self, parameter):
        return getattr(self, parameter)

    def alpha(self, wavelength):
        return 4 * math.pi * self.k(wavelength) / wavelength

    def alphaE(self, energy):
        return self.alpha(h * c / energy)

    def latex_string(self):
        s = self.material_string
        for element in self.composition:
            if self.__dict__[element] != 0:
                s = s.replace(element, "{}_{{{:.3f}}}").format(element, self.__dict__[element])
            else:
                s = s.replace(element, "")
        return "$\mathrm{{{}}}$".format(s)

    def plain_string(self):
        s = self.material_string
        for element in self.composition:
            if self.__dict__[element] != 0:
                s = s.replace(element, "{}{:.3f}").format(element, self.__dict__[element])
            else:
                s = s.replace(element, "")
        return "{}".format(s)

    def html_string(self):
        s = self.material_string
        for element in self.composition:
            if self.__dict__[element] != 0:
                s = s.replace(element, "{}<sub>{:.3f}</sub>").format(element, self.__dict__[element])
            else:
                s = s.replace(element, "")
        return "{}".format(s)

    def __repr__(self):
        parameters = ["{}={}".format(key, getattr(self, key)) for key in self.key_parameters]
        parameters_string = " ".join(parameters)
        return "<'{}' material{}>".format(
            self.material_string,
            "" if len(self.key_parameters) == 0 else " " + parameters_string)
