""" This module contains a collection of mathematical models to calculate the dielectric constant of a material. The
modelling is largely based in the equations and naming used by:

J. A. Woollam, Guide to using WVASE: Spectroscopic Ellipsometry Data Acquisition and Analysis Software. 2012, pp. 1–696.

"""
import numpy as np
from solcore.absorption_calculator.kramers_kronig import Epsi1
from solcore.interpolate import interp1d


class Poles:
    """ Basic oscillator model for fitting ellipsometry and RAT data. It a pole without broadening and therefore can
    only represent non-absorbing materials. Usually, the central energy must lie outside the spectral range of interest.

    .. math:: \\epsilon = \\frac {A_n E_n^2} {E_n^2 - E^2}

    with:

    - An = the amplitud of the oscillator (dimensionless)
    - En = the central energy of the pole (in eV)

    It is equivalent to Pol.2 in WVASE

    """
    name = 'poles'
    var = 2

    def __init__(self, A=0., Ec=5.):
        """ Constructor of the pole oscillator. By default, the amplitudes are set to zero meaning that there are
        no poles.

        :param An: Amplitud of the pole (dimensionless)
        :param En: Center energy of the pole (must lay outside the spectral range of interest (eV)).
        """
        self.A = A
        self.Ec = Ec

    def dielectric(self, x):
        """ Returns the dielectric function as modelled by this class at the given x input values (a single value or
        an array of values.

        :param x: Spectral position in nm
        :return: The dielectric constant at this spectral position
        """
        x_eV = 1240 / x  # We transform it into eV

        return self.A * self.Ec ** 2 / (self.Ec ** 2 - x_eV ** 2)

    def __repr__(self):
        out = 'Poles oscillator with:\n\t- An = {}\n\t- En = {} eV\n'.format(self.A, self.Ec)

        return out


class Lorentz:
    """ Classic Lorentz oscillator involving a central energy and a broadening responsible for the absorption.

    .. math:: \\epsilon = \\frac {A_n E_n^2} {E_n^2 - E^2 - i Br_n E}

    with:

    - An = the amplitud of the oscillator (dimensionless)
    - Brn = the broadening (eV)
    - En = the central energy (eV)

    It is equivalent to Lor.2 in WVASE.

    """
    name = 'lorentz'
    var = 3

    def __init__(self, An=0., En=2., Brn=0.1):
        """ Constructor of the Lorentz oscillator. By default, the amplitude and broadening are set to zero meaning
        that there are no oscillator.

        :param An: Amplitud of the oscillator (dimensionless)
        :param En: Center energy of the oscillator (in eV).
        :param Brn: Broadening (in eV)
        """

        self.An = An
        self.En = En
        self.Brn = Brn

    def dielectric(self, x):
        """ Returns the dielectric function as modelled by this class at the given x input values (a single value or
        an array of values.

        :param x: Spectral position in nm
        :return: The dielectric constant at this spectral position
        """
        x_eV = 1240 / x
        return self.An * self.En ** 2 / (self.En ** 2 - x_eV ** 2 - 1.j * self.Brn * x_eV)

    def __repr__(self):
        out = 'Lorentz oscillator with:\n\t- An = {}\n\t- En = {} eV\n\t- Brn = {} eV\n'.format(self.An, self.En,
                                                                                                self.Brn)

        return out


class Gauss:
    """ Gauss oscillator involving a central energy and a broadening responsible for the absorption. It is defined as:

    .. math:: \\epsilon = \\epsilon_1 + i \\epsilon_2
    .. math:: \\epsilon_2 = A \\exp^{-{\\left( \\frac{E - E_c}{\\sigma} \\right)}^2} - A \\exp^{-{\\left( \\frac{E + E_c}{\\sigma} \\right)}^2}
    .. math:: \\sigma = \\frac{Br}{2 \\sqrt{\\ln(2)}}

    Epsilon_1 is calculated using the Kramers-Kronig relationships

    with:

    - A = the amplitud of the oscillator (dimensionless)
    - Br = the broadening (eV)
    - Ec = the central energy (eV)

    It is equivalent to Gau.0 in WVASE.

    """
    name = 'gauss'
    var = 3

    def __init__(self, A=0., Ec=2., Br=0.1):
        """ Constructor of the Gauss oscillator. By default, the amplitude and broadening are set to zero meaning
        that there are no oscillator.

        :param A: Amplitud of the oscillator (dimensionless)
        :param Ec: Center energy of the oscillator (in eV).
        :param Br: Broadening (in eV)
        """

        self.A = A
        self.Ec = Ec
        self.Br = Br

    def dielectric(self, x):
        """ Returns the dielectric function as modelled by this class at the given x input values (a single value or
        an array of values.

        :param x: Spectral position in nm
        :return: The dielectric constant at this spectral position
        """
        x_eV = 1240 / x
        s = 0.5 * self.Br / np.sqrt(np.log(2))

        epsi2 = lambda xx: self.A * (np.exp(-((xx - self.Ec) / s) ** 2) - np.exp(-((xx + self.Ec) / s) ** 2))

        epsi1 = Epsi1(epsi2)

        e1 = epsi1(x_eV)
        e2 = epsi2(x_eV)

        return e1 + 1.j * e2

    def __repr__(self):
        out = 'Gauss oscillator with:\n\t- A = {}\n\t- Ec = {} eV\n\t- Br = {} eV\n'.format(self.A, self.Ec, self.Br)

        return out


class Drude:
    """ The classic Drude oscillator describes free carrier effects on the dielectric response. Its form is a Lorentz
    oscillator with zero center energy.

    .. math:: \\epsilon = - \\frac {A_n Br_n} {E^2 + i Br_n E}

    with:

    - An = the amplitud of the oscillator (eV)
    - Brn = the broadening (eV)

    It is equivalent to Drd.0 in WVASE.

    """
    name = 'drude'
    var = 2

    def __init__(self, An=0., Brn=0.1):
        """ Constructor of the Drude oscillator. By default, the amplitude is set to zero meaning that there are no
        oscillator.

        :param An: Amplitud of the oscillator (in eV)
        :param Brn: Broadening (in eV)
        """

        self.An = An
        self.Brn = Brn

    def dielectric(self, x):
        """ Returns the dielectric function as modelled by this class at the given x input values (a single value or
        an array of values.

        :param x: Spectral position in nm
        :return: The dielectric constant at this spectral position
        """
        x_eV = 1240 / x
        epsi = - self.An * self.Brn / (x_eV ** 2 + 1.j * self.Brn * x_eV)

        return epsi

    def __repr__(self):
        out = 'Drude oscillator with:\n\t- An = {} eV\n\t- Brn = {} eV\n'.format(self.An, self.Brn)

        return out


class Cauchy:
    """ Cauchy oscillator model for fitting ellipsometry data. It is generally used to fit non absorbing materials,
    although an Urbach tail is included to indicate the onset of absorption.

    .. math:: N = An + \\frac {Bn} {E^2} + \\frac {Cn} {E^4}
    .. math:: K = A_k \\exp(B_k (E-C_k))
    .. math:: \\epsilon = (N + iK)^2

    with:

    - An, Bn and Cn = the Cauchy coeficients (dimensionless, µm-2 and µm-4, respectively)
    - Ak = the amplitud of the Urbach tail (dimensionless)
    - Bk = the multiplicative factor of the Urbach tail (eV-1)
    - Ck = the energy offset of the Urbach tail (eV) - it is user selected, but redundant with Ak.
    with Ak.

    It is equivalent to Chy.0 in WVASE.

    """
    name = 'cauchy'
    var = 5

    def __init__(self, An=0., Bn=0., Cn=0., Ak=0., Bk=0., Ck=0):
        """ Constructor of the Cauchy oscillator.

        :param An: 1st Cauchy coefficient (dimensionless)
        :param Bn: 2nd Cauchy coefficient (µm-2)
        :param Cn: 3rd Cauchy coefficient (µm-4)
        :param Ak: Amplitude of the Urback tail (dimensionless)
        :param Bk: multiplicative factor of the Urbach tail (eV-1)
        :param Ck: energy offset (eV)
        """

        self.An = An
        self.Bn = Bn
        self.Cn = Cn
        self.Ak = Ak
        self.Bk = Bk
        self.Ck = Ck

    def dielectric(self, x):
        """ Returns the dielectric function as modelled by this class at the given x input values (a single value or
        an array of values.

        :param x: Spectral position in nm
        :return: The dielectric constant at this spectral position
        """
        x_eV = 1240 / x

        N = self.An + self.Bn / (x / 1000.) ** 2 + self.Cn * (x / 1000.) ** 4
        K = self.Ak * np.exp(self.Bk * (x_eV - self.Ck))

        return (N + 1.j * K) ** 2

    def __repr__(self):
        out = 'Cauchy oscillator with:\n\t- An = {}\n\t- Bn = {} µm-2\n\t- Cn = {} µm-4\n\t- Ak = {}\n\t' \
              '- Bk = {} eV-1\n\t- Ck = {} eV\n'.format(self.An, self.Bn, self.Cn, self.Ak, self.Bk, self.Ck)

        return out


class PolySegment:
    """ The PolySegment model aims to approximate the dielectric function with a finite number of linear segments. It is
    similar to a point-by-point fitting but likely using much less free variables. The segments are created for
    Epsilon_2 and Epsilon_1 is calculated using the Kramers-Kronig relationships.

    It allows for a maximum freedom for fitting data while reducing noise, but does not provide any physically
    relevant information.

    The input parameters are the energies of the endpoints of the segments and the corresponding Epsilon_2 values.

    """
    name = 'polysegment'

    def __init__(self, energy, e2):
        """ Constructor of the PolySegment oscillator. By default, the amplitude and broadening are set to zero meaning
        that there are no oscillator.

        :param energy: Energies (eV)
        :param e2: Epsilon_2 values at the corresponding energies (dimensionless)
        """

        self.energy = energy
        self.e2 = e2
        self.epsi2 = interp1d(x=energy, y=e2)
        self.var = len(e2)

    def dielectric(self, x):
        """ Returns the dielectric function as modelled by this class at the given x input values (a single value or
        an array of values.

        :param x: Spectral position in nm
        :return: The dielectric constant at this spectral position
        """
        x_eV = 1240 / x

        epsi2 = lambda xx: self.epsi2(x_eV)

        epsi1 = Epsi1(epsi2)

        e1 = epsi1(x_eV)
        e2 = epsi2(x_eV)

        return e1 + 1.j * e2

    def __repr__(self):
        out = 'PolySegment oscillator with endpoints:\n' \
              '\t- Emin = {}\n' \
              '\t- Emax = {} eV\n' \
              '\t- Epsilon2(Emin) = {}\n' \
              '\t- Epsilon2(Emax) = {}\n'.format(self.energy[0], self.energy[-1], self.e2[0], self.e2[-1])

        return out


# Oscillator class designed to build an optical constrant model structure for customised modelling...
class Oscillator:
    def __init__(self, oscillator_type, material_parameters=None, **kwargs):
        self.oscillator = oscillator_type
        self.material_parameters = material_parameters
        self.arguments = kwargs

        if self.material_parameters == None:

            self.Material = "Not Specified"
        else:

            self.Material = self.material_parameters["Material"]

    def __str__(self):
        return "< Oscillator Type :: " + self.oscillator + ", Material :: " + self.Material + " >"


class DielectricConstantModel(list):
    """ This class groups together several models for the dielectric constant applicable to a given layer. Based on
    all of them, it calculates the dielectric constant and the complex refractive index at any energy. It requires as
    input an offset, offen refered as epsilon_infinity.

    """

    def __init__(self, e_inf=1., oscillators=()):
        """ Creates the dielectric constant model. The default model is air (constant dielectric constant = 1 at all
          wavelengths.

        :param e_inf: The offset dielectric constant. Default=1
        :param oscillators: A list containing the oscillators. Default=()
        """
        list.__init__(self, oscillators)

        self.e_inf = e_inf

    def add_oscillator(self, name, **kwargs):
        """ Adds an oscillator to the dielectric constant model.

        :param name: The name of the oscillator
        :param kwargs: The arguments needed to create that oscillator.
        :return:
        """
        available_models = {'poles': Poles,
                            'lorentz': Lorentz,
                            'drude': Drude,
                            'cauchy': Cauchy,
                            'gauss': Gauss,
                            'polysegment': PolySegment}

        assert name in available_models.keys(), 'Error: unknown oscillator model. ' \
                                                'Valid names are: {}'.format(available_models.keys())

        self.append(available_models[name](**kwargs))

    def remove_oscillator(self, idx):
        """ Removes an oscillator from the model.

        :param idx: The number of the oscillator (starting in 1)
        :return: None"""

        assert 0 < idx <= len(self), 'Error: The available oscillators are:\n {}'.format(self.__repr__())

        print(self[idx - 1])
        a = input('Do you want to remove this oscillator from the model (Y/n)?')

        if a in 'Yy':
            self.pop(idx - 1)

    def dielectric_constants(self, x):
        """ Returns the complex dielectric constant in the form e1 + ie2

        :param x: The spectral possition in nm (it can be an array)
        :return: The complex dielectric constant.
        """
        # If the wavelength is very small... then we are in the SI units. Lets change it to nm.
        try:
            in_si = any(x < 0.001)
        except TypeError:
            in_si = x < 0.001

        if in_si:
            x = x * 1e9

        epsi = self.e_inf
        for model in self:
            epsi = epsi + model.dielectric(x)

        return epsi

    def n_and_k(self, x):
        """ Returns the complex refractive index in the form n + ik

        :param x: The spectral possition in nm (it can be an array)
        :return: The complex refractive index.
        """

        return np.sqrt(self.dielectric_constants(x))

    def __repr__(self):

        out = 'e_inf = {}\n'.format(self.e_inf)
        for i in range(len(self)):
            out += '{} - {}'.format(i + 1, self[i].__repr__())

        return out


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from solcore.absorption_calculator.transfer_matrix import calculate_rat

    E = 2 * np.logspace(2, 4, 200)

    a = Poles(A=0, Ec=2)
    loz = Lorentz(An=1, En=2, Brn=0.5)
    drud = Drude(An=24.317000, Brn=0.125740)
    gau1 = Gauss(A=2.8468, Ec=0.1299, Br=0.0111)
    gau2 = Gauss(A=7.192, Ec=0.058331, Br=0.01682)
    gau3 = Gauss(A=1.9737, Ec=0.13991, Br=0.02144)
    c = Cauchy()

    model = DielectricConstantModel(e_inf=2.9987, oscillators=[gau1, gau2, gau3])

    # out = calculate_rat([[300, model]], x)

    # plt.semilogx(x, out['R'], 'b', label='Reflexion')
    # plt.semilogx(x, out['A'], 'r', label='Absorption')
    # plt.semilogx(x, out['T'], 'g', label='Transmission')

    n = model.n_and_k(E)
    #
    plt.semilogx(E, np.real(n), 'b', label='n')
    plt.semilogx(E, np.imag(n), 'r', label='k')
    plt.legend()
    plt.show()
