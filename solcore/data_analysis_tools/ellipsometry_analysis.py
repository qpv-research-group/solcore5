import numpy as np
import matplotlib.pyplot as plt
import codecs
from solcore.absorption_calculator import OptiStack
from solcore.absorption_calculator import calculate_ellipsometry
from scipy.optimize import minimize
from typing import Optional, List, Dict, Collection


class EllipsometryData:
    """ This class is used to open and arrange in a sensible way ellipsometry datafiles created with the WVASE
    software."""

    def __init__(self, path: str) -> None:
        """ Creates an ellipsometry data object and loads data into it the ellipsometry data from a file.

        :param path: Path to the file
        """

        self.comment: str = ''
        self.meas_setup: str = ''
        self.units: str = ''
        self.angles: List[float] = []

        self.data: Dict = {}
        self.dpol: Dict = {}
        self._load_data(path)

    def _load_data(self, path: str) -> None:
        """ This is the function that actually loads the data

        :param path: Path to the file
        :return: None
        """

        with codecs.open(path, "r", encoding='utf-8', errors='ignore') as f:
            lines = f.readlines()
            self.comment = lines[0]
            self.meas_setup = lines[1]

            for l in lines[2:]:
                data = l.strip().split()

                if data[0] in ['nm', '1/cm']:
                    self.units = data[0]
                    continue
                elif data[0][0] in '0123456789':
                    s = 0
                    key = float(data[1])
                elif data[0] in ['E', 'dpolE']:
                    s = 1
                    key = float(data[2])
                else:
                    continue

                if key not in self.data.keys():
                    self.angles.append(key)
                    self.data[key] = []
                    self.dpol[key] = []

                if data[0] == 'E' or data[0][0] in '0123456789':
                    self.data[key].append([])
                    self.data[key][-1].append(float(data[0 + s]))  # wl
                    self.data[key][-1].append(float(data[2 + s]))  # Psi
                    self.data[key][-1].append(float(data[4 + s]))  # err Psi
                    self.data[key][-1].append(float(data[3 + s]))  # Delta
                    self.data[key][-1].append(float(data[5 + s]))  # err Delta

                elif data[0] == 'dpolE':
                    self.dpol[key].append([])
                    self.dpol[key][-1].append(float(data[1]))  # wl
                    self.dpol[key][-1].append(float(data[3]))  # dpol
                    self.dpol[key][-1].append(float(data[4]))  # err dpol

        for key in self.data.keys():
            self.data[key] = np.array(self.data[key]).T

            if 'cm' in self.units:
                self.data[key][0] = 10000. / self.data[key][0]
            else:
                self.data[key][0] /= 1000

        self.units = 'um'

    def append_data(self, path: str) -> None:
        """ Ellipsometry data from more than one source can be combined in a single object. This function does that.

        :param path: Path to the file with data we want to add to the current object.
        :return:
        """

        new_data = EllipsometryData(path)

        for key in new_data.data.keys():
            if key not in self.data.keys():
                self.data[key] = new_data.data[key]
                self.dpol[key] = new_data.dpol[key]
            else:
                self.data[key] = np.hstack((new_data.data[key], self.data[key]))
                self.data[key] = self.data[key][:, np.argsort(self.data[key][0], axis=0)]

                try:
                    self.dpol[key] = np.vstack((new_data.dpol[key], self.dpol[key]))
                    self.dpol[key] = self.dpol[key][:, np.argsort(self.dpol[key][0], axis=0)]
                except:
                    pass

    def plot(self, log: bool = True):
        """ Plots the ellipsometry data contained in this object in a sensible way.

        :param log: True/False for ploting logarithmic X scale.
        :return: None
        """
        fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True)

        if log:
            for key in self.data.keys():
                ax[0].semilogx(self.data[key][0], self.data[key][1], 'o', mfc='none')
                ax[1].semilogx(self.data[key][0], self.data[key][3], 'o', mfc='none')
        else:
            for key in self.data.keys():
                ax[0].plot(self.data[key][0], self.data[key][1], 'o', mfc='none')
                ax[1].plot(self.data[key][0], self.data[key][3], 'o', mfc='none')

        plt.xlabel('Wavelength (µm)')
        ax[0].set_ylabel(r'$\Psi$ (º)')
        ax[1].set_ylabel(r'$\Delta$ (º)')
        plt.show()


class EllipsometryFitter:
    """ Class that contains all the tools to fit the properties of a OptiStack model to ellipsometry data."""

    def __init__(self, data: EllipsometryData, stack: Optional[OptiStack] = None,
                 vars: Optional[Collection] = None, wavelength_range: Optional[np.ndarray] = None) -> None:

        self.data = data
        self.stack = stack
        self.vars = vars
        self.range = wavelength_range

    def set_variables(self, vars: Collection) -> None:
        """ Sets the variables that will be used during the fitting. The input variable "vars" is a list with as many
        elements as layers in the model. Each of these elemenets is, in turn, a list with the names of the free
        variables of that layer. If the layers are defined by a dielectric model - rather than experimental n and k data
        the variable names other than the thickness and e_inf must be inside a list.

        vars =  [   [],
                    ['width', ['A1', 'A2']],
                    [['A1', 'E1'], ['A1', 'Br1']]
                ]

        The above variable list indicate that there are no free variables in the first layer, that the width and
        two amplitudes of the first model are free variables in the second layer, and for the first layer, that
        the free variables are the amplitude and energy of the first model, and the amplitude and broadening of the
        second one.

        The models are read in order. If there are free variables in the second model but not in the first one,
        then an empty list must be used. The last layer in the above example will look, in this case:

        [[], ['A1', 'Br1']

        :param vars: The list of free variables.
        :return: None.
        """
        assert len(vars) == len(self.stack.models), 'Error: the lenght of the vars list must be equal to the number ' \
                                                    'of layers in the stack = {}'.format(self.stack.num_layers)
        self.vars = vars

    def update_variables(self, v: List[float]) -> None:
        """ Updates the stack with the values of the variables contained in v. The meaning of each of the elements in v
        is given by the definitions in self.vars.

        :param v: An array with the new values of the variables.
        :return: None
        """
        k = 0
        for i, layer in enumerate(self.vars):
            # Loop over the layers
            if len(layer) > 0:
                model = 0
                for j, param in enumerate(layer):
                    # Loop over the parameterts of the layer
                    if param == 'width':
                        self.stack.widths[i] = v[k]
                        k += 1
                    # Parameters of the DielectricModel of the layer
                    elif param == 'e_inf':
                        self.stack.models[i].e_inf = v[k]
                        k += 1
                    elif param == 'e_gap':
                        self.stack.models[i].e_gap = v[k]
                        k += 1
                    else:
                        for val in param:
                            # Loop over the parameters of each oscillator within the DielectricModel
                            self.stack.models[i][model].__dict__[val] = v[k]
                            k += 1
                        model += 1

    def get_variables(self) -> List[float]:
        """ Provides a list with the current values of the free variables.

        :return: The values of the free variables, in the order they appear in self.vars.
        """
        out = []

        for i, layer in enumerate(self.vars):
            # Loop over the layers
            if len(layer) > 0:
                model = 0
                for j, param in enumerate(layer):
                    # Loop over the parameterts of the layer
                    if param == 'width':
                        out.append(self.stack.widths[i])
                    # Parameters of the DielectricModel of the layer
                    elif param == 'e_inf':
                        out.append(self.stack.models[i].e_inf)
                    elif param == 'e_gap':
                        out.append(self.stack.models[i].e_gap)
                    else:
                        for val in param:
                            # Loop over the parameters of each oscillator within the DielectricModel
                            out.append(self.stack.models[i][model].__dict__[val])
                        model += 1

        return out

    def set_range(self, wavelength_range: List[float]) -> None:
        """ Set the wavelengths in which to make the ellipsometry calculations. Ideally, it should be dense enought to
        capture all the features of the dielectric functions but sparse enought so the calculations are fast.

        :param wavelength_range: An array with the wavelength range in nanometers.
        :return: None
        """

        self.range = wavelength_range

    def mse(self, v: List[float] = None) -> np.ndarray:
        """ Calculates the mean-squared error between the experimental data and the modelled one. If the vars input
        is given, the variables are first updated and then the MSE is calculated.

        Since the experimental data might have a lot of wavelength points, we use a reduced number of them as defined
        in self.range in order to speed up the calculation of the ellipsometric parameters. The modelled data is then
        interpolated to the experimental one in order to calculate the MSE. While this will lead to higher uncertainties
        when the experimental data varies rapidly, it makes the calculation much faster.

        :param v: The free variables to update.
        :return: The mean-squared error between the experimental data and the modelled one
        """
        assert self.range is not None, 'Error: A wavelength range must be defined before calculating the MSE'

        if v is not None:
            self.update_variables(v)
            num_var = len(v)
        else:
            num_var = 0

        out = 0.
        elements = 0.

        mod_data = calculate_ellipsometry(self.stack, self.range, angle=self.data.angles)

        for i, ang in enumerate(self.data.angles):
            new_wl = self.data.data[ang][0] * 1000.
            elements += len(new_wl)
            psi = np.interp(new_wl, self.range, mod_data['psi'][:, i])
            delta = np.interp(new_wl, self.range, mod_data['Delta'][:, i])

            p = (self.data.data[ang][1] - psi) / self.data.data[ang][1]
            d = (self.data.data[ang][3] - delta) / self.data.data[ang][4]

            out += np.sum(p ** 2 + d ** 2)

        out = np.sqrt(1. / (2 * elements - num_var) * out)

        return np.array(out)

    def fit(self, show: bool = True) -> None:
        """ Calculate the best fit of the model to the data by minimising the mean-squared error using as free varaibles
         those defined in self.vars. The hessian matrix is used to calculate the 90% confident limits and the
         correlation matrix.

        :param show: If the result of the fit must be printed.
        :return:
        """

        guess = np.array(self.get_variables())

        def fun(x):
            return self.mse(x)

        self.fit_information = minimize(fun, guess, method='BFGS', options={'gtol': 1e-02, 'disp': True})
        self.calculate_uncertainties()
        self.update_variables(self.fit_information['x'])

        if show:
            print(self.fit_information)

    def calculate_uncertainties(self) -> None:
        """ Calculates the uncertainties in the calculated optimum parameters. For that, it uses the inverse of the
         hessian (or curvature matrix) C provided by the fitting algorithm. Two quantities are defined:

        Error for parameter i (= 90% confidence limit):
        .. math:: Err_{ii} = 1.65\\sqrt{MSE*C_{ii} }

        Correlation matrix:
        .. math:: CM_{ij} = \\frac{C_{ij}} {\\sqrt{C_{ii}} \\sqrt{C_{jj}}}

        The diagonals of the CM are always equal to one and the off-diagonal terms represent how correlated are the
        corresponding variables. Highly correlated variables (absolute values close to 1) mean that the minimum of the
        MSE is not very sharp with respect to changes between those two variables.

        The above information is stored in self.errors and self.correlations, respectively.

        :return: None
        """

        C = self.fit_information['hess_inv']
        MSE = self.fit_information['fun']

        self.errors = 1.65 * np.sqrt(MSE * np.diag(C))
        self.correlations = np.zeros_like(C)

        for i in range(len(C[0])):
            for j in range(len(C[0])):
                self.correlations[i, j] = C[i, j] / (np.sqrt(C[i, i]) * np.sqrt(C[j, j]))

    def plot(self, log: bool = True) -> None:
        """ Plots the expeirmental and the modelled data.

        :param log: If the X scale should be log.
        :return: None
        """

        mod_data = calculate_ellipsometry(self.stack, self.range, angle=self.data.angles)

        fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True)
        min_wl = 0
        max_wl = 1e6

        if log:
            for i, ang in enumerate(self.data.angles):
                new_wl = self.data.data[ang][0]
                if min(new_wl) > min_wl: min_wl = min(new_wl)
                if max(new_wl) < max_wl: max_wl = max(new_wl)

                psi = np.interp(new_wl, self.range / 1000., mod_data['psi'][:, i])
                delta = np.interp(new_wl, self.range / 1000., mod_data['Delta'][:, i])

                ax[0].semilogx(new_wl, self.data.data[ang][1], 'o', mfc='none')
                ax[1].semilogx(new_wl, self.data.data[ang][3], 'o', mfc='none')
                ax[0].semilogx(new_wl, psi, 'r')
                ax[1].semilogx(new_wl, delta, 'b')

        else:
            for i, ang in enumerate(self.data.angles):
                new_wl = self.data.data[ang][0] * 1000.
                if min(new_wl) > min_wl: min_wl = min(new_wl)
                if max(new_wl) < max_wl: max_wl = max(new_wl)

                psi = np.interp(new_wl, self.range / 1000., mod_data['psi'][:, i])
                delta = np.interp(new_wl, self.range / 1000., mod_data['Delta'][:, i])

                ax[0].plot(new_wl, self.data.data[ang][1], 'o', mfc='none')
                ax[1].plot(new_wl, self.data.data[ang][3], 'o', mfc='none')
                ax[0].semilogx(new_wl, psi, 'r')
                ax[1].semilogx(new_wl, delta, 'b')

        ax[0].set_xlim((min_wl, max_wl))
        ax[1].set_xlim((min_wl, max_wl))

        x = [0.3, 0.4, 0.5, 0.7, 1, 2, 3, 4, 5, 7, 10, 20]
        plt.xticks(x, x)
        plt.xlabel('Wavelength (µm)')
        ax[0].set_ylabel(r'$\Psi$ (º)')
        ax[1].set_ylabel(r'$\Delta$ (º)')
        plt.show()

    def __repr__(self) -> str:

        out = ''
        for i, layer in enumerate(self.stack.models):
            out += '\n LAYER {}:\n'.format(i)
            out += 'Width:\t {} nm\n'.format(self.stack.widths[i])
            out += layer.__repr__()
            out += '\n'

        return out


if __name__ == '__main__':
    from solcore.absorption_calculator.dielectric_constant_models import Poles, Lorentz, Drude, Gauss, Cauchy
    from solcore.absorption_calculator import DielectricConstantModel

    import matplotlib.pyplot as plt

    # Experimental data
    filename = '/Users/diego/Dropbox/WorkIC/Data/Imperial/RTA_samples/Ellipsometry_IR/nSi_bare.dat'
    filename2 = '/Users/diego/Dropbox/WorkIC/Data/Imperial/RTA_samples/Ellipsometry_VIS/160920/ellipvis_nsi2.dat'
    data = EllipsometryData(filename)
    data.append_data(filename2)
    # data.plot()

    # Models
    silicon_filename = '/Users/diego/Dropbox/WorkIC/Data/Materials_database/FTIR_mat/silicon.mat'
    si_data = np.loadtxt(silicon_filename, skiprows=3, unpack=True)
    si_data[0] = 1.24 / si_data[0]
    n = np.sqrt(si_data[1] + si_data[2] * 1.j)
    si_data[1] = np.real(n)
    si_data[2] = np.imag(n)

    # model_substrate = [300, si_data[0]*1000, si_data[1], si_data[2]]
    #

    drud = Drude(An=2.26598948, Brn=0.049)
    model_substrate = [300, si_data[0] * 1000, si_data[1], si_data[2],
                       DielectricConstantModel(e_inf=11.71255715, oscillators=[drud]), [2000, 500, 0]]

    siO2_filename = '/Users/diego/Dropbox/WorkIC/Data/Materials_database/FTIR_mat/SiO2_thermal.mat'
    sio2_data = np.loadtxt(siO2_filename, skiprows=3, unpack=True)

    sio2_data[0] = 1.24 / sio2_data[0]
    n = np.sqrt(sio2_data[1] + sio2_data[2] * 1.j)
    sio2_data[1] = np.real(n)
    sio2_data[2] = np.imag(n)
    model_sio2 = [1.5, sio2_data[0] * 1000, sio2_data[1], sio2_data[2]]

    # gau1 = Gauss(A=2.8468, Ec=0.1299, Br=0.0111)
    # gau2 = Gauss(A=7.192, Ec=0.058331, Br=0.01682)
    # gau3 = Gauss(A=1.9737, Ec=0.13991, Br=0.02144)
    # cau = Cauchy(An=1.5, Bn=0.002, Ak=1, Bk=0.2580645161)

    # model_sio2 = [1.5, DielectricConstantModel(e_inf=3, oscillators=[cau])]

    stack = OptiStack([model_sio2, model_substrate])

    # Fitting stuff
    fit = EllipsometryFitter(data, stack)
    variables = [['width'], ['e_inf', ['An']]]
    fit.set_variables(variables)
    wl = np.logspace(np.log10(300), np.log10(20000), 200)
    fit.set_range(wl)

    fit.plot()

    # fit.fit()
    #
    # print(fit.errors)
    # print(fit.correlations)

    # Ploting the result
    # fit.plot()

    # fit.get_all_vars()
    #
    # n_fin = fit.stack.models[0].n_and_k(wl)
    #
    # plt.semilogx(wl, np.real(n_ini))
    # plt.semilogx(wl, np.real(n_fin))
    # plt.show()
