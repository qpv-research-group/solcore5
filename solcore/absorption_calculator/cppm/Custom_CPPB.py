# Custom_CPPB module for building customisable optical constant models based on any number of local oscillator functions

import numpy as np
import os, sys
import re
from solcore.science_tracker import science_reference
import solcore.constants as const


class Custom_CPPB():
    """
    Customisable Critical-Point Parabolic-Band model (CPPB) using expressions from Sadao Adachi and Charles Kim.

    The Custom_CPPB() class is designed to work with solcore's structure module, where now you can define a structure
    containing any number of oscillator functions and calculate the corresponding optical constants.

    """

    def __init__(self):

        # Define filepath to the MaterialParameters file import...
        DIR, SCRIPT = os.path.split(__file__)

        self.__PARAMS_PATH = os.path.join(DIR, 'Custom_CPPB_MaterialParamaters.txt')

        # Opens MaterialParameters file and dumps all content into a list...
        with open(self.__PARAMS_PATH, mode="r") as File:
            self.__MatParams = File.read().splitlines()

        # Single line lists all callable functions in the Custom_CPPB class, with the exception of __init__...
        self.__METHODS = [func for func in dir(Custom_CPPB)
                          if callable(getattr(Custom_CPPB, func)) and not func.startswith("__")]

    def Material_Params(self, material, *parameters):
        """
        Custom_CPPB.Material_Params() :: loads the required material parameters from the database file and
            returns a dictionary of all material parameters or individual variable for the required material.

        :param material: string variable, must match the available materials in Mod_CPPB_MaterialParameters.txt
        :param parameters: variable length argument to specify individual parameters if required.
                            Must be entered as strings.

        :return: Material Parameters
        """
        # Initialise variables
        ID = None
        Found = False
        PARAMS = {}
        PARAMS.setdefault("Material", material)

        for Lines in self.__MatParams:
            # Using regular expressions package, matching the material name with the correct section of the file.
            # If the name is found, ID is switched to True.
            if re.fullmatch("Material :: " + material, Lines) != None:
                ID = True
                Found = True

            # if the $END$ identifier is found in a line the ID is switched to False
            if "$END$" in Lines:
                ID = False

            # Whilst ID is true the parameters are collected and stored in a dictionary variable.
            # The additional conditional statements tell the code to only store 'variable' lines (containing =) and to ignore
            # comment lines (containing #).
            if ID is True and "=" in Lines and "$" not in Lines:
                split_line = Lines.split(" = ")

                PARAMS[split_line[0]] = float(split_line[1])

        if Found is False:
            raise ValueError("Material not found...")

        # If specific parameter or parameters are required, check to see if it exists and if so return it.
        if parameters.__len__() == 0:

            return PARAMS

        else:

            PARAMS_LIST = {}

            for p in parameters:

                if p in PARAMS:
                    PARAMS_LIST[p] = PARAMS[p]

                else:
                    raise ValueError("Material_Params ERROR :: Material parameter %s not found..." % p)

            return PARAMS_LIST

    def Broad(self, Gamma, Alpha, E0, energy):
        """
        Custom_CPPB.Broad() :: defines the frequency dependent gaussian broadening function proposed by C. Kim.

        :param Gamma: Broadening parameter (eV).
        :param Alpha: Parameter describing the transition from pure Lorentzian (Alpha=0) to approximated Gaussian
                (Alpha = 0.2) lineshape broadneing.
        :param E0: Critical point centre energy (eV).
        :param energy: Energy array (eV).

        :return: Frequency dependent broadening parameter.
        """

        # Scientific reference for this work...
        science_reference("Charles Kim, 'Modelling the optical dielectric function of semiconductors",
                          "C. C. Kim et al, 'Modelling the optical dielectric function of semiconductors: Extension of "
                          "the critical-point parabolic band approximation', Physical Review B 45(20) 11749, 1992")

        return Gamma * np.exp(-1 * Alpha * ((energy - E0) / Gamma) ** 2)

    def E0andE0_d0(self, energy, material_parameters=None, **kwargs):
        """
        Custom_CPPB.E0andE0_d0() = Intraband transition region, E0 and E0+delta0 transitions at the 3D M0 Critical
            Point.

        :param energy: energy: Energy array (eV).
        :param material_parameters: Parameter set imported using Material_Parameters() method. Not required as long as
            keyword arguments are specified.
        :param kwargs: These take in individual parameters for the model. Keywords should take the following form;
                    Gamma_E0
                    Alpha_E0
                    E0
                    E0_d0
                    A

        :return: E0 and E0_d0 critical point contributions to the complex dielectric function.
        """

        # Scientific reference for this work...
        science_reference("Sadao Adachi, Physical Properties of III-V Semiconductor Compounds",
                          "Adachi, S., Physical Properties of III-V Semiconductor Compounds, John Wiley & Sons (1992)")

        # Conditional statement determining where the input parameters come from...
        if material_parameters is None and bool(kwargs) is False:

            raise ValueError("No material parameters specified...")

        elif material_parameters is not None and bool(kwargs) is False:

            Params = material_parameters

        elif material_parameters is not None and bool(kwargs) is True:

            Params = material_parameters

            for key in kwargs:

                try:

                    Params[key] = kwargs[key]
                except KeyError:
                    print("Invalid material parameter...")

        elif bool(kwargs) is True:

            Params = {}

            for key in kwargs:

                try:

                    Params[key] = kwargs[key]
                except KeyError:
                    print("Invalid material parameter...")

        # Frequency dependent broadening...
        Gamma = self.Broad(Params["Gamma_E0"], Params["Alpha_E0"], Params["E0"], energy)

        # Oscillator model based on the form of an 3D M0 type critical point...
        ChiO = (energy + 1j * Gamma) / Params["E0"]
        ChiSO = (energy + 1j * Gamma) / Params["E0_d0"]

        F_ChiO = (ChiO ** (-2)) * (2 - (1 + ChiO) ** 0.5 - (1 - ChiO) ** 0.5)
        F_ChiSO = (ChiSO ** (-2)) * (2 - (1 + ChiSO) ** 0.5 - (1 - ChiSO) ** 0.5)

        Eps = (Params["A"] * (Params["E0"] ** (-1.5))) * (F_ChiO + 0.5 * ((Params["E0"] / Params["E0_d0"]) ** 1.5) \
                                                          * F_ChiSO)

        # Additional line to address change in phase of the imaginary signal.
        return Eps.real + 1j * abs(Eps.imag)

    def E1andE1_d1(self, energy, material_parameters, **kwargs):
        """
        Custom_CPPB.E1andE1_d1() = function describing the E1 and E1 + delta1 intraband transitions at the 3D M1 critical
            point.

        :param energy: energy: Energy array (eV).
        :param material_parameters: Parameter set imported using Material_Parameters() method. Not required as long as
            keyword arguments are specified.
        :param kwargs: These take in individual parameters for the model. Keywords should take the following form;
                    Gamma_E1
                    Alpha_E1
                    E1
                    E1_d1
                    B1
                    B1s

        :return: E1 and E1_d1 critical point contributions to the complex dielectric function.
        """

        # Scientific reference for this work...
        science_reference("Sadao Adachi, Physical Properties of III-V Semiconductor Compounds",
                          "Adachi, S., Physical Properties of III-V Semiconductor Compounds, John Wiley & Sons (1992)")

        # Conditional statement determining where the input parameters come from...
        if material_parameters is None and bool(kwargs) is False:

            raise ValueError("No material parameters specified...")

        elif material_parameters is not None and bool(kwargs) is False:

            Params = material_parameters

        elif material_parameters is not None and bool(kwargs) is True:

            Params = material_parameters

            for key in kwargs:

                try:

                    Params[key] = kwargs[key]
                except KeyError:
                    print("Invalid material parameter...")

        elif bool(kwargs) is True:

            Params = {}

            for key in kwargs:

                try:

                    Params[key] = kwargs[key]
                except KeyError:
                    print("Invalid material parameter...")

        # Frequency dependent broadening parameter...
        Gamma = self.Broad(Params["Gamma_E1"], Params["Alpha_E1"], Params["E1"], energy)

        # Oscillator function based on the form of an 2D M0 type critical point...
        Chi_1D = (energy + 1j * Gamma) / Params["E1"]
        Chi_1SD = (energy + 1j * Gamma) / Params["E1_d1"]

        Eps = (-1 * Params["B1"] * Chi_1D ** (-2)) * np.log(1 - Chi_1D ** 2) + (-1 * Params["B1s"] * Chi_1SD ** (-2)) \
                                                                               * np.log(1 - Chi_1SD ** 2)

        # Additional line to address change in phase of the imaginary signal.
        return Eps.real + 1j * abs(Eps.imag)

    def E_ID(self, energy, material_parameters, **kwargs):
        """
        Custom_CPPB.Eg_ID() :: From Adachi's formalism, contributions to the complex dielectric function from the indirect
            band-gap transitions.

        :param energy: energy: Energy array (eV).
        :param material_parameters: Parameter set imported using Material_Parameters() method. Not required as long as
            keyword arguments are specified.
        :param kwargs: These take in individual parameters for the model. Keywords should take the following form;
                    Gamma_Eg_ID
                    Alpha_Eg_ID
                    Eg_ED
                    Ec
                    D

        :return: Indirect band gap contributions to the complex dielectric function.
        """

        # Scientific reference for this work...
        science_reference("Sadao Adachi, Physical Properties of III-V Semiconductor Compounds",
                          "Adachi, S., Physical Properties of III-V Semiconductor Compounds, John Wiley & Sons (1992)")

        # Conditional statement determining where the input parameters come from...
        if material_parameters is None and bool(kwargs) is False:

            raise ValueError("No material parameters specified...")

        elif material_parameters is not None and bool(kwargs) is False:

            Params = material_parameters

        elif material_parameters is not None and bool(kwargs) is True:

            Params = material_parameters

            for key in kwargs:

                try:

                    Params[key] = kwargs[key]
                except KeyError:
                    print("Invalid material parameter...")

        elif bool(kwargs) is True:

            Params = {}

            for key in kwargs:

                try:

                    Params[key] = kwargs[key]
                except KeyError:
                    print("Invalid material parameter...")

        # Frequency dependent broadening parameter...
        Gamma = self.Broad(Params["Gamma_Eg_ID"], Params["Alpha_Eg_ID"], Params["Eg_ID"], energy)

        # Function describing the contributions of indirect transitions...
        term1 = -1 * (((Params["Eg_ID"] ** 2) / (energy + 1j * Gamma) ** 2) * np.log(Params["Ec"] / Params["Eg_ID"]))
        term2 = 0.5 * ((1 + (Params["Eg_ID"] / (energy + 1j * Gamma))) ** 2) * \
                np.log((energy + 1j * Gamma + Params["Ec"]) / (energy + 1j * Gamma + Params["Eg_ID"]))
        term3 = 0.5 * ((1 - (Params["Eg_ID"] / (energy + 1j * Gamma))) ** 2) * \
                np.log((energy + 1j * Gamma - Params["Ec"]) / (energy + 1j * Gamma - Params["Eg_ID"]))

        Eps = (2 * Params["D"] / np.pi) * (term1 + term2 + term3)

        # Additional line to address change in phase of the imaginary signal.
        return Eps.real + 1j * abs(Eps.imag)

    def E0_Exciton(self, energy, material_parameters, **kwargs):
        """
        Custom_CPPB.E0_Exciton() :: From Adachi's formalism, contributions to the complex dielectric function from bound
            excitons in the vicinity of the E0 fundamental band-gap.

        :param energy: energy: Energy array (eV).
        :param material_parameters: Parameter set imported using Material_Parameters() method. Not required as long as
            keyword arguments are specified.
        :param kwargs: These take in individual parameters for the model. Keywords should take the following form;
                    Gamma_Ex0
                    Alpha_Ex0
                    A_Ex
                    G_3D
                    n

        :return: Excitonic contributions at E0 to the complex dielectric function.
        """

        # Scientific reference for this work...
        science_reference("Sadao Adachi, Physical Properties of III-V Semiconductor Compounds",
                          "Adachi, S., Physical Properties of III-V Semiconductor Compounds, John Wiley & Sons (1992)")

        # Conditional statement determining where the input parameters come from...
        if material_parameters is None and bool(kwargs) is False:

            raise ValueError("No material parameters specified...")

        elif material_parameters is not None and bool(kwargs) is False:

            Params = material_parameters

        elif material_parameters is not None and bool(kwargs) is True:

            Params = material_parameters

            for key in kwargs:

                try:

                    Params[key] = kwargs[key]
                except KeyError:
                    print("Invalid material parameter...")

        elif bool(kwargs) is True:

            Params = {}

            for key in kwargs:

                try:

                    Params[key] = kwargs[key]
                except KeyError:
                    print("Invalid material parameter...")

        # Frequency dependent broadening parameter...
        Gamma = self.Broad(Params["Gamma_Ex0"], Params["Alpha_Ex0"], Params["E0"], energy)
        Eps_n = np.zeros((len(energy), int(Params["n"] + 1)), dtype="complex")

        for j in range(1, int(Params["n"] + 1)):
            Eps_n[:, j] = (Params["A_Ex"] / j ** 3) * ((Params["E0"] - (Params["G_3D"] / \
                                                                        (j ** 2)) - energy - 1j * Gamma) ** -1)

        Eps = np.sum(Eps_n, 1)

        # Additional line to address change in phase of the imaginary signal.
        return Eps.real + 1j * abs(Eps.imag)

    def E2(self, energy, material_parameters, **kwargs):
        """
        Custom_CPPB.E2() :: S. Adachi finds that the high energy E2 critical point is well approximated using a simple
            damped harmonic oscillator.

        :param energy: energy: Energy array (eV).
        :param material_parameters: Parameter set imported using Material_Parameters() method. Not required as long as
            keyword arguments are specified.
        :param kwargs: These take in individual parameters for the model. Keywords should take the following form;
                    Gamma_E2
                    Alpha_E2
                    E2
                    C

        :return: E2 contributions to the complex dielectric function.
        """

        # Scientific reference for this work...
        science_reference("Sadao Adachi, Physical Properties of III-V Semiconductor Compounds",
                          "Adachi, S., Physical Properties of III-V Semiconductor Compounds, John Wiley & Sons (1992)")

        # Conditional statement determining where the input parameters come from...
        if material_parameters is None and bool(kwargs) is False:

            raise ValueError("No material parameters specified...")

        elif material_parameters is not None and bool(kwargs) is False:

            Params = material_parameters

        elif material_parameters is not None and bool(kwargs) is True:

            Params = material_parameters

            for key in kwargs:

                try:

                    Params[key] = kwargs[key]
                except KeyError:
                    print("Invalid material parameter...")

        elif bool(kwargs) is True:

            Params = {}

            for key in kwargs:

                try:

                    Params[key] = kwargs[key]
                except KeyError:
                    print("Invalid material parameter...")

        # Frequency dependent broadning parameter...
        Gamma = self.Broad(Params["Gamma_E2"], Params["Alpha_E2"], Params["E2"], energy)

        # Damped harmonic oscillator function described by Adachi...
        Eps = Params["C"] / ((1 - (energy / Params["E2"]) ** 2) - 1j * (energy / Params["E2"]) * Gamma)

        # Additional line to address change in phase of the imaginary signal.
        return Eps.real + 1j * abs(Eps.imag)

    def Lorentz(self, energy, **kwargs):
        """
        Custom_CPPB.Lorentz() :: Classic Lorentz oscillator expression, taken from the J. A. Woollam ellipsometer
            manual.

        :param energy: Energy array (eV)
        :param kwargs: These take in individual parameters for the model. Keywords should take the following form;
                    Gamma
                    Alpha
                    E0
                    Amp

        :return: Lorentz oscillator contributions to the complex dielectric function.
        """

        # Scientific reference for this work...
        science_reference("J. A. Woollam, Guide to using WVASE", "J. A. Woollam Co. Inc., Copyright 1994-2012")

        # Remove material_parameters argument from the kwargs dictionary to avoid confusion...
        if "material_parameters" in kwargs:
            del kwargs["material_parameters"]

        Gamma = self.Broad(kwargs["Gamma"], kwargs["Alpha"], kwargs["E0"], energy)

        # Eps = (Params["C"] * Params["E0"]**2) / ((Params["E0"]**2 - energy**2) - 1j*energy*Gamma)
        Eps = (kwargs["Amp"] * kwargs["E0"]) / (kwargs["E0"] ** 2 - energy ** 2 - 1j * Gamma * energy)

        # Additional line to address change in phase of the imaginary signal.
        return Eps.real + 1j * abs(Eps.imag)

    def Sellmeier(self, energy, **kwargs):
        """
        Custom_CPPB.Sellmeier() :: Calculate the Sellmeier dispersion relation for the entirely real parts of the
            dielectric function. The Sellmeier relation is a sum of N terms.

        :param energy: Energy array (eV)
        :param kwargs: These take in individual parameters for the model. Parameters should be in units of "um".
            Keywords should take the following form;
                    An
                    Ln
                    NOTE: 'n' should take an integer from 1 to N where N is the number of Sellmeier terms.

        :return: Sellmeier contributions to the real part of the dielectric function.
        """

        # Scientific reference for this work...
        science_reference("Wilhelm Sellmeier's development of the Cauchy dispersion relation, taken from the J."
                          " A. Woollam WVASE manual", "J. A. Woollam Co. Inc., Copyright 1994-2012")

        # Remove material_parameters argument from the kwargs dictionary to avoid confusion...
        if "material_parameters" in kwargs:
            del kwargs["material_parameters"]

        # Work out number of Sellmeier terms from length of kwargs...
        if np.mod(len(kwargs), 2) != 0:

            raise ValueError("Insufficient number of Sellmeier coefficients :: No. kwargs == %d" % len(kwargs))

        else:

            terms = len(kwargs) / 2

        # Define empty array for individual Sellmeier terms...
        epsilon = np.zeros((int(terms), len(energy)), dtype=complex)

        # define the conversion constant to go between energy in eV and lambda in SI units...
        C = ((const.h * const.c) / const.q) * 1E6

        # calculate Sellmeier expression for the given coefficients...
        for i in range(0, int(terms)):
            epsilon[i] = (kwargs["A%d" % (i + 1)] * (C / energy) ** 2) / (
            (C / energy) ** 2 - kwargs["L%d" % (i + 1)] ** 2)

        return sum(epsilon, 0) + 1

    # def eps_calc(self, oscillator_structure, energy_array, components=False):
    #     """
    #     Custom_CPPB.eps_calc() :: Calculates the complex dielectric function of the presented oscillator structure.
    #
    #     :param oscillator_structure: Structure object containing information about the individual component functions.
    #     :param energy_array: Energy array (eV)
    #     :param components: Default=False, selects whether components are output along with the final result.
    #
    #     :return: Complex dielectric function and components from each oscillator function.
    #     """
    #
    #     # Statement checking whether the components argument has been specified as boolean...
    #     if isinstance(components, bool) == False:
    #         raise ValueError("'components' variable is invalid... state as either 'True or False")
    #
    #     # Find the length of oscillator structure and energy array and initialise a complex array for epsilon...
    #     epsilon = np.zeros((len(oscillator_structure), len(energy_array)), dtype="complex")
    #
    #     # Calculate the contributions from each specified oscillator...
    #     ind = 0
    #     for Osc_Num in oscillator_structure:
    #
    #         if Osc_Num.oscillator == "E0andE0_d0":
    #
    #             epsilon[ind] = self.E0andE0_d0(energy=energy_array, material_parameters=Osc_Num.material_parameters,
    #                                              **Osc_Num.arguments)
    #         elif Osc_Num.oscillator == "E1andE1_d1":
    #
    #             epsilon[ind] = self.E1andE1_d1(energy=energy_array, material_parameters=Osc_Num.material_parameters,
    #                                               **Osc_Num.arguments)
    #         elif Osc_Num.oscillator == "E_ID":
    #
    #             epsilon[ind] = self.E_ID(energy=energy_array, material_parameters=Osc_Num.material_parameters,
    #                                               **Osc_Num.arguments)
    #         elif Osc_Num.oscillator == "E0_Exciton":
    #
    #             epsilon[ind] = self.E_ID(energy=energy_array, material_parameters=Osc_Num.material_parameters,
    #                                         **Osc_Num.arguments)
    #         elif Osc_Num.oscillator == "E2":
    #
    #             epsilon[ind] = self.E2(energy=energy_array, material_parameters=Osc_Num.material_parameters,
    #                                      **Osc_Num.arguments)
    #         elif Osc_Num.oscillator == "Lorentz":
    #
    #             epsilon[ind] = self.Lorentz(energy=energy_array, **Osc_Num.arguments)
    #         else:
    #
    #             raise ValueError("Custom_CPPB() does not contain an oscillator function with name = " +
    #                              Osc_Num.oscillator)
    #
    #         ind += 1
    #
    #     # Sum all the individual contributions to produce the total complex dalectric function. The eps_inf parameter
    #     # is also added to the data.
    #     try:
    #         eps1_inf = oscillator_structure[0].material_parameters["eps1_inf"]
    #     except:
    #         eps1_inf = 0.0
    #
    #     epsilon_tot = sum(epsilon, 0) + eps1_inf
    #
    #     # Select whether component outputs are required...
    #     if components == False:
    #
    #         return epsilon_tot
    #
    #     elif components == True:
    #
    #         return (epsilon_tot, epsilon)

    def eps_calc(self, oscillator_structure, energy_array):
        """
        Custom_CPPB.eps_calc() :: Calculates the complex dielectric function of the presented oscillator structure.

        :param oscillator_structure: Structure object containing information about the individual component functions.
        :param energy_array: Energy array (eV)

        :return: A dictionary containing the complex dielectric function and components from each oscillator function.
            "eps" :: Sum of all contributions.
            "components" :: Individual components in a numpy array
        """

        # Find the length of oscillator structure and energy array and initialise a complex array for epsilon...
        epsilon = np.zeros((len(oscillator_structure), len(energy_array)), dtype="complex")

        # Calculate the contributions from each specified oscillator...
        i = 0

        for Osc in oscillator_structure:

            try:

                getattr(self, Osc.oscillator)
            except KeyError:

                print("Custom_CPPB() does not contain an oscillator function with name = " + Osc.oscillator)

            method = getattr(self, Osc.oscillator)

            epsilon[i] = method(energy=energy_array, material_parameters=Osc.material_parameters,
                                **Osc.arguments)

            i += 1

        # Sum all the individual contributions to produce the total complex dalectric function. The eps_inf parameter
        # is also added to the data.
        try:

            eps1_inf = oscillator_structure[0].material_parameters["eps1_inf"]
        except:

            eps1_inf = 0.0

        epsilon_tot = sum(epsilon, 0) + eps1_inf

        # A dictionary containing the total of all contributions to epsilon and also the components...
        return {"eps": epsilon_tot, "components": epsilon}

    def alpha_calc(self, oscillator_structure, energy_array):
        """
        Custom_CPPB.alpha_calc() :: Calculates the absorption coefficient from the complex dielectric function.

        :param oscillator_structure: Structure object containing information about the individual component functions.
        :param energy_array: Energy array (eV)

        :return: A dictionary containing absorption coeffcient and components from individual oscillator functions.
            "alpha" :: returns the total contributions to the absorption coefficient.
            "components" :: returns the individual components of alpha from each oscillator function.
        """

        # Define lambda functions for calculating k and alpha quickly...
        k_func = lambda Epsilon: np.sqrt(((Epsilon.real ** 2 + Epsilon.imag ** 2) ** 0.5 - Epsilon.real) / 2)
        alpha_func = lambda energy, k: ((4 * const.pi) / ((const.h * const.c) / (energy * const.q))) * k

        # Calculate the complex dielectruc function...
        epsilon = self.eps_calc(oscillator_structure, energy_array)

        # Calculate alpha for each component...
        Alpha_components = alpha_func(energy_array, k_func(epsilon["components"]))
        Alpha = alpha_func(energy_array, k_func(epsilon["eps"]))

        return {"alpha": Alpha, "components": Alpha_components}

    def nk_calc(self, oscillator_structure, energy_array):
        """
        Custom_CPPB.nk_calc() :: Calculates the refractive index and extinction coefficient from the complex dielectric
            function.

        :param oscillator_structure: Structure object containing information about the individual component functions.
        :param energy_array: Energy array (eV)

        :return: A dictionary containing refractive index, extinction coefficient and their components from individual
            oscillator functions.
            "n" :: returns refractive index.
            "n_components" :: returns components of refractive index from individual oscillator functions.
            "k" :: returns extinction coefficient.
            "k_components" :: returns components of extinction coefficient from individual oscillator functions.
        """

        # Define lambda function for calculating n...
        n_func = lambda Epsilon: np.sqrt(((Epsilon.real ** 2 + Epsilon.imag ** 2) ** 0.5 + Epsilon.real) / 2)
        # Define lambda function for calculating k...
        k_func = lambda Epsilon: np.sqrt(((Epsilon.real ** 2 + Epsilon.imag ** 2) ** 0.5 - Epsilon.real) / 2)

        # Calculate the complex dielectruc function...
        epsilon = self.eps_calc(oscillator_structure, energy_array)

        # Calculate n for each component...
        n_components = n_func(epsilon["components"])
        n = n_func(epsilon["eps"])

        # Calculate k for each component...
        k_components = k_func(epsilon["components"])
        k = k_func(epsilon["eps"])

        return {"n": n, "n_components": n_components, "k": k, "k_components": k_components}
