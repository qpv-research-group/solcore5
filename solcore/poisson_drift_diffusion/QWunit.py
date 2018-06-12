import numpy as np
from scipy.special import erfc
from scipy.interpolate import interp1d
from scipy.integrate import quad
import sys

from solcore import constants, si
from solcore.structure import Layer, Structure
from solcore.quantum_mechanics import schrodinger
from solcore.quantum_mechanics.kp_bulk import kp8x8_bulk
from solcore.quantum_mechanics.strain import strain_calculation_parameters
from solcore.absorption_calculator import adachi_alpha
from solcore.science_tracker import science_reference

q = constants.q
pi = constants.pi
h = constants.h
kb = constants.kb
m0 = constants.electron_mass
vacuum_permittivity = constants.vacuum_permittivity
c = constants.c


class QWunit(Structure):
    """ Assembles a group of layers as a quantum well, calculating its properties as a whole
    
    - It needs a minimum of 3 layers
    - There must be exactly one layer with the role = "well".
    - Top and bottom layers are assumed to be the barriers.
    - Anything not a barrier not a well will be considered as "interlayer".    
    
    Since the Poisson - Drift-Diffusion solver can not work with superlattices, there is no point of considering the solution of more than one QW.
    All QWs entering into the solver are independent, regardless of the width of the barriers.
    
    To account for the possible presence of nearby QWs, perdiodic boundary conditions can be used to solve the Schrodinger equation. 
    If barriers + interlayers are thick enought, this make no difference but if they are thin, it affects the number of confined levels. 
    It does not 'create' minibands, though, so the results with thin barriers must be used with caution.  
    
    """

    def __init__(self, *args, **kwargs):
        """ Initialises the QWunit, checking that it has all the necessary layers """
        super().__init__(*args, **kwargs)

        if len(self) < 3:
            print("ERROR creating a QWunit: a minimum of 3 layers is necessary to make a unit. Only {0} found. ".format(
                len(self)))
            sys.exit()

        self.QW_number = 0
        self.QW_width = 0

        self[0].__dict__["role"] = "barrier"
        self[-1].__dict__["role"] = "barrier"
        for index in range(1, len(self) - 1):

            if self[index].role in ["well", "qw", "QW", "Qw", "Well"]:
                self.QW_number += 1
                self.QW_width = self.QW_width + self[index].width
                self[index].role = "well"
            else:
                self[index].__dict__["role"] = "interlayer"

        # if self.QW_number != 1:
        #     print("ERROR creating a QWunit: There must be exactly one layer with the role = 'well'. {0} found. ".format(
        #         self.QW_number))
        #     sys.exit()

        if 'substrate' not in kwargs.keys():
            print("ERROR creating a QWunit: There must be a substrate defined with the keyword 'substrate'. ")
            sys.exit()

        # We update the structure of labels since that is what the quantum mechanics modules uses to identify the regions, not the role. 
        for index in range(len(self)):
            self.labels[index] = self[index].role

        self.T = kwargs['T'] if 'T' in kwargs.keys() else 298.0

    # ------
    def RecalculateKP(self, index, mode):
        """ Recalculates the band structure of the QW unit to account for modifications due to strain. 
        
        This is already done in the Schrodinger solver, but it does not produce that information as output. We do it so here.
        """
        if mode == 'kp8x8_bulk':
            kp_result = kp8x8_bulk(material=self[index].material, host_material=self.substrate, averaging_points=6)
            c, hh, lh, so, me, mhh, mlh, mso = kp_result

            # We recalculate the bandgap adn the electron afinity
            Ec = self[index].material.valence_band_offset + self[index].material.band_gap
            Eb = c - max(hh, lh)
            Xb = self[index].material.electron_affinity - c + Ec
        else:
            me = self[index].material.eff_mass_electron_Gamma * m0
            mhh = self[index].material.eff_mass_hh_z * m0
            mlh = self[index].material.eff_mass_lh_z * m0

            # Idem
            strain_parameters = strain_calculation_parameters(self.substrate, self[index].material, should_print=False)
            Eb = self[index].material.band_gap + strain_parameters.delta_Ec + max(
                strain_parameters.delta_Ehh, strain_parameters.delta_Elh)
            Xb = self[index].material.electron_affinity - strain_parameters.delta_Ec

        Ncc = 2 * (2 * pi * me * kb * self.T / h ** 2) ** 1.5
        Nvhh = 2 * (2 * pi * mhh * kb * self.T / h ** 2) ** 1.5
        Nvlh = 2 * (2 * pi * mlh * kb * self.T / h ** 2) ** 1.5

        return Eb, Xb, me, mhh, mlh, Ncc, Nvhh, Nvlh
        # ------

    def get_average_params(self, x, y):
        """ Calculates the average value of a parameter y in each layer of the structure.

        :param x: x values
        :param y: the parameter of interest at positions x
        :return: the average value of y in layer i
        """

        widths = np.array([0] + [layer.width for layer in self])
        cumul = np.cumsum(widths)

        aver = np.zeros(len(self))
        for i in range(len(self)):
            ind = np.searchsorted(x, cumul[i:i + 2])
            yy = y[ind[0]:ind[1] - 1]
            xx = x[ind[0]:ind[1] - 1]
            aver[i] = np.trapz(yy, xx) / (xx[-1] - xx[0])

        return aver

    def get_probability_per_layer(self, x, y):
        """ Calculates the average value of a parameter y in each layer of the structure.

        :param x: x values
        :param y: the parameter of interest at positions x
        :return: the average value of y in layer i
        """

        widths = np.array([0] + [layer.width for layer in self])
        cumul = np.cumsum(widths)

        aver = np.zeros(len(self))
        for i in range(len(self)):
            ind = np.searchsorted(x, cumul[i:i + 2])
            yy = y[ind[0]:ind[1] - 1]
            xx = x[ind[0]:ind[1] - 1]
            aver[i] = np.trapz(yy, xx)

        return aver

    def RecalculateBandEdges(self, mode, SR):
        """The effective bandgap and effective affinity of each layer depends on the energy level, e.g.: the minimum energy for electrons is not the band edge of the conduction band, but the ground level. In this function we modify that creating an effective electron afinity and band gaps for all layers in the QW.

        For the barriers, the electron afinity and band gap is the same than in bulk, modified by the kp calculation, if necesary."""

        eff_c_shift = np.where(SR["E"]["Ee"][0] > SR["potentials"]["Ve"], SR["E"]["Ee"][0],
                               SR["potentials"]["Ve"]) - SR["potentials"]["Ve"]

        max_VB = np.where(SR["potentials"]["Vhh"] > SR["potentials"]["Vlh"], SR["potentials"]["Vhh"],
                          SR["potentials"]["Vlh"])

        max_hole_level = max(SR["E"]["Ehh"][0], SR["E"]["Elh"][0])

        eff_h_shift = np.where(max_hole_level < max_VB, max_hole_level, max_VB) - max_VB

        shift_e = self.get_average_params(SR["x"], eff_c_shift)
        shift_h = self.get_average_params(SR["x"], eff_h_shift)

        new_values = []
        for index in range(len(self)):
            new_values.append(self.RecalculateKP(index, mode))

        for index in range(len(self)):
            Eb, Xb, me, mhh, mlh, Ncc, Nvhh, Nvlh = new_values[index]

            self[index].material.eff_mass_electron_Gamma = me
            self[index].material.eff_mass_hh_z = mhh
            self[index].material.eff_mass_lh_z = mlh

            self[index].material.__dict__["Ncc"] = 2 * (2 * pi * me * kb * self.T / h ** 2) ** 1.5
            self[index].material.__dict__["Nvhh"] = 2 * (2 * pi * mhh * kb * self.T / h ** 2) ** 1.5
            self[index].material.__dict__["Nvlh"] = 2 * (2 * pi * mlh * kb * self.T / h ** 2) ** 1.5

            self[index].band_gap = Eb
            self[index].electron_affinity = Xb

            self[index].__dict__["eff_band_gap"] = Eb + shift_e[index] - shift_h[index]
            self[index].__dict__["eff_electron_affinity"] = Xb - shift_e[index]

    # ------
    def RecalculateDensityOfStates(self, SR):
        """ Calculates the effective density of states for each layer in the QW. The general rule is:

        - Barriers have the bulk density of states
        - QW have ALL the density of states asociated with the confined states + bulk density of states above the barrier
        - Interlayers have only the bulk density of states above the barrier

        This simplification is similar to that in J. Nelson, M. Paxman, K. W. J. Barnham, J. S. Roberts, and C. Button, “Steady-state carrier escape from single quantum wells,” IEEE J. Quantum Electron., vol. 29, no. 6, pp. 1460–1468, 1993.

        From a physical porint of view, it can certainly be done much better"""

        science_reference("Density of state for QWs in the PDD solver",
                          "J. Nelson, M. Paxman, K. W. J. Barnham, J. S. Roberts, and C. Button, “Steady-state carrier escape from single quantum wells,” IEEE J. Quantum Electron., vol. 29, no. 6, pp. 1460–1468, 1993")

        # We use the average barrier electron afinity (considering only the barriers at the ends) as the electron afinity of the barrier
        Xb = (self[0].electron_affinity + self[-1].electron_affinity) / 2
        Egb = (self[0].band_gap + self[-1].band_gap) / 2

        b = q / (kb * self.T)

        e_prob = []
        hh_prob = []
        lh_prob = []
        for i in range(len(self.elevels)):
            e_prob.append(self.get_probability_per_layer(SR["x"], SR["wavefunctions"]['psi_e'][i] ** 2))
        for i in range(len(self.hhlevels)):
            hh_prob.append(self.get_probability_per_layer(SR["x"], SR["wavefunctions"]['psi_hh'][i] ** 2))
        for i in range(len(self.lhlevels)):
            lh_prob.append(self.get_probability_per_layer(SR["x"], SR["wavefunctions"]['psi_lh'][i] ** 2))

        for index in range(len(self)):
            Ncc = self[index].material.Ncc
            Nvhh = self[index].material.Nvhh
            Nvlh = self[index].material.Nvlh

            #   1- Contribution from bulk:
            Vconf = self[index].electron_affinity - Xb
            self[index].material.__dict__["Nc"] = Ncc * erfc(max(0, b * Vconf) ** 0.5)

            #   2- Contribution from confined levels
            Nqw = 0
            for i, eLevel in enumerate(self.elevels):
                Nqw = Nqw + np.exp(-b * eLevel) * e_prob[i][index]

            Nqw = 2. / self[index].width * (Ncc / 2.) ** (2. / 3.) * Nqw
            self[index].material.Nc = self[index].material.Nc + Nqw

            # For holes ------------------------------------------------------------------

            #   1- Contribution from bulk:
            Vconf = Xb + Egb - self[index].electron_affinity - self[index].band_gap
            self[index].material.__dict__["Nv"] = (Nvhh + Nvlh) * erfc(max(0, b * Vconf) ** 0.5)

            #   2- Contribution from heavy holes confined levels
            Nqw = 0
            for i, hhLevel in enumerate(self.hhlevels):
                Nqw = Nqw + (Nvhh / 2.) ** (2. / 3.) * np.exp(-b * hhLevel) * hh_prob[i][index]

                #   3- Contribution from light holes confined levels
            for i, lhLevel in enumerate(self.lhlevels):
                Nqw = Nqw + (Nvlh / 2.) ** (2. / 3.) * np.exp(-b * lhLevel) * lh_prob[i][index]

            Nqw = 2. / self[index].width * Nqw
            self[index].material.Nv = self[index].material.Nv + Nqw

    # ------
    def CalculateAbsorption(self, use_Adachi, SR):
        """ If required, this function calculates the absorption of the QW, putting together the absorption of the confined levels and the absorption of the bulk. As with the density of states, the rules are:

        - Barriers have the bulk absorption
        - Interlayers have the bulk absorption from the barrier energy and zero below that
        - Wells have the absorption of the confined levels below the barrier energy and of the bulk above it. The calculation is similar to: C. I. Cabrera, J. C. Rimada, J. P. Connolly, and L. Hernandez, “Modelling of GaAsP/InGaAs/GaAs strain-balanced multiple-quantum well solar cells,” J. Appl. Phys., vol. 113, no. 2, p. 024512, Jan. 2013."""

        science_reference("Absorption for QWs in the PDD solver",
                          "C. I. Cabrera, J. C. Rimada, J. P. Connolly, and L. Hernandez, “Modelling of GaAsP/InGaAs/GaAs strain-balanced multiple-quantum well solar cells,” J. Appl. Phys., vol. 113, no. 2, p. 024512, Jan. 2013.")

        edge = (self[0].band_gap + self[-1].band_gap) / 2 / q
        edge_index = np.abs(self.wl - 1240e-9 / edge).argmin()

        # We set the energy levels and absorption coefficient of each layer of the QW unit.
        self.absorption = 0 * self.wl
        for index in range(len(self)):
            if use_Adachi:
                try:
                    self[index].material.absorption = \
                        adachi_alpha.create_adachi_alpha(SolcoreMaterialToStr(self[index].material), T=self.T,
                                                         wl=self.wl)[
                            3]  # 0 = Energy, 1 = n, 2 = k, 3 = Absorption
                except:
                    self[index].material.absorption = self[index].material.alpha(self.wl)
                    # self[index].material.absorption[self.wl>1240e-9/sc.asUnit(self[index].material.band_gap, 'eV' )] = 0

            else:
                try:
                    self[index].material.absorption = self[index].material.alpha(self.wl)
                    # self[index].material.absorption[self.wl>1240e-9/sc.asUnit(self[index].material.band_gap, 'eV' )] = 0
                except:
                    print("Warning: Using Adachi calculation to estimate the absorption coefficient of layer: ",
                          self[index])
                    self[index].material.absorption = \
                        adachi_alpha.create_adachi_alpha(SolcoreMaterialToStr(self[index].material), T=self.T,
                                                         wl=self.wl)[
                            3]  # 0 = Energy, 1 = n, 2 = k, 3 = Absorption

            # self[index].material.absorption[ edge_index: ] = 0
            if self.labels[index] is "well":

                self[index].material.absorption = self[index].material.absorption - self[index].material.absorption[
                    edge_index]
                self[index].material.absorption[edge_index:] = 0
                self[index].material.absorption[self[index].material.absorption < 0] = 0
                self[index].material.absorption = self[index].material.absorption + SR["alpha"][1] * self[
                    index].width / self.QW_width

            elif self.labels[index] is "interlayer":
                self[index].material.absorption[edge_index:] = 0

            else:
                self[index].material.absorption[edge_index:] = 0

        self.absorption += self[index].material.absorption

    # ------
    def solve(self, Efield=0, WLsteps=(300e-9, 1100e-9, 801), wavelengths=None, periodic=True, filter_strength=0.0,
              blur=None, blurmode="left", offset=0, mode='kp8x8_bulk', use_Adachi=False, calculate_absorption=True,
              alpha_params=None, T=293):
        """ Solves the structure, calculating the energy levels, the absorption, etc. """

        if alpha_params == None:
            if wavelengths is None:
                self.wl = np.linspace(*WLsteps)
            else:
                self.wl = wavelengths

            E = 1240 / (self.wl * 1e9) * q

            alpha_params = {
                "well_width": self.QW_width,
                "theta": 0,
                "eps": 12.9 * vacuum_permittivity,
                "espace": E,
                "hwhm": si("4meV"),
                "dimensionality": 0.2,
                "line_shape": "Gauss"
            }
        else:
            self.wl = 1240 / (alpha_params["espace"] * 1e-9) * q

        SR, bands = schrodinger(self, mode=mode, periodic=True, calculate_absorption=calculate_absorption,
                                Efield=Efield, blur=None, blurmode=blurmode, alpha_params=alpha_params,
                                filter_strength=filter_strength)

        self.schrodinger = SR
        self.bands = bands

        self.__dict__["elevels"] = (SR["E"]["Ee"][:] - min(SR["potentials"]["Ve"])) / q
        self.__dict__["hhlevels"] = (max(SR["potentials"]["Vhh"]) - SR["E"]["Ehh"][:]) / q
        self.__dict__["lhlevels"] = (max(SR["potentials"]["Vlh"]) - SR["E"]["Elh"][:]) / q

        self.RecalculateBandEdges(mode, SR)
        self.RecalculateDensityOfStates(SR)
        if calculate_absorption:
            self.CalculateAbsorption(use_Adachi, SR)

        return SR

    def GetEffectiveQW(self, calculate_absorption=True, wavelengths=None, periodic=True, filter_strength=0.0, blur=None,
                       blurmode="left", mode='kp8x8_bulk', use_Adachi=False, alpha_params=None):
        """ Considers the device as a QW and solves its properties, including the modification of the bandeges due to strain, the efective mases and the absorption coefficient. The output is a list of layers made with materials with the effective properties after considering all these effects in addition to the quantum confinement.

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

        self.solve(calculate_absorption=calculate_absorption, wavelengths=wavelengths,
                   T=self.T, periodic=periodic, filter_strength=filter_strength, blur=blur,
                   blurmode=blurmode, mode=mode, use_Adachi=use_Adachi, alpha_params=alpha_params)

        for i in range(len(self)):
            self[i].material.band_gap = self[i].eff_band_gap
            self[i].material.electron_affinity = self[i].eff_electron_affinity
            self[i].material.ni = np.sqrt(
                self[i].material.Nc * self[i].material.Nv * np.exp(-self[i].eff_band_gap / (kb * self.T)))

        # Finally, we re-build a list of layers with the effective properties
        new_QW = []
        for i in range(len(self)):
            # We recover the composition and thickness
            layer_mat = self[i].material
            width = self[i].width

            # In the end, we convert the absorption coefficient in extinction coefficient
            kk = self[i].material.absorption * self.wl / 4 / np.pi
            layer_mat.k = interp1d(self.wl, kk, bounds_error=False, fill_value=(0, 0))

            # And the radiative recombination parameter
            inter = lambda E: 1.0 / layer_mat.ni ** 2 * 2 * pi / (h ** 3 * c ** 2) * layer_mat.n(
                1240e-9 / (E / q)) ** 2 * layer_mat.alphaE(E) * np.exp(-E / (kb * self.T)) * E ** 2
            Br = -np.trapz(np.nan_to_num(inter(1240e-9 / self.wl * q)), 1240e-9 / self.wl * q)
            layer_mat.radiative_recombination = Br

            # And add the layer to the list of layers
            new_QW.append(Layer(width=width, material=layer_mat))

        # As the QW might be actually a MQW, we repeat this as many times as needed
        new_QW = self.repeat * new_QW
        return new_QW


# ------
def SolcoreMaterialToStr(material_input):
    """Translate a solcore material composition into a string that is what the Adachi calculator needs."""

    material_string = material_input.__str__().strip('<>').split(" ")
    material_name = material_string[0].strip("'")
    composition = {'material': material_name}
    if len(material_name) > 4:
        material_composition = material_string[2].split("=")
        for i, comp in enumerate(material_composition):
            if comp in material_name:
                composition['element'] = material_composition[i]
                composition['fraction'] = float(material_composition[i + 1])

    output = "%s {'%s':%s}" % (composition['material'], composition['element'], composition['fraction'])

    return output
