import numpy as np
from scipy.special import erfc
from solcore import constants, material, si
from solcore.structure import Layer, Structure
from solcore.quantum_mechanics import schrodinger
from solcore.quantum_mechanics.kp_bulk import kp8x8_bulk
from solcore.quantum_mechanics.strain import strain_calculation_parameters
from solcore.absorption_calculator import adachi_alpha
import sys

q = constants.q
pi = constants.pi
h = constants.h
kb = constants.kb
m0 = constants.electron_mass
vacuum_permittivity = constants.vacuum_permittivity


class QWunit(Structure):
    """ Asembles a group of layers as a quantum well, calculating its properties as a whole
    
    - It needs a minimum of 3 layers
    - There must be exactly one layer with the role = "well".
    - Top and bottom layers are assumed to be the barriers.
    - Anything not a barrier not a well will be considered as "interlayer".    
    
    Since the Poisson - Drift-Diffusion solver can not work with superlatices, there is no point of considering the solution of more than one QW. 
    All QWs entering into the solver are independent, regardless of the width of the barriers.
    
    To account for the possible presence of nearby QWs, perdiodic boundary conditions can be used to solve the Schrodinger equation. 
    If barriers + interlayers are thick enought, this make no difference but if they are thin, it affects the number of confined levels. 
    It does not 'create' minibands, though, so the results with thin barriers must be used with caution.  
    
    """

    def __init__(self, *args, **kwargs):
        """ Initialises the QWunit, checking that it has all the necessary layers """
        super().__init__(*args, **kwargs)

        if len(self) < 3:
            print("ERROR creating a QWunit: a minimum of 3 layers is necesarry to make a unit. Only {0} found. ".format(
                len(self)))
            sys.exit()

        self.QW_number = 0
        self.QW_width = 0

        self[0].__dict__["role"] = "barrier"
        self[-1].__dict__["role"] = "barrier"
        for index in range(1, len(self) - 1):
            layer = self[index]

            if self[index].role in ["well", "qw", "QW", "Qw", "Well"]:
                self.QW_number += 1
                self.QW_width = self[index].width
                self[index].role = "well"
            else:
                self[index].__dict__["role"] = "interlayer"

        if self.QW_number != 1:
            print("ERROR creating a QWunit: There must be exactly one layer with the role = 'well'. {0} found. ".format(
                self.QW_number))
            sys.exit()

        if 'substrate' not in kwargs.keys():
            print("ERROR creating a QWunit: There must be a substrate defined with the keyword 'substrate'. ")
            sys.exit()

        # We update the structure of labels since that is what the quantum mechanics modules uses to identify the regions, not the role. 
        for index in range(len(self)):
            self.labels[index] = self[index].role

    # ------
    def RecalculateKP(self, index, mode):
        """ Recalculates the band structure of the QW unit to account for modifications due to strain. 
        
        This is already done in the Schrodinger solver, but it does not produce that information as output. We do it so here.
        """
        if mode == 'kp8x8_bulk':
            kp_result = kp8x8_bulk(material=self[index].material, host_material=self.substrate, averaging_points=6)
            c, hh, lh, so, me, mhh, mlh, mso = kp_result

            self[index].material.eff_mass_electron_Gamma = me
            self[index].material.eff_mass_hh_z = mhh
            self[index].material.eff_mass_lh_z = mlh

            # We recalculate the bandgap adn the electron afinity
            Ec = self[index].material.valence_band_offset + self[index].material.band_gap
            self[index].band_gap = c - max(hh, lh)
            self[index].electron_affinity = self[index].material.electron_affinity - c + Ec
        else:
            me = self[index].material.eff_mass_electron_Gamma * m0
            mhh = self[index].material.eff_mass_hh_z * m0
            mlh = self[index].material.eff_mass_lh_z * m0

            # Idem
            strain_parameters = strain_calculation_parameters(self.substrate, self[index].material, should_print=False)
            self[index].electron_affinity = self[index].material.electron_affinity - strain_parameters.delta_Ec
            self[index].band_gap = self[index].material.band_gap + strain_parameters.delta_Ec + max(
                strain_parameters.delta_Ehh, strain_parameters.delta_Elh)

        # Calculating this properly for QWs is done in the "RecalculateDensityOfStates" routine      
        self[index].material.__dict__["Ncc"] = 2 * (2 * pi * me * kb * self.T / h ** 2) ** 1.5
        self[index].material.__dict__["Nvhh"] = 2 * (2 * pi * mhh * kb * self.T / h ** 2) ** 1.5
        self[index].material.__dict__["Nvlh"] = 2 * (2 * pi * mlh * kb * self.T / h ** 2) ** 1.5

        # ------

    def RecalculateBandEdges(self, mode, SR):
        """The effective bandgap and effective affinity of each layer depends on the energy level, e.g.: the minimum energy for electrons is not the band edge of the conduction band, but the ground level. In this function we modify that creating an effective electron afinity and band gaps for all layers in the QW.
        # 
        # For the barriers, the electron afinity and band gap is the same than in bulk, modified by the kp calculation, if necesary."""

        for index in range(len(self)):
            self.RecalculateKP(index, mode)

            self[index].__dict__["elevels"] = []
            self[index].__dict__["hhlevels"] = []
            self[index].__dict__["lhlevels"] = []

            self[index].__dict__["eff_band_gap"] = self[index].band_gap
            self[index].__dict__["eff_electron_affinity"] = self[index].electron_affinity

            if self.labels[index] is "well":

                self[index].elevels = self.elevels
                self[index].hhlevels = self.hhlevels
                self[index].lhlevels = self.lhlevels

                # We create an effective bandgap and electron afinitiy to account for the confinement energies of electrons and holes
                self[index].eff_band_gap = self[index].eff_band_gap + self.elevels[0] * q
                self[index].eff_electron_affinity = self[index].eff_electron_affinity - self.elevels[0] * q

                # if ( max(schrodinger_result["potentials"]["Vhh"]) > max(schrodinger_result["potentials"]["Vlh"]) ): 
                if (SR["E"]["Ehh"][0] > SR["E"]["Elh"][0]):
                    self[index].eff_band_gap = self[index].eff_band_gap + self.hhlevels[0] * q
                else:
                    self[index].eff_band_gap = self[index].eff_band_gap + self.lhlevels[0] * q

                well_eff_electron_afinity = self[index].eff_electron_affinity
                well_eff_band_gap = self[index].eff_band_gap

        # We re-run the loop to correct the electron affinity of the interlayers    
        for index, layer in enumerate(self):
            if self.labels[index] is "interlayer":
                # The effective electron afinity is whatever is lower, its own afinity or the effective affinity of the well
                self[index].eff_electron_affinity = min(self[index].eff_electron_affinity, well_eff_electron_afinity)
                # The bandgap is corrected so the valence band in the interlayer is always equal or lower than in the QW
                gap_offset = (well_eff_electron_afinity + well_eff_band_gap) - (
                    self[index].eff_electron_affinity + self[index].eff_band_gap)
                if gap_offset > 0:
                    self[index].eff_band_gap = self[index].eff_band_gap + gap_offset

    # ------
    def RecalculateDensityOfStates(self):
        """ Calculates the effective density of states for each layer in the QW. The general rule is:
        #   1- Barriers have the bulk density of states
        #   2- QW have ALL the density of states asociated with the confined states + bulk density of states above the barrier
        #   3- Interlayers have only the bulk density of states above the barrier
        #
        # This simplification is similar to that in J. Nelson, M. Paxman, K. W. J. Barnham, J. S. Roberts, and C. Button, “Steady-state carrier escape from single quantum wells,” IEEE J. Quantum Electron., vol. 29, no. 6, pp. 1460–1468, 1993.
        #
        # From a physical porint of view, it can certainly be done much better"""

        # We use the average barrier electron afinity (considering only the barriers at the ends) as the electron afinity of the barrier
        Xb = (self[0].electron_affinity + self[-1].electron_affinity) / 2
        Egb = (self[0].band_gap + self[-1].band_gap) / 2

        b = q / (kb * self.T)

        for index in range(len(self)):
            Ncc = self[index].material.Ncc
            Nvhh = self[index].material.Nvhh
            Nvlh = self[index].material.Nvlh
            if self.labels[index] is "barrier":
                # Barriers have the bulk density of states
                self[index].material.__dict__["Nc"] = Ncc
                self[index].material.__dict__["Nv"] = Nvhh + Nvlh

            elif self.labels[index] is "interlayer":
                # For electrons  -------------------------------------------------------------
                Vconf = self[index].electron_affinity - Xb
                self[index].material.__dict__["Nc"] = Ncc * erfc((b * Vconf) ** 0.5)

                # For holes ------------------------------------------------------------------
                Vconf = Xb + Egb - self[index].electron_affinity - self[index].band_gap
                self[index].material.__dict__["Nv"] = (Nvhh + Nvlh) * erfc((b * Vconf) ** 0.5)

            else:
                # For electrons  -------------------------------------------------------------

                #   1- Contribution from bulk:
                Vconf = self[index].electron_affinity - Xb
                self[index].material.__dict__["Nc"] = Ncc * erfc((b * Vconf) ** 0.5)

                #   2- Contribution from confined levels
                Nqw = 0
                for eLevel in self[index].elevels:
                    Nqw = Nqw + np.exp(-b * eLevel)

                Nqw = 2. / self[index].width * (Ncc / 2.) ** (2. / 3.) * Nqw
                self[index].material.Nc = self[index].material.Nc + Nqw

                # For holes ------------------------------------------------------------------

                #   1- Contribution from bulk:
                Vconf = Xb + Egb - self[index].electron_affinity - self[index].band_gap
                self[index].material.__dict__["Nv"] = (Nvhh + Nvlh) * erfc((b * Vconf) ** 0.5)

                #   2- Contribution from heavy holes confined levels
                Nqw = 0
                for hhLevel in self[index].hhlevels:
                    Nqw = Nqw + (Nvhh / 2.) ** (2. / 3.) * np.exp(-b * hhLevel)

                    #   3- Contribution from light holes confined levels
                for lhLevel in self[index].lhlevels:
                    Nqw = Nqw + (Nvlh / 2.) ** (2. / 3.) * np.exp(-b * lhLevel)

                Nqw = 2. / self[index].width * Nqw
                self[index].material.Nv = self[index].material.Nv + Nqw

    # ------
    def CalculateAbsorption(self, use_Adachi, SR):
        """ If required, this function calculates the absorption of the QW, putting together the absorption of the confined levels and the absorption of the bulk. As with the density of states, the rules are:
        # 
        #   1- Barriers have the bulk absorption
        #   2- Interlayers have the bulk absorption from the barrier energy and zero below that
        #   3- Wells have the absorption of the confined levels below the barrier energy and of the bulk above it. The calculation is similar to: C. I. Cabrera, J. C. Rimada, J. P. Connolly, and L. Hernandez, “Modelling of GaAsP/InGaAs/GaAs strain-balanced multiple-quantum well solar cells,” J. Appl. Phys., vol. 113, no. 2, p. 024512, Jan. 2013."""

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
    def solve(self, Efield=0, WLsteps=(300e-9, 1100e-9, 201), wavelengths=None, periodic=True, filter_strength=0.0,
              blur=None, blurmode="left", offset=0, mode='kp8x8_bulk', use_Adachi=False, calculate_absorption=True,
              alpha_params=None, T=293):
        """ Solves the structure, calculating the energy levels, the absorption, etc. """

        self.T = T

        if calculate_absorption:
            # If we want the absorption but have not provided enough data for that, we make up some default values
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
                                Efield=Efield,
                                blur=None, blurmode=blurmode, alpha_params=alpha_params, filter_strength=filter_strength)
        self.__dict__["elevels"] = (SR["E"]["Ee"][:] - min(SR["potentials"]["Ve"])) / q
        self.__dict__["hhlevels"] = (max(SR["potentials"]["Vhh"]) - SR["E"]["Ehh"][:]) / q
        self.__dict__["lhlevels"] = (max(SR["potentials"]["Vlh"]) - SR["E"]["Elh"][:]) / q

        self.RecalculateBandEdges(mode, SR)
        self.RecalculateDensityOfStates()
        if calculate_absorption: self.CalculateAbsorption(use_Adachi, SR)

        return SR

        # ------


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
