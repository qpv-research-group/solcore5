import numpy as np

from solcore import material, si

import matplotlib.pyplot as plt

from solcore.optics.tmm import OptiStack
from solcore.optics.tmm import calculate_rat

# Import the DE implementation
from solcore.optimization import PDE, DE
from solcore.light_source import LightSource

from solcore.solar_cell import SolarCell
from solcore.structure import Junction, Layer
from solcore.solar_cell_solver import solar_cell_solver
from solcore.constants import q, kb
from solcore.material_system import create_new_material
from solcore.absorption_calculator import search_db

# The "if __name__ == "__main__" construction is used to avoid issues with parallel processing on Windows.
# The issue arises because the multiprocessing module uses a different process on Windows than on UNIX
# systems which will throw errors if this construction is not used.

# Optimizing a four-junction cell using SiGeSn as the second-lowest bandgap material.

# First, using a purely optical TMM simulation to calculate the photogenerated current in each sub-cell. The thing to
# optimize is then the current of the current-limiting cell in the structure; in other words we want to *maximize* the
# lowest sub-cell current, to achieve current-matching with the highest possible current.
# Since differential evolution does a minimization, we are actually minimizing the negative of
# this value.

# Once we have good initial values for the total layer thicknesses, use full electrical simulations to determine the
# n and p type layer thicknesses to calculate a maximum possible efficiency for the 4J device.

# To use yabox/the optimization module for the DE, we need to define a class which sets up the problem and has an 'evaluate' function, which
# will actually calculate the value we are trying to minimize for a set of parameters.

# add SiGeSn optical constants to the database
# create_new_material('SiGeSn', 'SiGeSn_n.txt', 'SiGeSn_k.txt', 'SiGeSn_params.txt') # Note: comment out this line after the material
# has been added to avoid being asked if you want to overwrite it.

# Search the refractiveindex.info for the Ta2O5 data we want to use and save the pageid to use later
Ta2O5_pageid = str(search_db("Ta2O5/Rodriguez-de Marcos")[0][0])

# define class for the optimization:
class calc_min_Jsc():

    def __init__(self):
        # initialize an instance of the class; set some information which will be used in each iteration of the calculation:
        # materials, wavelengths, the light source

        T = 298
        wl = np.linspace(300, 1900, 800)

        # materials
        SiGeSn = material('SiGeSn')(T=T)

        GaAs = material('GaAs')(T=T)
        InGaP = material('GaInP')(In=0.5, T=T)
        Ge = material('Ge')(T=T)

        # make these attributes of 'self' so they can be accessed by the class object
        # here I am also creating lists of wavelengths and corresponding n and k data from
        # the Solcore materials - the reason for this is that there is currently an issue with using the Solcore
        # material class in parallel computations. Thus the information for the n and k data is saved here as a list
        # rather than a material object (see the documentation of OptiStack for the different acceptable formats
        # to pass optical constants for an OptiStack

        self.wl = wl
        self.SiGeSn = [self.wl, SiGeSn.n(self.wl*1e-9), SiGeSn.k(self.wl*1e-9)]
        self.Ge = [self.wl, Ge.n(self.wl*1e-9), Ge.k(self.wl*1e-9)]

        self.InGaP = [self.wl, InGaP.n(self.wl*1e-9), InGaP.k(self.wl*1e-9)]
        self.GaAs = [self.wl, GaAs.n(self.wl*1e-9), GaAs.k(self.wl*1e-9)]
        self.MgF2 = [self.wl, material('MgF2')().n(self.wl*1e-9), material('MgF2')().k(self.wl*1e-9)]
        self.Ta2O5 = [self.wl, material(Ta2O5_pageid,
                                        nk_db=True)().n(self.wl*1e-9), material(Ta2O5_pageid,
                                                                                nk_db=True)().k(self.wl*1e-9)]

        # assuming an AM1.5G spectrum
        self.spectr = LightSource(source_type='standard', version='AM1.5g', x=self.wl,
                           output_units='photon_flux_per_nm', concentration=1).spectrum(self.wl)[1]





    def evaluate(self, x):

        # x[0] = MgF2 thickness (anti-reflection coating)
        # x[1] = Ta2O5 thickness (anti-reflection coating)
        # x[2]  = InGaP (top junction) thickness
        # x[3] = GaAs (second junction) thickness
        # x[4] = SiGeSn (third junction) thickness

        # keep the thickness of the bottom cell constant; from an optical point of view, this should be infinitely thick

        SC = [[x[0]] + self.MgF2, [x[1]] + self.Ta2O5, [x[2]] + self.InGaP, [x[3]] + self.GaAs, [x[4]] + self.SiGeSn, [300e3]+self.Ge]#, [x[5]] + self.Ge]

        # create the OptiStack
        full_stack = OptiStack(SC, no_back_reflection=False)

        # calculate reflection, transmission, and absorption in each layer. We are specifying that the last layer,
        # a very thick Ge substrate, should be treated incoherently, otherwise we would see very narrow, unphysical oscillations
        # in the R/A/T spectra.
        RAT = calculate_rat(full_stack, self.wl, no_back_reflection=False, coherent=False,
                            coherency_list=['c', 'c', 'c', 'c', 'c', 'i'])

        # extract absorption per layer
        A_InGaP = RAT['A_per_layer'][3]
        A_GaAs = RAT['A_per_layer'][4]
        A_SiGeSn = RAT['A_per_layer'][5]
        A_Ge = RAT['A_per_layer'][6]

        ## calculate photo-generated currents using the AM1.5 G spectrum for each layer
        Jsc_InGaP = 0.1 * q * np.trapz(A_InGaP* self.spectr, self.wl)
        Jsc_GaAs = 0.1 * q * np.trapz(A_GaAs * self.spectr, self.wl)
        Jsc_SiGeSn = 0.1 * q * np.trapz(A_SiGeSn* self.spectr, self.wl)
        Jsc_Ge = 0.1 * q* np.trapz(A_Ge * self.spectr, self.wl)

        # find the limiting current by checking which junction has the lowest current. Then take the negative since
        # we need to minimize (not maximize)
        limiting_Jsc = -min([Jsc_InGaP, Jsc_GaAs, Jsc_SiGeSn, Jsc_Ge])

        return limiting_Jsc

    def plot(self, x):

        # this does basically what evaluate() does, but plots the results
        SC = [[x[0]] + self.MgF2, [x[1]] + self.Ta2O5, [x[2]] + self.InGaP,
              [x[3]] + self.GaAs, [x[4]] + self.SiGeSn, [300e3] + self.Ge]

        full_stack = OptiStack(SC, no_back_reflection=False)

        RAT = calculate_rat(full_stack, self.wl, no_back_reflection=False,
                            coherent=False, coherency_list=['c', 'c', 'c', 'c', 'c', 'i'])

        A_InGaP = RAT['A_per_layer'][3]
        A_GaAs = RAT['A_per_layer'][4]
        A_SiGeSn = RAT['A_per_layer'][5]
        A_Ge = RAT['A_per_layer'][6]

        plt.figure()
        plt.plot(self.wl, A_InGaP, label='InGaP')
        plt.plot(self.wl, A_GaAs, label='A_GaAs')
        plt.plot(self.wl, A_SiGeSn, label='SiGeSn')
        plt.plot(self.wl, A_Ge, label = 'Ge')
        plt.plot(self.wl, RAT['T'], label='T')
        plt.plot(self.wl, RAT['R'], label='R')
        plt.legend()
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('R/A/T')
        plt.show()


class calc_min_Jsc_DA():

    def __init__(self, ARC_thickness):
        self.ARC = ARC_thickness

    def make_cell(self, x):

        #x[0]: total InGaP thickness
        #x[1]: total InGaAs thickness
        #x[2]: total SiGeSn thickness
        #x[3]: total Ge thickness

        #x[4]: InGaP n thickness
        #x[5]: InGaAs n thickness
        #x[6]: SiGeSn n thickness
        #x[7]: Ge n thickness

        e_charge = si('1eV')

        # materials
        SiGeSn = material('SiGeSn')

        GaAs = material('GaAs')
        InGaP = material('GaInP')
        Ge = material('Ge')
        MgF2 = material('MgF2')()
        Ta2O5 = material(Ta2O5_pageid, nk_db=True)()
        AlInP = material("AlInP")

        window_material = AlInP(Al=0.52)

        GaInP_mobility_h = 0.03  #
        GaInP_lifetime_h = 1e-8
        GaInP_D_h = GaInP_mobility_h * kb * 300 / e_charge
        GaInP_L_h = np.sqrt(GaInP_D_h * GaInP_lifetime_h)
        GaInP_mobility_e = 0.015
        GaInP_lifetime_e = 1e-8
        GaInP_D_e = GaInP_mobility_e * kb * 300 / e_charge
        GaInP_L_e = np.sqrt(GaInP_D_e * GaInP_lifetime_e)

        top_cell_n_material = InGaP(In=0.5, Nd=si("2e18cm-3"), hole_diffusion_length=GaInP_L_h, hole_mobility=GaInP_mobility_h)
        top_cell_p_material = InGaP(In=0.5, Na=si("2e17cm-3"), electron_diffusion_length=GaInP_L_e, electron_mobility=GaInP_mobility_e)

        # MID CELL  - GaAs

        GaAs_mobility_h = 0.85  #
        GaAs_lifetime_h = 1e-8
        GaAs_D_h = GaAs_mobility_h * kb * 300 / e_charge
        GaAs_L_h = np.sqrt(GaAs_D_h * GaAs_lifetime_h)
        GaAs_mobility_e = 0.08
        GaAs_lifetime_e = 1e-8
        GaAs_D_e = GaAs_mobility_e * kb * 300 / e_charge
        GaAs_L_e = np.sqrt(GaAs_D_e * GaAs_lifetime_e)

        mid_cell_n_material = GaAs(Nd=si("1e18cm-3"), hole_diffusion_length=GaAs_L_h, hole_mobility=GaAs_mobility_h)
        mid_cell_p_material = GaAs(Na=si("1e17cm-3"), electron_diffusion_length=GaAs_L_e, electron_mobility=GaAs_mobility_e)


        SiGeSn.band_gap = si('0.77eV') # from PL measurement
        SiGeSn_L_h = si('0.35um')
        SiGeSn_L_e = si('5um')
        SiGeSn_lifetime_e = 1e-6
        SiGeSn_lifetime_h = 1e-6
        SiGeSn_mobility_h = SiGeSn_L_h ** 2 * e_charge / (SiGeSn_lifetime_h * kb * 300)
        SiGeSn_mobility_e = SiGeSn_L_e ** 2 * e_charge / (SiGeSn_lifetime_e * kb * 300)

        pen_cell_n_material = SiGeSn(Nd=si("1e18cm-3"), hole_diffusion_length=SiGeSn_L_h,
                                     relative_permittivity=16, hole_mobility=SiGeSn_mobility_h)
        pen_cell_p_material = SiGeSn(Na=si("1e17cm-3"), electron_diffusion_length=SiGeSn_L_e,
                                     relative_permittivity=16, electron_mobility=SiGeSn_mobility_e)

        # Ge_mobility_h = 0.38 #
        Ge_lifetime_h = 1e-6
        # Ge_D_h = Ge_mobility_h*kb*300/e_charge
        # Ge_L_h = np.sqrt(Ge_D_h*Ge_lifetime_h)
        Ge_L_h = si('500nm')
        Ge_mobility_h = Ge_L_h ** 2 * e_charge / (Ge_lifetime_h * kb * 300)
        Ge_mobility_e = 0.18
        Ge_lifetime_e = 1e-6
        Ge_D_e = Ge_mobility_e * kb * 300 / e_charge
        Ge_L_e = np.sqrt(Ge_D_e * Ge_lifetime_e)


        bot_cell_n_material = Ge(Nd=si("2e18cm-3"), hole_diffusion_length=Ge_L_h, hole_mobility=Ge_mobility_h)
        bot_cell_p_material = Ge(Na=si("1e17cm-3"), electron_diffusion_length=Ge_L_e, electron_mobility=Ge_mobility_e)


        # And, finally, we put everything together, adding also the surface recombination velocities. We also add some shading
        # due to the metallisation of the cell = 8%, and indicate it has an area of 0.7x0.7 mm2 (converted to m2)
        solar_cell = SolarCell([
            # Layer(si('110nm'), material = MgF2), Layer(si('55nm'), material = ZnS),
            Layer(si(self.ARC[0], 'nm'), material=MgF2), Layer(si(self.ARC[1], 'nm'), material=Ta2O5),
            Junction([Layer(si(25, 'nm'), material=window_material, role='window'),
                      Layer(si(x[4], 'nm'), material=top_cell_n_material, role='emitter'),
                      Layer(si(x[0]-x[4], 'nm'), material=top_cell_p_material, role='base'),
                      ], sn=1, sp=1, kind='DA'),
            Junction([Layer(si(x[5], 'nm'), material=mid_cell_n_material, role='emitter'),
                      Layer(si(x[1]-x[5], 'nm'), material=mid_cell_p_material, role='base'),
                      ], sn=1, sp=1, kind='DA'),
            Junction([Layer(si(x[6], 'nm'), material=pen_cell_n_material, role='emitter'),
                      Layer(si(x[2]-x[6], 'nm'), material=pen_cell_p_material, role='base'),
                      ], sn=1, sp=1, kind='DA'),
            Junction([Layer(si(x[7], 'nm'), material=bot_cell_n_material, role='emitter'),
                      Layer(si(x[3]-x[7], 'nm'), material=bot_cell_p_material, role='base'),
                      ], sn=1, sp=1, kind='DA'),
        ], shading=0.0, substrate=bot_cell_n_material)

        return solar_cell

    def evaluate(self, x):

        light_source = LightSource(source_type='standard', version='AM1.5g')

        wl = np.linspace(300, 1850, 500) * 1e-9

        solar_cell = self.make_cell(x)

        position = [1e-10] * 10 + [5e-8]

        V = np.linspace(0, 3.5, 300)
        solar_cell_solver(solar_cell, 'iv',
                          user_options={'voltages': V, 'light_iv': True, 'wavelength': wl, 'mpp': True,
                                        'light_source': light_source,
                                        'optics_method': 'TMM', 'BL_correction': True, 'position': position})

        efficiency = solar_cell.iv["Eta"]

        # print('Efficiency =', efficiency)

        return -efficiency

    def plot(self, x):
        light_source = LightSource(source_type='standard', version='AM1.5g')

        wl = np.linspace(300, 1850, 500) * 1e-9

        solar_cell = self.make_cell(x)

        position = [1e-10] * 10 + [5e-8]

        V = np.linspace(0, 3.5, 300)
        solar_cell_solver(solar_cell, 'iv',
                          user_options={'voltages': V, 'light_iv': True, 'wavelength': wl, 'mpp': True,
                                        'light_source': light_source,
                                        'optics_method': 'TMM', 'BL_correction': True, 'position': position})

        efficiency = solar_cell.iv["Eta"]
        pmax = solar_cell.iv["Pmpp"]
        ff = solar_cell.iv["FF"]
        voc = solar_cell.iv["Voc"]
        isc = solar_cell.iv["Isc"]

        plt.figure()

        plt.plot(V, solar_cell.iv['IV'][1] / 10, 'k', linewidth=3, label='Total')
        plt.plot(V, -solar_cell[2].iv(V) / 10, 'b', label='GaInP')
        plt.plot(V, -solar_cell[3].iv(V) / 10, 'g', label='GaAs')
        plt.plot(V, -solar_cell[4].iv(V) / 10, 'r', label='SiGeSn')
        plt.plot(V, -solar_cell[5].iv(V) / 10, 'y', label='Ge')
        plt.text(2, 10, '$\eta = $' + str(round(efficiency * 100, 1)) + '%')
        plt.text(2, 8,'Pmax='+str(round(pmax,1))+'W/m$^2$')
        plt.text(2, 9, 'FF = ' + str(round(ff * 100, 1)) + '%')
        plt.text(2,7,'Voc='+str(round(voc,1))+'V')
        plt.text(2,6, 'Jsc='+str(round(0.1*isc,1))+'mA/cm$^2$')

        plt.legend()
        plt.ylim(0, 18)
        plt.xlim(0, 3.5)
        plt.ylabel('Current (mA/cm$^2$)')
        plt.xlabel('Voltage (V)')

        plt.show()

        # print('Efficiency =', efficiency)

        solar_cell_solver(solar_cell, 'qe',
                         user_options={'wavelength': wl, 'optics_method': 'TMM', 'BL_correction': True, 'position': position})

        plt.figure()
        plt.plot(wl * 1e9, solar_cell[2].eqe(wl) * 100, 'b', label='InGaP')
        plt.plot(wl * 1e9, solar_cell[3].eqe(wl) * 100, 'g', label='InGaAs')
        plt.plot(wl * 1e9, solar_cell[4].eqe(wl) * 100, 'r', label='SiGeSn')
        plt.plot(wl * 1e9, solar_cell[5].eqe(wl) * 100, 'y', label='Ge')
        plt.plot(wl * 1e9, solar_cell.absorbed * 100, 'k--', label='Absorption')
        # plt.plot(wl * 1e9, solar_cell[5].eqe(wl)*100, 'y', label='Ge')

        plt.legend(loc='upper right')
        plt.xlim(290, 1850)
        plt.ylim(0, 100)
        plt.ylabel('EQE (%)')
        plt.xlabel('Wavelength (nm)')
        plt.show()


def main():
    # number of iterations for Differential Evolution optimization of the optical stack
    maxiters=300

    # make an instance of the class the DE algorithm is going to use, as defined above
    DE_class = calc_min_Jsc()

    # pass the function which will be minimized to the PDE (parallel differential evolution) solver. PDE calculates the
    # results for each population in parallel to speed up the overall process. The bounds argument sets upper and lower bounds
    # for each parameter. PDE_obj contains all the information to run the DE but does not actually invoke the calculation....
    PDE_obj = PDE(DE_class.evaluate, bounds=[[10,150], [10,105], [200, 1000], [500, 10000], [500, 10000]], maxiters=maxiters)

    # this will run the calculation in parallel, with all the cores available. If you don't want this, use 'DE' instead of 'PDE'

    # to run the DE, use the .solve() function of the PDE object class
    res = PDE_obj.solve()

    # PDE_obj.solve() returns 5 things:
    # res[0] is a list of the parameters which gave the minimized value
    # res[1] is that minimized value
    # res[2] is the evolution of the best population (the best population from each iteration
    # res[3] is the evolution of the minimized value, i.e. the fitness over each iteration
    # res[4] is the evolution of the mean fitness over the iterations

    # best population:
    best_pop = res[0]

    print('parameters for best result:', best_pop, '\n', 'optimized Jsc value (mA/cm2):', -res[1])

    # plot the result at these best parameters
    DE_class.plot(best_pop)

    best_pop_evo = res[2]
    best_fitn_evo = res[3]
    mean_fitn_evo = res[4]
    final_fitness = res[1]

    # plot evolution of the fitness of the best population per iteration

    plt.figure()
    plt.plot(-best_fitn_evo, '-k')
    plt.xlabel('iteration')
    plt.ylabel('fitness')
    plt.title('Best fitness')
    plt.show()

    # plot evolution of the mean fitness of the population per iteration

    plt.figure()
    plt.plot(-mean_fitn_evo, '-k')
    plt.xlabel('iteration')
    plt.ylabel('fitness')
    plt.title('Mean fitness')
    plt.show()

    # these plots show that the fitness of the best population 'jumps' every few iterations as a new best population is found,
    # while the mean fitness converges slowly as the whole population gradually improves


    ## Now that the layer thicknesses have been optimized from an optical point of view, we want to design the device (or
    # at least a simplified version, by calculating a more realistic EQE. Obviously additional parameters like the doping of the
    # layers could be varied too.

    # x[0]: total InGaP thickness
    # x[1]: total InGaAs thickness
    # x[2]: total SiGeSn thickness
    # x[3]: total Ge thickness

    # x[4]: InGaP n thickness
    # x[5]: InGaAs n thickness
    # x[6]: SiGeSn n thickness
    # x[7]: Ge n thickness

    # keep the ARC thicknesses fixed at the values obtained in the optical simulation

    # generate upper and lower bounds: total layer thickness between 75% and 125% of values fitted in TMM calculation. Ge
    # starting value is 200 um
    starting_params = np.append(best_pop[2:], [200000])

    lower = 0.75*starting_params
    upper = 1.25*starting_params

    # upper and lower bounds for the n-type (highly doped) layers
    lower_ntype = [20, 20, 20, 20]

    upper_ntype = [200, 300, 300, 500]

    all_lower = np.append(lower, lower_ntype)

    all_upper = np.append(upper, upper_ntype)

    # full list of bounds
    all_bounds = np.stack((all_lower, all_upper)).T

    # DE calculation for the electrical simulation
    maxiters_DA = 10
    DE_class_DA = calc_min_Jsc_DA(best_pop[0:2])


    # default population size = 5*number of params
    PDE_obj_DA = PDE(DE_class_DA.evaluate, bounds=all_bounds, maxiters=maxiters_DA)


    # solve, i.e. minimize the problem
    res_DA = PDE_obj_DA.solve()

    best_pop_DA = res_DA[0]

    print('parameters for best result:', best_pop_DA, 'optimized efficiency (%)', res_DA[1]*100)

    # plot the result at these best parameters
    DE_class_DA.plot(best_pop_DA)

    best_pop_evo = res_DA[2]
    best_fitn_evo = res_DA[3]
    mean_fitn_evo = res_DA[4]
    final_fitness = res_DA[1]

    # plot evolution of the fitness of the best population per iteration, and the mean fitness

    plt.figure()
    plt.plot(-best_fitn_evo, '-k')
    plt.xlabel('iteration')
    plt.ylabel('fitness')
    plt.title('Best fitness')
    plt.show()

    plt.figure()
    plt.plot(-mean_fitn_evo, '-k')
    plt.xlabel('iteration')
    plt.ylabel('fitness')
    plt.title('Mean fitness')
    plt.show()


if __name__ == '__main__':
    main()