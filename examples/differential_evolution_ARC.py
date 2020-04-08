import numpy as np

from solcore import material
import matplotlib.pyplot as plt

from solcore.optics.tmm import OptiStack
from solcore.optics.tmm import calculate_rat

# Import the DE implementations
from solcore.optimization import PDE, DE
from solcore.light_source import LightSource

# Optimizing a double-layer MgF2/Ta2O5 anti-reflection coating for "infinitely-thick" GaAs. Minimize reflection * AM0 spectrum
# (weighted reflectance).

# To use yabox for the DE, we need to define a class which sets up the problem and has an 'evaluate' function within it, which
# will actually calculate the value we are trying to minimize for each set of parameters.

class calc_R_diff():

    def __init__(self):

        # wavelength, materials.
        # make these attributes of 'self' so they can be accessed by the class object
        # here I am also creating lists of wavelengths and corresponding n and k data from
        # the Solcore materials - the reason for this is that there is currently an issue with using the Solcore
        # material class in parallel computations. Thus the information for the n and k data is saved here as a list
        # rather than a material object.
        self.wl = np.linspace(250, 900, 700)
        self.MgF2 = [self.wl, material('MgF2')().n(self.wl*1e-9), material('MgF2')().k(self.wl*1e-9)]
        self.Ta2O5 = [self.wl, material('TAOX1')().n(self.wl*1e-9), material('TAOX1')().k(self.wl*1e-9)] # this is Ta2O5 from the SOPRA database included in Solcore
        self.GaAs = [1000, self.wl, material('GaAs')().n(self.wl*1e-9), material('GaAs')().k(self.wl*1e-9)]

        # assuming an AM0 spectrum
        spectr = LightSource(source_type='standard', version='AM1.5g', x=self.wl,
                           output_units='photon_flux_per_m', concentration=1).spectrum(self.wl*1e-9)[1]

        self.spectrum = spectr/max(spectr) # only want to use spectrum to weight reflectivity result so don't care about
        # absolute values



    def evaluate(self, x):
        # create a list of layers with the format [thickness, wavelengths, n_data, k_data] for each layer.
        # This is one of the acceptable formats in which OptiStack can take information (look at the Solcore documentation
        # or at the OptiStack code for more info

        ARC = [[x[0]] + self.MgF2, [x[1]] + self.Ta2O5, self.GaAs]

        # create the OptiStack. We set no_back_reflection to True because we DO  NOT want to include reflection at the back surface
        # (assume GaAs is infinitely thick)

        full_stack = OptiStack(ARC, no_back_reflection=True)

        # calculate reflection, transmission, and absorption in each layer.

        RAT = calculate_rat(full_stack, self.wl, no_back_reflection=True)

        R_weighted = RAT['R']*self.spectrum

        # 'evaluate' should return a single number. The DE algorithm tries to minimize this number.

        return sum(R_weighted)

    def plot(self, x):
        # this does basically what evaluate() does, but plots the reflectivity

        ARC = [[x[0]] + self.MgF2, [x[1]] + self.Ta2O5, self.GaAs]

        full_stack = OptiStack(ARC, no_back_reflection=True)

        RAT = calculate_rat(full_stack, self.wl, no_back_reflection=True)
        plt.figure()
        plt.plot(self.wl, RAT['R'])
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('R')
        plt.show()

    def plot_weighted(self, x):
        # this does basically what evaluate() does, but plots the reflectivity weighted by the AM0 solar spectrum
        ARC = [[x[0]] + self.MgF2, [x[1]] + self.Ta2O5, self.GaAs]

        full_stack = OptiStack(ARC, no_back_reflection=True)

        RAT = calculate_rat(full_stack, self.wl, no_back_reflection=True)
        spect_n = self.spectrum/np.max(self.spectrum)
        plt.figure()
        plt.plot(self.wl, RAT['R']*spect_n)
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('R weighted by AM0')
        plt.show()


# number of iterations for Differential Evolution
maxiters=100

# class the DE algorithm is going to use, as defined above
PDE_class = calc_R_diff()

# pass the function which will be minimized to the PDE (parallel differential evolution) solver. PDE calculates the
# results for each population in parallel to speed up the overall process

PDE_obj = PDE(PDE_class.evaluate, bounds=[[0,400], [0,400]], maxiters=maxiters)

# solve, i.e. minimize the problem
res = PDE_obj.solve()

# PDE_obj.solve() returns 5 things:
# res[0] is a list of the parameters which gave the minimized value
# res[1] is that minimized value
# res[2] is the evolution of the best population (the best population from each iteration
# res[3] is the evolution of the minimized value, i.e. the fitness over each iteration
# res[4] is the evolution of the mean fitness over the iterations

best_pop = res[0]
print('parameters for best result:', best_pop, res[1])

PDE_class.plot(best_pop)
PDE_class.plot_weighted(best_pop)


best_pop_evo = res[2]
best_fitn_evo = res[3]
mean_fitn_evo = res[4]
final_fitness = res[1]


# plot evolution of the fitness of the best population per iteration

plt.figure()
plt.plot(best_fitn_evo, '-k')
plt.xlabel('iteration')
plt.ylabel('fitness')
plt.title('Best fitness')
plt.show()

plt.figure()
plt.plot(mean_fitn_evo, '-k')
plt.xlabel('iteration')
plt.ylabel('fitness')
plt.title('Mean fitness')
plt.show()

