"""
This is completely based on Pablo R. Mier's yaboc package, with modifications to slightly change what is retured by
DE.solve to be able to visualize the convergence of the DE without having to use de.iterator, and just return
a list of optimized parameters rather than a reshaped array.

All credit of the algorithm, testing, etc. goes to Pablo R. Mier. For more details, visit:
- https://github.com/pablormier/yabox
- https://pablormier.github.io/2017/09/05/a-tutorial-on-differential-evolution-with-python/
"""

from yabox.algorithms.de import DEIterator, PDEIterator
from yabox.algorithms.base import *

class DE:
    def __init__(self, fobj, bounds, mutation=(0.5, 1.0), crossover=0.7, maxiters=1000,
                 self_adaptive=False, popsize=None, seed=None):
        self.adaptive = self_adaptive
        # Indicates the number of extra parameters in an individual that are not used for evaluating
        # If extra_params = d, discards the last d elements from an individual prior to evaluation.
        self.extra_params = 0
        # Convert crossover param to an interval, as in mutation. If min/max values in the interval are
        # different, a dither mechanism is used for crossover (although this is not recommended, but still supported)
        # TODO: Clean duplicate code

        self.crossover_bounds = crossover
        self.mutation_bounds = mutation

        if getattr(crossover, '__len__', None) is None:
            self.crossover_bounds = [crossover, crossover]

        if getattr(mutation, '__len__', None) is None:
            self.mutation_bounds = [mutation, mutation]

        # If self-adaptive, include mutation and crossover as two new variables
        bnd = list(bounds)
        if self_adaptive:
            bnd.append(self.mutation_bounds)
            bnd.append(self.crossover_bounds)
            self.extra_params = 2
        self._MIN, self._MAX = np.asarray(bnd, dtype='f8').T
        self._DIFF = np.fabs(self._MAX - self._MIN)
        self.dims = len(bnd)
        self.fobj = fobj
        self.maxiters = maxiters
        if popsize is None:
            self.popsize = self.dims * 5
        else:
            self.popsize = popsize
        self.initialize_random_state(seed)
        self.name = 'DE'

    @staticmethod
    def initialize_random_state(seed):
        np.random.seed(seed)

    @staticmethod
    def crossover(target, mutant, probability):
        return binomial_crossover(target, mutant, probability)

    @staticmethod
    def mutate(target_idx, population, f):
        return rand1(target_idx, population, f)

    @staticmethod
    def repair(x):
        return random_repair(x)

    def init(self):
        return random_init(self.popsize, self.dims)

    def denormalize(self, population):
        return denormalize(self._MIN, self._DIFF, population)

    def mutant(self, target_idx, population, f, cr):
        # Create a mutant using a base vector
        trial = self.mutate(target_idx, population, f)
        # Repair the individual if a gene is out of bounds
        mutant = self.repair(self.crossover(population[target_idx], trial, cr))
        return mutant

    def evaluate(self, P):
        # Denormalize population matrix to obtain the scaled parameters
        PD = self.denormalize(P)
        if self.extra_params > 0:
            PD = PD[:, :-self.extra_params]
        return self.evaluate_denormalized(PD)

    def evaluate_denormalized(self, PD):
        return [self.fobj(ind) for ind in PD]

    def iterator(self):
        return iter(DEIterator(self))

    def geniterator(self):
        it = self.iterator()
        iteration = 0
        for step in it:
            if step.iteration != iteration:
                iteration = step.iteration
                yield step

    def solve(self, show_progress=False):
        best_pop_evo = []
        best_fitn_evo = []
        mean_fitn_evo = []
        if show_progress:
            from tqdm import tqdm
            iterator = tqdm(self.iterator(), total=self.maxiters, desc='Optimizing ({0})'.format(self.name))
        else:
            iterator = self.iterator()
        for step in iterator:
            idx = step.best_idx
            P = step.population

            fitness = step.fitness
            if step.iteration > self.maxiters:
                if show_progress:
                    iterator.n = self.maxiters
                    iterator.refresh()
                    iterator.close()
                break
            best_pop_evo.append(P[idx])
            best_fitn_evo.append(fitness[idx])
            mean_fitn_evo.append(np.mean(fitness))
        return self.denormalize(P[idx]), fitness[idx], \
               self.denormalize(best_pop_evo[self.popsize-1::self.popsize]), \
               np.array(best_fitn_evo[self.popsize-1::self.popsize]), \
               np.array(mean_fitn_evo[self.popsize-1::self.popsize])




class PDE(DE):
    def __init__(self, fobj, bounds, mutation=(0.5, 1.0), crossover=0.7, maxiters=1000,
                 self_adaptive=False, popsize=None, seed=None, processes=None, chunksize=None):
        super().__init__(fobj, bounds, mutation, crossover, maxiters, self_adaptive, popsize, seed)
        from multiprocessing import Pool
        self.processes = processes
        self.chunksize = chunksize
        self.name = 'Parallel DE'
        self.pool = None
        if processes is None or processes > 0:
            self.pool = Pool(processes=self.processes)

    def iterator(self):
        it = PDEIterator(self)
        try:
            for data in it:
                yield data
        finally:
            if self.pool is not None:
                self.pool.terminate()

    def evaluate_denormalized(self, PD):
        if self.pool is not None:
            return list(self.pool.map(self.fobj, PD, chunksize=self.chunksize))
        else:
            return super().evaluate_denormalized(PD)
