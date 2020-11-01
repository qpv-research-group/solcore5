"""
This is completely based on Pablo R. Mier's yabox package, with modifications to slightly change what is returned by
DE.solve to be able to visualize the convergence of the DE without having to use de.iterator, and just return
a list of optimized parameters rather than a reshaped array.

All credit of the algorithm, testing, etc. goes to Pablo R. Mier. For more details, visit:
- https://github.com/pablormier/yabox
- https://pablormier.github.io/2017/09/05/a-tutorial-on-differential-evolution-with-python/
"""

from yabox.algorithms.de import PDEIterator
from yabox.algorithms.de import DE as DE_yabox
import numpy as np


class DE(DE_yabox):
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
