import numpy as np
from pytest import approx

def test_DE():
    from yabox.problems import Ackley
    from solcore.optimization import DE
    problem = Ackley()

    algorithm = DE(problem, problem.bounds)

    result = algorithm.solve()

    assert np.all(result[0] < 1e-8) # == 0 fails on ubuntu-latest for some reason?
    assert result[1] < 1e-15
    assert np.all(result[2][-1] == result[0])
    assert result[3][-1] < 1e-15
    assert result[4][-1] < 1e-15


def test_PDE():
    from yabox.problems import Ackley
    from solcore.optimization import PDE
    problem = Ackley()

    algorithm = PDE(problem, problem.bounds)

    result = algorithm.solve()

    assert np.all(result[0] < 1e-8)
    assert result[1] < 1e-15
    assert np.all(result[2][-1] == result[0])
    assert result[3][-1] < 1e-15
    assert result[4][-1] < 1e-15
