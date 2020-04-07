import numpy as np

def test_DE():
    from yabox.problems import Ackley
    from solcore.optimization import DE
    problem = Ackley()

    algorithm = DE(problem, problem.bounds)

    result = algorithm.solve()

    assert np.all(result[0] == np.array([0, 0]))
    assert result[1] < 1e-15
    assert np.all(result[2][-1] == np.array([0, 0]))
    assert result[3][-1] < 1e-15
    assert result[4][-1] < 1e-15


def test_PDE():
    from yabox.problems import Ackley
    from solcore.optimization import PDE
    problem = Ackley()

    algorithm = PDE(problem, problem.bounds)

    result = algorithm.solve()

    assert np.all(result[0] == np.array([0, 0]))
    assert result[1] < 1e-15
    assert np.all(result[2][-1] == np.array([0, 0]))
    assert result[3][-1] < 1e-15
    assert result[4][-1] < 1e-15