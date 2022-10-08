from pytest import approx, mark


@mark.flaky(reruns=5)
def test_DE():
    from yabox.problems import Ackley
    from solcore.optimization import DE

    problem = Ackley()

    algorithm = DE(problem, problem.bounds)

    result = algorithm.solve()

    assert result[0] == approx(0.0)
    assert result[1] == approx(0.0)
    assert result[2][-1] == approx(result[0])
    assert result[3][-1] == approx(0.0)
    assert result[4][-1] == approx(0.0)


@mark.flaky(reruns=5)
def test_PDE():
    from yabox.problems import Ackley
    from solcore.optimization import PDE

    problem = Ackley()

    algorithm = PDE(problem, problem.bounds, processes=0)

    result = algorithm.solve()

    assert result[0] == approx(0.0)
    assert result[1] == approx(0.0)
    assert result[2][-1] == approx(result[0])
    assert result[3][-1] == approx(0.0)
    assert result[4][-1] == approx(0.0)
