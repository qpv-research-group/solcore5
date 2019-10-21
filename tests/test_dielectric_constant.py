""" Dielectric constant model tests
"""
from pytest import approx
import numpy as np
import random

from solcore.absorption_calculator.dielectric_constant_models import (
    Cauchy,
    Drude,
    Gauss,
    Lorentz,
    Poles,
    PolySegment,
)
from solcore.absorption_calculator.dielectric_constant_models import (
    DielectricConstantModel,
)
from solcore.interpolate import interp1d


def test_poles():
    A = random.uniform(0, 10)
    Ec = random.uniform(0, 10)
    pole1 = Poles(A, Ec)
    assert pole1.A == A
    assert pole1.Ec == Ec

    pole2 = Poles()
    assert pole2.dielectric(random.uniform(0.1, 10)) == 0

    pole3 = Poles(3, 2)
    assert pole3.dielectric(1240) == 4


def test_lorentz():
    An = random.uniform(0, 10)
    En = random.uniform(0, 10)
    Brn = random.uniform(0, 10)
    lorentz1 = Lorentz(An, En, Brn)
    assert lorentz1.An == An
    assert lorentz1.En == En
    assert lorentz1.Brn == Brn

    lorentz2 = Lorentz(3, 2, 3)
    assert lorentz2.dielectric(1240) == (2 + 2j)


def test_gauss():
    A = random.uniform(0, 10)
    Ec = random.uniform(0, 10)
    Br = random.uniform(0, 10)
    gauss1 = Gauss(A, Ec, Br)
    assert gauss1.A == A
    assert gauss1.Ec == Ec
    assert gauss1.Br == Br

    gauss2 = Gauss(random.uniform(0.1, 10), 0, random.uniform(0.1, 10))
    assert gauss2.dielectric(random.randint(1, 100000)) == 0

    gauss3 = Gauss(A=2.8468, Ec=0.1299, Br=0.0111)
    assert gauss3.dielectric(1240) == approx(-0.00282958252877878)

    gauss4 = Gauss(A=7.192, Ec=0.058331, Br=0.01682)
    assert gauss4.dielectric(1240) == approx(-0.004798822090878001)

    gauss5 = Gauss(A=1.9737, Ec=0.13991, Br=0.02144)
    assert gauss5.dielectric(1240) == approx(-0.004093225570635009)


def test_drude():
    An = random.uniform(0, 10)
    Brn = random.uniform(0, 10)
    drude1 = Drude(An, Brn)
    assert drude1.An == An
    assert drude1.Brn == Brn

    drude2 = Drude(random.uniform(0.1, 10), 0)
    assert drude2.dielectric(random.uniform(0.1, 10)) == 0

    drude3 = Drude(2, 1)
    assert drude3.dielectric(1240) == (-1 + 1j)


def test_cauchy():
    An = random.uniform(0, 10)
    Bn = random.uniform(0, 10)
    Cn = random.uniform(0, 10)
    Ak = random.uniform(0, 10)
    Bk = random.uniform(0, 10)
    Ck = random.uniform(0, 10)
    cauchy1 = Cauchy(An, Bn, Cn, Ak, Bk, Ck)
    assert cauchy1.An == An
    assert cauchy1.Bn == Bn
    assert cauchy1.Cn == Cn
    assert cauchy1.Ak == Ak
    assert cauchy1.Bk == Bk
    assert cauchy1.Ck == Ck

    cauchy2 = Cauchy()
    assert cauchy2.dielectric(random.randint(0, 100000)) == 0

    cauchy3 = Cauchy(1, 2, 3, 1, np.log(2), 0.24)
    assert cauchy3.dielectric(1000) == (32 + 24j)


def test_poly_segment():
    energy = np.array([1, 100])
    e2 = np.array([1, 100])
    poly1 = PolySegment(energy, e2)
    assert all([input == output for input, output in zip(energy, poly1.energy)])
    assert all([input == output for input, output in zip(e2, poly1.e2)])
    assert isinstance(poly1.epsi2, interp1d)
    assert poly1.var == 2
    interp_value = random.uniform(1, 100)
    assert poly1.epsi2(interp_value) == interp_value  # linear interpolation
    assert poly1.dielectric(1240) == approx(1.008422827635958 + 1j)


def test_dielectric_constant_model(mocker):
    mocker.patch("builtins.input", return_value="Y")
    cauchy = Cauchy()
    model = DielectricConstantModel(e_inf=0, oscillators=[cauchy])
    assert model.dielectric_constants(1000) == 0
    model.add_oscillator("drude", An=2, Brn=1)
    assert model.dielectric_constants(1240) == (-1 + 1j)
    model.add_oscillator("lorentz", An=3, En=2, Brn=3)
    assert model.dielectric_constants(1240) == (1 + 3j)
    model.remove_oscillator(2)
    assert model.dielectric_constants(1240) == (2 + 2j)
