from pytest import approx

import solcore


def test_constants_correctly_imported():
    q = 1.60217646e-19
    assert q == solcore.constants.q


def test_material_correctly_imported():
    GaAs = solcore.material("GaAs")(T=300)
    beta_Gamma_GaAs = 204
    assert GaAs.beta_Gamma == beta_Gamma_GaAs


def test_units_correctly_calculated():
    a_nm = 1239.8417166827828
    assert a_nm == approx(solcore.eVnm(1))
