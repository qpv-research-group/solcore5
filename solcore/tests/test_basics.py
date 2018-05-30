from unittest import TestCase

import solcore


class TestConstants(TestCase):
    def test_01_constants_correctly_imported(self):
        q = 1.60217646e-19
        self.assertTrue(q == solcore.constants.q)


class TestMaterial(TestCase):
    def test_02_material_correctly_imported(self):
        GaAs = solcore.material('GaAs')(T=300)
        beta_Gamma_GaAs = 204
        self.assertTrue(GaAs.beta_Gamma == beta_Gamma_GaAs)


class TestUnits(TestCase):
    def test_03_units_correctly_calculated(self):
        a_nm = 1239.8417166827828
        self.assertAlmostEqual(a_nm, solcore.eVnm(1))
