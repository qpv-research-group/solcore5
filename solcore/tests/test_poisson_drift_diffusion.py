from unittest import TestCase

import solcore
import solcore.poisson_drift_diffusion as PDD
from solcore.structure import Layer
from solcore.light_source import LightSource
from solcore import material
import tempfile
import os
import numpy as np

# We create the other materials we need for the device
window = material('AlGaAs')(T=300, Na=1e24, Al=0.8)
p_GaAs = material('GaAs')(T=300, Na=1e24)
i_GaAs = material('GaAs')(T=300)
n_GaAs = material('GaAs')(T=300, Nd=1e23)
bsf = material('AlGaAs')(T=300, Nd=1e24, Al=0.4)

# And finally we create another p-i-n structure incorporating the QWs in the intrinsic region.
MyDevice = PDD.CreateDeviceStructure('TestDevice', T=300, layers=[
    Layer(width=30e-9, material=window, role="Window"),
    Layer(width=400e-9, material=p_GaAs, role="Emitter"),
    Layer(width=500e-9, material=i_GaAs, role="Intrinsic"),
    Layer(width=2000e-9, material=n_GaAs, role="Base"),
    Layer(width=200e-9, material=bsf, role="BSF")
])

my_answer = {'layers': [{'class': 'Layer', 'repeat': 1, 'numlayers': 1,
                         'properties': {'sp': 1000000.0, 'electron_auger_recombination': 1e-42,
                                        'Nv': 9.0651047816548187e+24, 'band_gap': 3.2412807434271126e-19,
                                        'sn': 1000000.0, 'hole_minority_lifetime': 2.5e-07, 'Na': 1e+24,
                                        'Nc': 1.2226626898157938e+24, 'electron_mobility': 0.012113993200582249,
                                        'radiative_recombination': 7.2e-16, 'permittivity': 12.9,
                                        'electron_affinity': 6.3461431932128865e-19,
                                        'electron_minority_lifetime': 3e-06, 'eff_mass_hh_z': 0.44812181631113623,
                                        'eff_mass_lh_z': 0.15540371244551329,
                                        'composition': {'element': 'Al', 'fraction': 0.8, 'material': 'AlGaAs'},
                                        'eff_mass_electron_Gamma': 0.1334, 'width': 3e-08, 'Nd': 1,
                                        'hole_auger_recombination': 1e-42, 'ni': 33845254.328695573,
                                        'hole_mobility': 0.0047195310636572903}, 'label': 'Window', 'group': None},
                        {'class': 'Layer', 'repeat': 1, 'numlayers': 1,
                         'properties': {'sp': 1000000.0, 'electron_auger_recombination': 1e-42,
                                        'Nv': 5.6577936548188829e+24, 'band_gap': 2.2790674040560711e-19,
                                        'sn': 1000000.0, 'hole_minority_lifetime': 2.5e-07, 'Na': 1e+24,
                                        'Nc': 4.3519622483962564e+23, 'electron_mobility': 0.27085675648546026,
                                        'radiative_recombination': 7.2e-16, 'permittivity': 12.9,
                                        'electron_affinity': 6.62903371354393e-19, 'electron_minority_lifetime': 3e-06,
                                        'eff_mass_hh_z': 0.34129421032745344, 'eff_mass_lh_z': 0.087937066805091599,
                                        'composition': {'material': 'GaAs'}, 'eff_mass_electron_Gamma': 0.067,
                                        'width': 4e-07, 'Nd': 1, 'hole_auger_recombination': 1e-42,
                                        'ni': 1767480124457.7329, 'hole_mobility': 0.017374281562790292},
                         'label': 'Emitter', 'group': None}, {'class': 'Layer', 'repeat': 1, 'numlayers': 1,
                                                              'properties': {'sp': 1000000.0,
                                                                             'electron_auger_recombination': 1e-42,
                                                                             'Nv': 5.6577936548188829e+24,
                                                                             'band_gap': 2.2790674040560711e-19,
                                                                             'sn': 1000000.0,
                                                                             'hole_minority_lifetime': 2.5e-07, 'Na': 1,
                                                                             'Nc': 4.3519622483962564e+23,
                                                                             'electron_mobility': 0.9399999987600498,
                                                                             'radiative_recombination': 7.2e-16,
                                                                             'permittivity': 12.9,
                                                                             'electron_affinity': 6.62903371354393e-19,
                                                                             'electron_minority_lifetime': 3e-06,
                                                                             'eff_mass_hh_z': 0.34129421032745344,
                                                                             'eff_mass_lh_z': 0.087937066805091599,
                                                                             'composition': {'material': 'GaAs'},
                                                                             'eff_mass_electron_Gamma': 0.067,
                                                                             'width': 5e-07, 'Nd': 1,
                                                                             'hole_auger_recombination': 1e-42,
                                                                             'ni': 1767480124457.7329,
                                                                             'hole_mobility': 0.04914999990380032},
                                                              'label': 'Intrinsic', 'group': None},
                        {'class': 'Layer', 'repeat': 1, 'numlayers': 1,
                         'properties': {'sp': 1000000.0, 'electron_auger_recombination': 1e-42,
                                        'Nv': 5.6577936548188829e+24, 'band_gap': 2.2790674040560711e-19,
                                        'sn': 1000000.0, 'hole_minority_lifetime': 2.5e-07, 'Na': 1,
                                        'Nc': 4.3519622483962564e+23, 'electron_mobility': 0.4503690283162032,
                                        'radiative_recombination': 7.2e-16, 'permittivity': 12.9,
                                        'electron_affinity': 6.62903371354393e-19, 'electron_minority_lifetime': 3e-06,
                                        'eff_mass_hh_z': 0.34129421032745344, 'eff_mass_lh_z': 0.087937066805091599,
                                        'composition': {'material': 'GaAs'}, 'eff_mass_electron_Gamma': 0.067,
                                        'width': 2e-06, 'Nd': 1e+23, 'hole_auger_recombination': 1e-42,
                                        'ni': 1767480124457.7329, 'hole_mobility': 0.027327813913248372},
                         'label': 'Base', 'group': None}, {'class': 'Layer', 'repeat': 1, 'numlayers': 1,
                                                           'properties': {'sp': 1000000.0,
                                                                          'electron_auger_recombination': 1e-42,
                                                                          'Nv': 7.9627295597221238e+24,
                                                                          'band_gap': 3.0008009508278419e-19,
                                                                          'sn': 1000000.0,
                                                                          'hole_minority_lifetime': 2.5e-07, 'Na': 1,
                                                                          'Nc': 7.959293095895547e+23,
                                                                          'electron_mobility': 0.013014278277956979,
                                                                          'radiative_recombination': 7.2e-16,
                                                                          'permittivity': 12.9,
                                                                          'electron_affinity': 6.2469615762921586e-19,
                                                                          'electron_minority_lifetime': 3e-06,
                                                                          'eff_mass_hh_z': 0.42218491524687168,
                                                                          'eff_mass_lh_z': 0.12273039939829694,
                                                                          'composition': {'element': 'Al',
                                                                                          'fraction': 0.4,
                                                                                          'material': 'AlGaAs'},
                                                                          'eff_mass_electron_Gamma': 0.1002,
                                                                          'width': 2e-07, 'Nd': 1e+24,
                                                                          'hole_auger_recombination': 1e-42,
                                                                          'ni': 466523271.02800423,
                                                                          'hole_mobility': 0.0061990364205366678},
                                                           'label': 'BSF', 'group': None}], 'reflection': None,
             'repeat': 1, 'numlayers': 5, 'substrate': {'material': 'GaAs'}, 'role': 'device', 'T': 300,
             'name': 'TestDevice', 'comments': ''}


class TestPDD(TestCase):
    def test_91_device_construccion(self):
        self.assertTrue(MyDevice == my_answer)

    def test_92_light_iv(self):
        with tempfile.TemporaryDirectory(prefix="tmp", suffix="_sc3TESTS") as working_directory:
            filename = os.path.join(working_directory, 'solcore_log.txt')
            PDD.log(filename)
            IV = PDD._IV(MyDevice, vfin=1.2, vstep=0.05, light=True)
            Jsc = 194.25993896711984
        self.assertAlmostEqual(Jsc, IV['IV']['Jsc'])

    def test_93_qe(self):
        with tempfile.TemporaryDirectory(prefix="tmp", suffix="_sc3TESTS") as working_directory:
            filename = os.path.join(working_directory, 'solcore_log.txt')
            PDD.log(filename)
            QE = PDD._QE(MyDevice)
            qe_wl = 0.94589342955180866
        self.assertAlmostEqual(qe_wl, QE['QE']['IQE'][130])

    def test_94_tmm_optics_calculator(self):
        with tempfile.TemporaryDirectory(prefix="tmp", suffix="_sc3TESTS") as working_directory:
            filename = os.path.join(working_directory, 'solcore_log.txt')
            PDD.log(filename)
            wavelength = np.linspace(300, 1000, 200) * 1e-9
            output = PDD.calculate_optics(MyDevice, wavelength)
            expected = 9833723.1702278145
        self.assertAlmostEqual(expected, output['absorption'][10, 10])
