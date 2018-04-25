"""Quantum mechanics tests"""
from unittest import TestCase

import solcore
from solcore import material, si
from solcore.constants import vacuum_permittivity, q
from solcore.structure import Layer, Structure
import solcore.quantum_mechanics as QM
import numpy as np

# Energies in meV
my_energies = {'Elh': np.array([-770., -802., -803., -810., -812., -818., -823., -827.]),
               'Ee': np.array([577., 629., 630., 640., 643., 648., 658., 663.]),
               'Ehh': np.array([-677., -731., -801., -801., -803., -803., -806., -807., -811., -813., -815., -816.,
                                -820., -824., -828., -832., -838., -841., -845., -849., ])}

my_absorption = [1.267056856187291, 988133.9170541381]  # Energy in meV and absorption coeficent in m-1


class TestQM(TestCase):
    def test_61_kp_bulk_solvers(self):
        GaAs = solcore.material("GaAs")(T=300)
        GaAsP = solcore.material("GaAsP")(P=0.3, T=300)
        InGaAs = solcore.material("InGaAs")(In=0.2, T=300)

        # Testing QM.kp_bands
        bands_GaAsP = QM.kp_bands(GaAsP, GaAs, graph=False, fit_effective_mass=True, effective_mass_direction="X",
                                  return_so=True)
        bands_InGaAs = QM.kp_bands(InGaAs, GaAs, graph=False, fit_effective_mass=True, effective_mass_direction="X",
                                   return_so=True)

        expected_bands_GaAsP = (1.2168382480631407e-19, -1.5452519153253004e-19, -1.4042149045828435e-19,
                                -1.9192138182935611e-19, 8.0093555784846597e-32, 1.2472835955929216e-31,
                                2.6423749777877535e-31, 1.2393634061184521e-31)
        expected_bands_InGaAs = (8.6764014773233634e-20, -1.0573103504669798e-19, -1.1984351916698905e-19,
                                 -1.6993543036257329e-19, 6.6193386922591731e-32, 1.3576713980579555e-31,
                                 8.0904387208083259e-32, 1.2268187169919973e-31)

        for i in range(len(bands_GaAsP)):
            self.assertAlmostEqual(bands_GaAsP[i], expected_bands_GaAsP[i])
            self.assertAlmostEqual(bands_InGaAs[i], expected_bands_InGaAs[i])

        # Testing QM.KPbands and QM.fit_effective_masses
        edges_GaAsP = QM.KPbands(GaAsP, GaAs, fraction=0.2, return_edges_only=True)
        bands_GaAsP = QM.KPbands(GaAsP, GaAs, fraction=0.2)
        masses_GaAsP = QM.fit_effective_masses(bands_GaAsP, GaAsP, GaAs, plot_result=False)

        edges_InGaAs = QM.KPbands(InGaAs, GaAs, fraction=0.2, return_edges_only=True)
        bands_InGaAs = QM.KPbands(InGaAs, GaAs, fraction=0.2)
        masses_InGaAs = QM.fit_effective_masses(bands_InGaAs, InGaAs, GaAs, plot_result=False)

        expected_edges_GaAsP = (1.2168382480631407e-19, -1.5452519153253004e-19, -1.4042149045828435e-19,
                                -1.9192138182935611e-19)

        expected_edges_InGaAs = (8.6764014773233634e-20, -1.0573103504669798e-19, -1.1984351916698905e-19,
                                 -1.6993543036257329e-19)

        expected_masses_GaAsP = (8.049577422084102e-32, 1.2627430248682043e-31, 2.6577242586172804e-31,
                                 1.2305748108835472e-31)

        expected_masses_InGaAs = (6.6895885457875e-32, 1.3994390560400583e-31, 8.142667105522975e-32,
                                  1.2060355194525871e-31)

        for i in range(len(edges_GaAsP)):
            self.assertAlmostEqual(edges_GaAsP[i], expected_edges_GaAsP[i])
            self.assertAlmostEqual(edges_InGaAs[i], expected_edges_InGaAs[i])
            self.assertAlmostEqual(masses_GaAsP[i], expected_masses_GaAsP[i])
            self.assertAlmostEqual(masses_InGaAs[i], expected_masses_InGaAs[i])

        # Testing QM.kp8x8_bulk
        bands_GaAsP = QM.kp8x8_bulk(GaAsP, GaAs)
        bands_InGaAs = QM.kp8x8_bulk(InGaAs, GaAs)

        expected_bands_GaAsP = (1.2168382480631407e-19, -1.5452519153253004e-19, -1.4042149045828435e-19,
                                -1.9192138182935611e-19, 8.138445281947437e-32, 1.3543674202428507e-31,
                                2.848952594319033e-31, 1.145125159048442e-31)

        expected_bands_InGaAs = (8.6764014773233634e-20, -1.0573103504669798e-19, -1.1984351916698905e-19,
                                 -1.6993543036257329e-19, 6.77957483053393e-32, 1.6988821114765817e-31,
                                 7.820551493038613e-32, 1.1461138067300424e-31)

        for i in range(len(bands_GaAsP)):
            self.assertAlmostEqual(bands_GaAsP[i], expected_bands_GaAsP[i])
            self.assertAlmostEqual(bands_InGaAs[i], expected_bands_InGaAs[i])

    def test_64_quantum_mechanics_schrodinger(self):
        bulk = material("GaAs")(T=293)
        barrier = material("GaAsP")(T=293, P=0.1)

        bulk.strained = False
        barrier.strained = True

        top_layer = Layer(width=si("30nm"), material=bulk)
        inter = Layer(width=si("3nm"), material=bulk)
        barrier_layer = Layer(width=si("15nm"), material=barrier)
        bottom_layer = top_layer

        E = np.linspace(1.15, 1.5, 300) * q
        alfas = np.zeros((len(E), 6))

        alfas[:, 0] = E / q

        alpha_params = {
            "well_width": si("7.2nm"),
            "theta": 0,
            "eps": 12.9 * vacuum_permittivity,
            "espace": E,
            "hwhm": si("6meV"),
            "dimensionality": 0.16,
            "line_shape": "Gauss"
        }

        QW = material("InGaAs")(T=293, In=0.2)
        QW.strained = True
        well_layer = Layer(width=si("7.2nm"), material=QW)

        my_structure = Structure(
            [top_layer, barrier_layer, inter, well_layer, inter, barrier_layer, inter, bottom_layer], substrate=bulk)

        band_edge, bands = QM.schrodinger(my_structure, quasiconfined=0,
                                          num_eigenvalues=20, alpha_params=alpha_params, calculate_absorption=True)

        for key in band_edge['E']:
            for i in range(len(band_edge['E'][key])):
                band_edge['E'][key][i] = solcore.asUnit(band_edge['E'][key][i], 'eV') * 1000
                band_edge['E'][key][i] = round(band_edge['E'][key][i])

        Ehh = np.all(np.equal(band_edge['E']['Ehh'], my_energies['Ehh']))
        Elh = np.all(np.equal(band_edge['E']['Elh'], my_energies['Elh']))
        Ee = np.all(np.equal(band_edge['E']['Ee'], my_energies['Ee']))

        idx = 100
        out = [band_edge['alpha'][0][idx] / q, band_edge['alpha'][1][idx]]

        # Test over the energies
        self.assertTrue(Ehh and Elh and Ee)
        # Test over the absorption coefficent at a given energy
        for i, data in enumerate(out):
            self.assertAlmostEqual(out[i], my_absorption[i])
