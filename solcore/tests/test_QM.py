from unittest import TestCase

import solcore
from solcore import material, si
from solcore.constants import vacuum_permittivity, q
from solcore.structure import Layer
import solcore.quantum_mechanics as QM
import solcore.quantum_mechanics.structure_utilities as su
import numpy as np

# Energies in meV
my_energies = {'Elh': np.array([-755., -806., -836.]), 'Ee': np.array([551., 616., 672.]),
               'Ehh': np.array([-644., -675., -717., -752., -785., -833., -860., -861.])}

my_absorption = [1.3031847133757961, 350268.06722092471]  # Energy in meV and absorption coefficent in m-1


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
        GaAs = material("GaAs")
        InGaAs = material("InGaAs")
        GaAsP = material("GaAsP")
        GaP = material("GaP")

        InGaAs.strained = True
        GaAsP.strained = True
        GaP.strained = True
        GaAs.strained = False

        bulkMaterial = GaAs(T=300)
        wellMaterial = InGaAs(In=0.245, T=300)
        barrierMaterial = GaAsP(P=0.1, T=300)
        interlayer_material = InGaAs(In=0.14, T=300)

        my_structure = QM.structure_utilities.assemble_qw_structure(
            repeats=1,
            well=Layer(si("7.2nm"), wellMaterial),
            bulk_l_top=Layer(si("1nm"), bulkMaterial),
            bulk_l_bottom=Layer(si("1nm"), bulkMaterial),
            barrier=Layer(si("28nm"), barrierMaterial),
            well_interlayer=Layer(si("3nm"), interlayer_material)
        )
        my_structure.substrate = bulkMaterial

        wl = np.linspace(800, 1300, 100)
        E = 1240 / wl * q

        alpha_params = {
            "well_width": my_structure.width(),
            "theta": 0,
            "eps": 12.9 * vacuum_permittivity,
            "espace": E,
            "hwhm": si("4meV"),
            "dimensionality": 0.2,
        }

        band_edge, bands = QM.schrodinger(my_structure, kpoints=30, krange=2e9, num_eigenvalues=20, plot_bands=False,
                                          symmetric=False, quasiconfined=0, blur=None, blurmode="right", Efield=0,
                                          mode='strain', smallest_feature_steps=10, filter_strength=0,
                                          calculate_absorption=True, alpha_params=alpha_params)

        for key in band_edge['E']:
            band_edge['E'][key] = solcore.asUnit(band_edge['E'][key], 'eV') * 1000
            for i in range(len(band_edge['E'][key])):
                band_edge['E'][key][i] = round(band_edge['E'][key][i])

        Ehh = np.all(np.equal(band_edge['E']['Ehh'], my_energies['Ehh']))
        Elh = np.all(np.equal(band_edge['E']['Elh'], my_energies['Elh']))
        Ee = np.all(np.equal(band_edge['E']['Ee'], my_energies['Ee']))

        idx = 30
        out = [band_edge['alpha'][0][idx] / q, band_edge['alpha'][1][idx]]

        # Test over the energies
        self.assertTrue(Ehh and Elh and Ee)
        # Test over the absorption coefficent at a given energy
        for i, data in enumerate(out):
            self.assertAlmostEqual(out[i], my_absorption[i])
